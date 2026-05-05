"""Export an external perturbation-response dictionary for DevVCell.

The expected primary input is an H5AD converted from a public perturbation
dataset. If the H5AD is not present yet, the script writes a metadata manifest
from ``config/external_datasets.json`` so the acquisition state is explicit.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from devvcell.io import load_json, resolve_project_path  # noqa: E402
from devvcell.tables import write_table  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/perturbation_transfer.json")
    parser.add_argument("--allow-missing", action="store_true", help="Write dataset manifest instead of failing when H5AD is absent.")
    return parser.parse_args()


def as_csr_float32(matrix) -> sparse.csr_matrix:
    if sparse.issparse(matrix):
        return matrix.tocsr().astype(np.float32)
    return sparse.csr_matrix(np.asarray(matrix, dtype=np.float32))


def selected_indices(obs: pd.DataFrame, cfg: dict, rng: np.random.Generator) -> np.ndarray:
    schema = cfg["external_schema"]
    max_cells = int(cfg["encoding"]["max_cells"])
    per_condition = int(cfg["encoding"]["max_cells_per_condition"])
    condition_col = schema["condition_column"]
    perturb_col = schema["perturbation_column"]
    if condition_col not in obs.columns:
        condition_col = perturb_col
    if perturb_col not in obs.columns:
        raise ValueError(f"External H5AD obs is missing perturbation column: {perturb_col}")

    sampled: list[np.ndarray] = []
    for _, group in obs.groupby([condition_col, perturb_col], observed=True, sort=False):
        take = min(per_condition, len(group))
        sampled.append(rng.choice(group.index.to_numpy(), size=take, replace=False))
    selected = np.concatenate(sampled) if sampled else np.array([], dtype=int)
    if len(selected) > max_cells:
        selected = rng.choice(selected, size=max_cells, replace=False)
    return np.sort(selected.astype(int))


def fit_latent(X: sparse.csr_matrix, latent_dim: int, seed: int) -> tuple[np.ndarray, TruncatedSVD, StandardScaler]:
    n_components = min(int(latent_dim), X.shape[0] - 1, X.shape[1] - 1)
    if n_components < 2:
        raise RuntimeError("Need at least 3 cells and 3 genes to fit external latent responses.")
    svd = TruncatedSVD(n_components=n_components, random_state=seed)
    raw = svd.fit_transform(X)
    scaler = StandardScaler()
    latent = scaler.fit_transform(raw).astype(np.float32)
    return latent, svd, scaler


def write_missing_manifest(cfg: dict) -> None:
    registry = load_json(cfg["input"]["external_dataset_registry"])
    rows: list[dict[str, object]] = []
    for section in ["primary_candidates", "normal_manifold_validation", "gonad_niche_validation"]:
        for item in registry.get(section, []):
            base = {
                "registry_section": section,
                "dataset_id": item.get("dataset_id"),
                "title": item.get("title"),
                "organism": item.get("organism"),
                "role": item.get("role"),
                "priority": item.get("priority"),
                "status": item.get("status", "planned"),
                "landing_page": item.get("landing_page"),
            }
            files = item.get("processed_files") or [{}]
            for file_item in files:
                rows.append({**base, **{f"file_{k}": v for k, v in file_item.items()}})
    output = cfg["output"]["external_response_metadata"]
    path = write_table(pd.DataFrame(rows), output)
    print(f"External H5AD is absent; wrote acquisition manifest: {path}")


def main() -> None:
    args = parse_args()
    cfg = load_json(args.config)
    rng = np.random.default_rng(int(cfg["seed"]))
    h5ad_path = resolve_project_path(cfg["input"]["external_h5ad"])
    if not h5ad_path.exists():
        if args.allow_missing:
            write_missing_manifest(cfg)
            return
        raise FileNotFoundError(f"External perturbation H5AD is missing: {h5ad_path}")

    adata = ad.read_h5ad(h5ad_path, backed="r")
    obs = adata.obs.reset_index(drop=True)
    schema = cfg["external_schema"]
    selected = selected_indices(obs, cfg, rng)
    sample_obs = obs.iloc[selected].copy().reset_index(drop=True)
    X = as_csr_float32(adata.X[selected, :])
    latent, svd, _ = fit_latent(X, int(cfg["encoding"]["latent_dim"]), int(cfg["seed"]))
    latent_cols = [f"latent_{i + 1:02d}" for i in range(latent.shape[1])]
    latent_df = pd.DataFrame(latent, columns=latent_cols)
    sample = pd.concat([sample_obs, latent_df], axis=1)

    perturb_col = schema["perturbation_column"]
    condition_col = schema["condition_column"] if schema["condition_column"] in sample.columns else perturb_col
    cell_type_col = schema["cell_type_column"] if schema["cell_type_column"] in sample.columns else None
    control_label = str(schema["control_label"])

    grouping = [condition_col]
    if cell_type_col is not None:
        grouping.append(cell_type_col)
    centroids = sample.groupby(grouping, as_index=False, observed=True)[latent_cols].mean()
    counts = sample.groupby(grouping, as_index=False, observed=True).size().rename(columns={"size": "n_cells"})
    centroids = centroids.merge(counts, on=grouping, how="left")

    control_mask = centroids[condition_col].astype(str).str.lower().eq(control_label.lower())
    if not control_mask.any():
        control_mask = centroids[condition_col].astype(str).str.contains("control|ntc|untreated", case=False, regex=True)
    if not control_mask.any():
        raise RuntimeError("Could not identify a control condition in external perturbation H5AD.")

    control_centroids = centroids[control_mask].copy()
    perturb_centroids = centroids[~control_mask].copy()
    rows: list[dict[str, object]] = []
    for _, pert in perturb_centroids.iterrows():
        if cell_type_col is not None:
            ctrl_pool = control_centroids[control_centroids[cell_type_col].astype(str) == str(pert[cell_type_col])]
            if ctrl_pool.empty:
                ctrl_pool = control_centroids
        else:
            ctrl_pool = control_centroids
        ctrl = ctrl_pool.sort_values("n_cells", ascending=False).iloc[0]
        response = pert[latent_cols].astype(float).to_numpy() - ctrl[latent_cols].astype(float).to_numpy()
        row = {
            "perturbation": str(pert[condition_col]),
            "control_condition": str(ctrl[condition_col]),
            "external_cell_type": str(pert[cell_type_col]) if cell_type_col is not None else "all",
            "n_perturbed_cells": int(pert["n_cells"]),
            "n_control_cells": int(ctrl["n_cells"]),
            "response_norm": float(np.linalg.norm(response)),
            "svd_explained_variance_ratio_sum": float(np.sum(svd.explained_variance_ratio_)),
        }
        row.update({f"response_latent_{i + 1:02d}": float(value) for i, value in enumerate(response)})
        row.update({f"control_latent_{i + 1:02d}": float(ctrl[col]) for i, col in enumerate(latent_cols)})
        rows.append(row)

    dictionary = pd.DataFrame(rows).sort_values(["perturbation", "external_cell_type"])
    metadata = pd.DataFrame(
        [
            {
                "external_h5ad": str(h5ad_path),
                "sampled_cells": int(len(sample)),
                "n_responses": int(len(dictionary)),
                "latent_dim": int(len(latent_cols)),
                "svd_explained_variance_ratio_sum": float(np.sum(svd.explained_variance_ratio_)),
            }
        ]
    )
    dict_path = write_table(dictionary, cfg["output"]["external_response_dictionary"])
    meta_path = write_table(metadata, cfg["output"]["external_response_metadata"])
    print(f"Wrote response dictionary: {dict_path}")
    print(f"Wrote response metadata: {meta_path}")


if __name__ == "__main__":
    main()
