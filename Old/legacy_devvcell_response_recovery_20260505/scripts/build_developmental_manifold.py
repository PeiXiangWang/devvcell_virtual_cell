"""Build the embryo developmental manifold used by DevVCell response-recovery.

The script is designed for the 11M-cell embryo H5AD in this repository. It
opens the file in backed mode, samples cells by stage and cell type, learns a
TruncatedSVD latent space, writes sampled cell embeddings and stage/cell-type
centroids, and reports classifier-quality metrics.
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
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from devvcell.io import load_json, resolve_project_path  # noqa: E402
from devvcell.stages import canonical_stage, stage_number  # noqa: E402
from devvcell.tables import write_table  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/response_recovery.json")
    parser.add_argument("--max-cells", type=int, default=None)
    return parser.parse_args()


def as_csr_float32(matrix) -> sparse.csr_matrix:
    if sparse.issparse(matrix):
        return matrix.tocsr().astype(np.float32)
    return sparse.csr_matrix(np.asarray(matrix, dtype=np.float32))


def assign_system(obs: pd.DataFrame, config: dict) -> pd.Series:
    cell_type_col = config["input"]["cell_type_column"]
    major_col = config["input"]["major_cluster_column"]
    text = (
        obs[cell_type_col].astype(str).str.lower().fillna("")
        + " "
        + obs[major_col].astype(str).str.lower().fillna("")
    )
    labels = pd.Series("unassigned", index=obs.index, dtype=object)
    for system_name, keywords in config["manifold"]["systems"].items():
        mask = np.zeros(len(obs), dtype=bool)
        for keyword in keywords:
            mask |= text.str.contains(str(keyword).lower(), regex=False).to_numpy()
        labels.loc[mask & (labels == "unassigned")] = system_name
    return labels


def stratified_sample(obs: pd.DataFrame, config: dict, rng: np.random.Generator, max_cells: int | None) -> np.ndarray:
    stage_col = config["input"]["stage_column"]
    cell_type_col = config["input"]["cell_type_column"]
    per_group = int(config["manifold"]["max_cells_per_stage_celltype"])
    total_limit = int(max_cells or config["manifold"]["max_cells_total"])

    obs = obs.copy()
    obs["devvcell_system"] = assign_system(obs, config)
    obs = obs[obs["devvcell_system"] != "unassigned"]
    if obs.empty:
        raise RuntimeError("No cells matched configured DevVCell systems.")

    sampled: list[np.ndarray] = []
    group_cols = ["devvcell_system", stage_col, cell_type_col]
    for _, group in obs.groupby(group_cols, observed=True, sort=False):
        take = min(per_group, len(group))
        if take <= 0:
            continue
        sampled.append(rng.choice(group.index.to_numpy(), size=take, replace=False))

    selected = np.concatenate(sampled) if sampled else np.array([], dtype=int)
    if len(selected) > total_limit:
        selected = rng.choice(selected, size=total_limit, replace=False)
    return np.sort(selected.astype(int))


def fit_latent(X: sparse.csr_matrix, latent_dim: int, seed: int) -> tuple[np.ndarray, TruncatedSVD, StandardScaler]:
    n_components = min(int(latent_dim), X.shape[1] - 1, X.shape[0] - 1)
    if n_components < 2:
        raise RuntimeError("Need at least 3 cells and 3 genes to fit a latent manifold.")
    svd = TruncatedSVD(n_components=n_components, random_state=seed)
    raw = svd.fit_transform(X)
    scaler = StandardScaler()
    latent = scaler.fit_transform(raw).astype(np.float32)
    return latent, svd, scaler


def classifier_accuracy(latent: np.ndarray, labels: pd.Series, test_fraction: float, seed: int) -> float:
    valid = labels.astype(str)
    counts = valid.value_counts()
    keep = valid.isin(counts[counts >= 3].index)
    if keep.sum() < 20 or valid[keep].nunique() < 2:
        return float("nan")
    x_train, x_test, y_train, y_test = train_test_split(
        latent[keep.to_numpy()],
        valid[keep].to_numpy(),
        test_size=test_fraction,
        random_state=seed,
        stratify=valid[keep].to_numpy(),
    )
    model = LogisticRegression(max_iter=800)
    model.fit(x_train, y_train)
    return float(accuracy_score(y_test, model.predict(x_test)))


def main() -> None:
    args = parse_args()
    config = load_json(args.config)
    rng = np.random.default_rng(int(config["seed"]))

    h5ad_path = resolve_project_path(config["input"]["embryo_h5ad"])
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Embryo H5AD is missing: {h5ad_path}")

    adata = ad.read_h5ad(h5ad_path, backed="r")
    obs = adata.obs.reset_index(drop=True)
    selected = stratified_sample(obs, config, rng, args.max_cells)
    X = as_csr_float32(adata.X[selected, :])
    sample_obs = obs.iloc[selected].copy().reset_index(drop=True)
    sample_obs["source_cell_index"] = selected
    sample_obs["stage"] = sample_obs[config["input"]["stage_column"]].map(canonical_stage)
    sample_obs["stage_num"] = sample_obs["stage"].map(stage_number)
    sample_obs["cell_type"] = sample_obs[config["input"]["cell_type_column"]].astype(str)
    sample_obs["devvcell_system"] = assign_system(sample_obs, config).to_numpy()

    latent, svd, _ = fit_latent(X, int(config["manifold"]["latent_dim"]), int(config["seed"]))
    latent_cols = [f"latent_{i + 1:02d}" for i in range(latent.shape[1])]
    latent_df = pd.DataFrame(latent, columns=latent_cols)

    keep_cols = [
        "source_cell_index",
        "stage",
        "stage_num",
        "cell_type",
        "devvcell_system",
        config["input"]["donor_column"],
        config["input"]["sex_column"],
        config["input"]["major_cluster_column"],
    ]
    keep_cols = [col for col in keep_cols if col in sample_obs.columns]
    cells = pd.concat([sample_obs[keep_cols], latent_df], axis=1)

    min_cells = int(config["manifold"]["min_cells_per_centroid"])
    centroids = (
        cells.groupby(["devvcell_system", "stage", "stage_num", "cell_type"], as_index=False, observed=True)
        .agg(n_cells=("source_cell_index", "size"), **{col: (col, "mean") for col in latent_cols})
    )
    centroids = centroids[centroids["n_cells"] >= min_cells].reset_index(drop=True)

    quality = pd.DataFrame(
        [
            {
                "metric": "sampled_cells",
                "value": float(len(cells)),
                "description": "Cells sampled into the DevVCell developmental manifold.",
            },
            {
                "metric": "centroids",
                "value": float(len(centroids)),
                "description": "Stage/cell-type centroids with enough sampled cells.",
            },
            {
                "metric": "svd_explained_variance_ratio_sum",
                "value": float(np.sum(svd.explained_variance_ratio_)),
                "description": "Variance captured by the latent SVD components.",
            },
            {
                "metric": "stage_classifier_accuracy",
                "value": classifier_accuracy(
                    latent,
                    cells["stage"],
                    float(config["manifold"]["classifier_test_fraction"]),
                    int(config["seed"]),
                ),
                "description": "Held-out logistic-regression accuracy for Theiler stage labels.",
            },
            {
                "metric": "cell_type_classifier_accuracy",
                "value": classifier_accuracy(
                    latent,
                    cells["cell_type"],
                    float(config["manifold"]["classifier_test_fraction"]),
                    int(config["seed"]),
                ),
                "description": "Held-out logistic-regression accuracy for cell type labels.",
            },
        ]
    )

    paths = config["output"]
    cells_path = write_table(cells, paths["embryo_manifold_cells"])
    centroids_path = write_table(centroids, paths["stage_celltype_centroids"])
    metrics_path = write_table(quality, paths["manifold_quality_metrics"])
    print(f"Wrote manifold cells: {cells_path}")
    print(f"Wrote centroids: {centroids_path}")
    print(f"Wrote quality metrics: {metrics_path}")


if __name__ == "__main__":
    main()
