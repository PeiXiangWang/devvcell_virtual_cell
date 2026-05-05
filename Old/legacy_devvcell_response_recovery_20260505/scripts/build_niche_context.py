"""Build niche-aware context features and the Tbx4-Glis3 case study."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from devvcell.io import load_json, resolve_project_path  # noqa: E402
from devvcell.stages import canonical_stage, stage_number  # noqa: E402
from devvcell.tables import write_table  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/niche_context.json")
    parser.add_argument("--max-cells", type=int, default=None)
    return parser.parse_args()


def contains_keywords(series: pd.Series, keywords: list[str]) -> pd.Series:
    text = series.astype(str).str.lower()
    mask = pd.Series(False, index=series.index)
    for keyword in keywords:
        mask |= text.str.contains(str(keyword).lower(), regex=False)
    return mask


def sample_obs(obs: pd.DataFrame, cfg: dict, rng: np.random.Generator, max_cells: int | None) -> np.ndarray:
    stage_col = cfg["input"]["stage_column"]
    cell_type_col = cfg["input"]["cell_type_column"]
    total = int(max_cells or cfg["sampling"]["max_cells_total"])
    per_group = int(cfg["sampling"]["max_cells_per_stage_celltype"])
    sampled: list[np.ndarray] = []
    for _, group in obs.groupby([stage_col, cell_type_col], observed=True, sort=False):
        take = min(per_group, len(group))
        sampled.append(rng.choice(group.index.to_numpy(), size=take, replace=False))
    selected = np.concatenate(sampled) if sampled else np.array([], dtype=int)
    if len(selected) > total:
        selected = rng.choice(selected, size=total, replace=False)
    return np.sort(selected.astype(int))


def gene_indices(adata: ad.AnnData, genes: list[str]) -> dict[str, int]:
    var = adata.var.copy()
    lookup: dict[str, int] = {}
    candidates = [adata.var_names.astype(str)]
    for column in ["gene_short_name", "feature_name", "gene_name"]:
        if column in var.columns:
            candidates.append(var[column].astype(str))
    for gene in genes:
        target = gene.lower()
        found = None
        for values in candidates:
            matches = np.where(values.str.lower().to_numpy() == target)[0]
            if len(matches):
                found = int(matches[0])
                break
        if found is not None:
            lookup[gene] = found
    return lookup


def matrix_to_dense(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix)


def safe_corr(x: pd.Series, y: pd.Series) -> float:
    values = pd.concat([pd.to_numeric(x, errors="coerce"), pd.to_numeric(y, errors="coerce")], axis=1).dropna()
    if len(values) < 3 or values.iloc[:, 0].nunique() < 2 or values.iloc[:, 1].nunique() < 2:
        return float("nan")
    return float(values.iloc[:, 0].corr(values.iloc[:, 1], method="spearman"))


def main() -> None:
    args = parse_args()
    cfg = load_json(args.config)
    rng = np.random.default_rng(int(cfg["seed"]))
    h5ad_path = resolve_project_path(cfg["input"]["embryo_h5ad"])
    if not h5ad_path.exists():
        raise FileNotFoundError(h5ad_path)

    adata = ad.read_h5ad(h5ad_path, backed="r")
    obs = adata.obs.reset_index(drop=True)
    selected = sample_obs(obs, cfg, rng, args.max_cells)
    sample = obs.iloc[selected].copy().reset_index(drop=True)
    sample["source_cell_index"] = selected
    sample["stage"] = sample[cfg["input"]["stage_column"]].map(canonical_stage)
    sample["stage_num"] = sample["stage"].map(stage_number)
    sample["cell_type"] = sample[cfg["input"]["cell_type_column"]].astype(str)
    sample["donor_id"] = sample[cfg["input"]["donor_column"]].astype(str)
    sample["major_cluster"] = sample[cfg["input"]["major_cluster_column"]].astype(str)

    marker_cfg = cfg["case_study"]
    genes = [marker_cfg["niche_marker_gene"], marker_cfg["target_gene"]]
    indices = gene_indices(adata, genes)
    if indices:
        X = matrix_to_dense(adata.X[selected, list(indices.values())])
        for pos, gene in enumerate(indices):
            sample[f"{gene}_expr"] = X[:, pos].astype(float)
    for gene in genes:
        if f"{gene}_expr" not in sample.columns:
            sample[f"{gene}_expr"] = np.nan

    niche_mask = contains_keywords(sample["cell_type"], marker_cfg["niche_cell_type_keywords"])
    target_mask = contains_keywords(sample["cell_type"], marker_cfg["target_cell_type_keywords"])
    sample["is_candidate_niche_cell"] = niche_mask
    sample["is_candidate_target_cell"] = target_mask

    q = float(marker_cfg["high_expression_quantile"])
    tbx4_col = f"{marker_cfg['niche_marker_gene']}_expr"
    glis3_col = f"{marker_cfg['target_gene']}_expr"
    tbx4_threshold = sample.loc[niche_mask, tbx4_col].quantile(q) if niche_mask.any() else np.nan
    glis3_threshold = sample.loc[target_mask, glis3_col].quantile(q) if target_mask.any() else np.nan
    sample["tbx4_high_niche"] = niche_mask & (sample[tbx4_col] >= tbx4_threshold)
    sample["glis3_high_target"] = target_mask & (sample[glis3_col] >= glis3_threshold)
    sample["tbx4_expr_in_niche"] = sample[tbx4_col].where(niche_mask)
    sample["glis3_expr_in_target"] = sample[glis3_col].where(target_mask)

    group_cols = ["stage", "stage_num", "donor_id"]
    total = sample.groupby(group_cols, as_index=False, observed=True).size().rename(columns={"size": "n_cells"})
    niche = (
        sample.groupby(group_cols, as_index=False, observed=True)
        .agg(
            niche_cells=("is_candidate_niche_cell", "sum"),
            target_cells=("is_candidate_target_cell", "sum"),
            tbx4_high_niche_cells=("tbx4_high_niche", "sum"),
            glis3_high_target_cells=("glis3_high_target", "sum"),
            mean_tbx4_niche=("tbx4_expr_in_niche", "mean"),
            mean_glis3_target=("glis3_expr_in_target", "mean"),
        )
        .merge(total, on=group_cols, how="left")
    )
    niche["niche_fraction"] = niche["niche_cells"] / niche["n_cells"].replace(0, np.nan)
    niche["target_fraction"] = niche["target_cells"] / niche["n_cells"].replace(0, np.nan)
    niche["tbx4_high_niche_fraction"] = niche["tbx4_high_niche_cells"] / niche["niche_cells"].replace(0, np.nan)
    niche["glis3_high_target_fraction"] = niche["glis3_high_target_cells"] / niche["target_cells"].replace(0, np.nan)
    niche["tbx4_glis3_niche_cooccurrence"] = (
        niche["tbx4_high_niche_fraction"].fillna(0.0) * niche["glis3_high_target_fraction"].fillna(0.0)
    )

    case = niche.sort_values("tbx4_glis3_niche_cooccurrence", ascending=False)
    autonomous_corr = safe_corr(sample.loc[target_mask, tbx4_col], sample.loc[target_mask, glis3_col])
    niche_corr = safe_corr(niche["tbx4_high_niche_fraction"], niche["glis3_high_target_fraction"])
    scores = pd.DataFrame(
        [
            {
                "case": "Tbx4-Glis3",
                "cell_autonomous_spearman": autonomous_corr,
                "niche_mediated_spearman": niche_corr,
                "interpretation": "niche_mediated_candidate" if pd.notna(niche_corr) and (pd.isna(autonomous_corr) or niche_corr > autonomous_corr) else "cell_autonomous_or_unresolved",
                "sampled_cells": int(len(sample)),
                "tbx4_gene_found": bool(marker_cfg["niche_marker_gene"] in indices),
                "glis3_gene_found": bool(marker_cfg["target_gene"] in indices),
            }
        ]
    )

    outputs = cfg["output"]
    niche_path = write_table(niche, outputs["niche_signature_by_stage_donor"])
    case_path = write_table(case, outputs["tbx4_glis3_niche_case"])
    score_path = write_table(scores, outputs["cell_autonomous_vs_niche"])
    print(f"Wrote niche signatures: {niche_path}")
    print(f"Wrote Tbx4-Glis3 case table: {case_path}")
    print(f"Wrote autonomous/niche scores: {score_path}")


if __name__ == "__main__":
    main()
