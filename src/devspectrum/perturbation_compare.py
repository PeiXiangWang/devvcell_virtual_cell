"""Control-vs-perturbation spectral comparison and endpoint projection."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from devspectrum.modules import load_module_registry, score_modules


def _subset_by_classified_cells(adata, classes: pd.DataFrame):
    ids = classes["cell_id"].astype(str).tolist()
    if "cell_id" not in adata.obs.columns:
        raise ValueError("AnnData obs must contain cell_id for endpoint projection.")
    local = adata.copy()
    local.obs_names = local.obs["cell_id"].astype(str)
    missing = sorted(set(ids) - set(local.obs_names.astype(str)))
    if missing:
        raise ValueError(f"Missing {len(missing)} classified cells in {adata}: {missing[:5]}")
    return local[ids].copy()


def score_classified_endpoint(
    h5ad_path: str | Path,
    classes_csv: str | Path,
    *,
    cohort: str,
    module_registry: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    classes = pd.read_csv(classes_csv, low_memory=False)
    adata = ad.read_h5ad(h5ad_path)
    adata = _subset_by_classified_cells(adata, classes)
    modules = load_module_registry(module_registry)
    scores, coverage = score_modules(adata, modules)
    metadata_cols = ["cell_id", "sample_id", "time_point", "time_numeric", "lineage", "normality_class", "perturbation_name"]
    metadata = classes[[column for column in metadata_cols if column in classes.columns]].copy()
    metadata.index = adata.obs_names.astype(str)
    long = (
        scores.reset_index(names="cell_id")
        .melt(id_vars="cell_id", var_name="module_name", value_name="feature_value")
        .merge(metadata.reset_index(drop=True), on="cell_id", how="left")
    )
    long["cohort"] = cohort
    return long, coverage


def endpoint_residuals(
    control_scores: pd.DataFrame,
    perturbation_scores: pd.DataFrame,
    *,
    min_cells_per_group: int = 10,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    group_cols = ["cohort", "lineage", "module_name", "normality_class"]
    pert = (
        perturbation_scores.groupby(group_cols, dropna=False, observed=True)
        .agg(
            observed_value=("feature_value", "mean"),
            observed_std=("feature_value", "std"),
            n_cells=("cell_id", "size"),
            n_samples=("sample_id", "nunique"),
        )
        .reset_index()
    )
    pert = pert[pert["n_cells"] >= int(min_cells_per_group)].copy()
    control = (
        control_scores.groupby(["lineage", "module_name"], dropna=False, observed=True)
        .agg(expected_control_value=("feature_value", "mean"), control_std=("feature_value", "std"), control_cells=("cell_id", "size"))
        .reset_index()
    )
    residuals = pert.merge(control, on=["lineage", "module_name"], how="left")
    global_control = (
        control_scores.groupby(["module_name"], dropna=False, observed=True)
        .agg(global_expected_control_value=("feature_value", "mean"), global_control_std=("feature_value", "std"))
        .reset_index()
    )
    residuals = residuals.merge(global_control, on="module_name", how="left")
    residuals["expected_control_value"] = residuals["expected_control_value"].fillna(residuals["global_expected_control_value"])
    residuals["control_std"] = residuals["control_std"].fillna(residuals["global_control_std"]).replace(0, np.nan)
    residuals["spectral_residual"] = residuals["observed_value"] - residuals["expected_control_value"]
    residuals["standardized_residual"] = residuals["spectral_residual"] / residuals["control_std"].fillna(1.0)
    residuals["absolute_spectral_residual"] = residuals["spectral_residual"].abs()
    summary = (
        residuals.groupby(["cohort", "lineage"], dropna=False, observed=True)
        .agg(
            global_spectral_distance=("spectral_residual", lambda values: float(np.linalg.norm(values))),
            mean_absolute_residual=("absolute_spectral_residual", "mean"),
            max_absolute_residual=("absolute_spectral_residual", "max"),
            n_modules=("module_name", "nunique"),
            n_cells=("n_cells", "sum"),
        )
        .reset_index()
    )
    return residuals.sort_values(["cohort", "absolute_spectral_residual"], ascending=[True, False]), summary


def perturbation_fingerprint(residuals: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cohort, group in residuals.groupby("cohort", dropna=False, observed=True):
        module_abs = group.groupby("module_name", dropna=False, observed=True)["absolute_spectral_residual"].mean()
        weighted_mean = float((group["absolute_spectral_residual"] * group["n_cells"]).sum() / max(group["n_cells"].sum(), 1))
        class_records = group[["lineage", "normality_class", "n_cells"]].drop_duplicates()
        total_cells = float(class_records["n_cells"].sum())
        failure_cells = float(
            class_records.loc[
                class_records["normality_class"].astype(str).isin(["fate_deviation", "abnormal_off_normal"]),
                "n_cells",
            ].sum()
        )
        failure_fraction = failure_cells / max(total_cells, 1.0)
        rows.append(
            {
                "perturbation_name": cohort,
                "global_spectral_distance": float(np.linalg.norm(group["spectral_residual"].fillna(0).to_numpy(dtype=float))),
                "mean_absolute_residual": float(group["absolute_spectral_residual"].mean()),
                "cell_weighted_mean_absolute_residual": weighted_mean,
                "low_frequency_suppression": float((-group["spectral_residual"]).clip(lower=0).mean()),
                "high_frequency_burst_score": float(group["absolute_spectral_residual"].quantile(0.9)),
                "spectral_entropy_increase": float(module_abs.std() / max(module_abs.mean(), 1e-12)),
                "devguard_failure_fraction": failure_fraction,
                "devguard_linked_spectral_failure_burden": weighted_mean * failure_fraction,
                "top_residual_module": str(module_abs.sort_values(ascending=False).index[0]) if not module_abs.empty else "",
                "n_lineage_module_classes": int(group.shape[0]),
            }
        )
    return pd.DataFrame(rows).sort_values("global_spectral_distance", ascending=False)
