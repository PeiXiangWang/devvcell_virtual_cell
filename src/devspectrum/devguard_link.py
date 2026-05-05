"""Link DevSpectrum residuals to DevGuard failure modes."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def class_summary_from_cells(frame: pd.DataFrame) -> pd.DataFrame:
    counts = frame.groupby(["perturbation_name", "lineage", "normality_class"], dropna=False, observed=True).size().reset_index(name="n_cells")
    totals = counts.groupby(["perturbation_name", "lineage"], dropna=False, observed=True)["n_cells"].transform("sum")
    counts["fraction"] = counts["n_cells"] / totals.where(totals > 0, np.nan)
    return counts


def failure_mode_spectral_signature(residuals: pd.DataFrame) -> pd.DataFrame:
    return (
        residuals.groupby(["cohort", "normality_class", "module_name"], dropna=False, observed=True)
        .agg(
            mean_residual=("spectral_residual", "mean"),
            mean_abs_residual=("absolute_spectral_residual", "mean"),
            median_abs_residual=("absolute_spectral_residual", "median"),
            max_abs_residual=("absolute_spectral_residual", "max"),
            n_records=("module_name", "size"),
            n_cells=("n_cells", "sum"),
        )
        .reset_index()
        .sort_values(["cohort", "normality_class", "mean_abs_residual"], ascending=[True, True, False])
    )


def correlate_dti_with_spectral_residuals(residual_summary: pd.DataFrame, dti_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for cohort, dti in dti_tables.items():
        if dti.empty:
            continue
        left = residual_summary[residual_summary["cohort"].eq(cohort)].copy()
        merged = left.merge(dti, on="lineage", how="inner")
        for metric in ["global_spectral_distance", "mean_absolute_residual", "max_absolute_residual"]:
            if metric not in merged.columns or "DTI" not in merged.columns or merged.shape[0] < 3:
                rho = np.nan
                p_value = np.nan
            else:
                result = spearmanr(merged["DTI"], merged[metric], nan_policy="omit")
                rho = float(result.statistic)
                p_value = float(result.pvalue)
            rows.append(
                {
                    "cohort": cohort,
                    "devguard_metric": "DTI",
                    "spectral_metric": metric,
                    "spearman_rho": rho,
                    "p_value": p_value,
                    "n_lineages": int(merged.shape[0]),
                }
            )
    return pd.DataFrame(rows)


def write_link_summary(correlation: pd.DataFrame, signature: pd.DataFrame) -> str:
    top = signature.sort_values("mean_abs_residual", ascending=False).head(10)
    lines = [
        "# DevGuard-DevSpectrum Link Summary",
        "",
        "## DTI correlations",
        "",
        correlation.to_markdown(index=False) if not correlation.empty else "No correlations were computed.",
        "",
        "## Top failure-mode spectral signatures",
        "",
        top.to_markdown(index=False) if not top.empty else "No failure-mode signatures were computed.",
        "",
    ]
    return "\n".join(lines)
