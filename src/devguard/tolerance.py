"""Developmental Tolerance Index summaries."""

from __future__ import annotations

import numpy as np
import pandas as pd

CLASS_COLUMNS = [
    "within_stage_normal",
    "developmental_delay",
    "developmental_acceleration",
    "fate_deviation",
    "abnormal_off_normal",
]


def compute_class_rates(frame: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    counts = frame.groupby(group_cols + ["normality_class"], dropna=False).size().reset_index(name="n_cells")
    pivot = counts.pivot_table(
        index=group_cols,
        columns="normality_class",
        values="n_cells",
        fill_value=0,
        aggfunc="sum",
    ).reset_index()
    for column in CLASS_COLUMNS:
        if column not in pivot.columns:
            pivot[column] = 0
    pivot["total_cells"] = pivot[CLASS_COLUMNS].sum(axis=1)
    for column in CLASS_COLUMNS:
        pivot[f"R_{column}"] = pivot[column] / pivot["total_cells"].where(pivot["total_cells"] > 0, 1)
    return pivot


def compute_developmental_tolerance_index(frame: pd.DataFrame, group_cols: list[str] | None = None) -> pd.DataFrame:
    if group_cols is None:
        group_cols = ["time_point", "lineage", "perturbation_name"]
    rates = compute_class_rates(frame, group_cols)
    rates["DTI"] = (
        rates["R_within_stage_normal"]
        - rates["R_developmental_delay"]
        - rates["R_developmental_acceleration"]
        - rates["R_fate_deviation"]
        - rates["R_abnormal_off_normal"]
    )
    return rates.sort_values("DTI", ascending=True).reset_index(drop=True)


def vulnerable_windows(dti: pd.DataFrame, n: int = 20) -> pd.DataFrame:
    return dti.sort_values("DTI", ascending=True).head(n).reset_index(drop=True)


def _single_group_dti(sample: pd.DataFrame) -> float:
    total = max(int(sample.shape[0]), 1)
    counts = sample["normality_class"].value_counts()
    return float(
        (
            counts.get("within_stage_normal", 0)
            - counts.get("developmental_delay", 0)
            - counts.get("developmental_acceleration", 0)
            - counts.get("fate_deviation", 0)
            - counts.get("abnormal_off_normal", 0)
        )
        / total
    )


def bootstrap_dti_ci(
    frame: pd.DataFrame,
    group_cols: list[str] | None = None,
    *,
    sample_col: str = "sample_id",
    n_bootstrap: int = 200,
    ci: float = 0.95,
    seed: int = 42,
) -> pd.DataFrame:
    """Estimate DTI confidence intervals using sample-level bootstrap when possible."""

    if group_cols is None:
        group_cols = ["time_point", "lineage", "perturbation_name"]
    rng = np.random.default_rng(seed)
    rows = []
    lower_q = (1.0 - ci) / 2.0
    upper_q = 1.0 - lower_q
    for keys, group in frame.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        if sample_col in group.columns and group[sample_col].nunique(dropna=False) > 1:
            sample_values = group[sample_col].astype("string").fillna("NA").astype(str)
            units = sample_values.unique()
            bootstrap_values = []
            for _ in range(n_bootstrap):
                sampled_units = rng.choice(units, size=len(units), replace=True)
                pieces = [group[sample_values.eq(unit)] for unit in sampled_units]
                bootstrap_values.append(_single_group_dti(pd.concat(pieces, axis=0)))
            bootstrap_level = "sample"
            n_bootstrap_units = int(len(units))
        else:
            bootstrap_values = []
            positions = np.arange(group.shape[0])
            for _ in range(n_bootstrap):
                sampled_positions = rng.choice(positions, size=positions.size, replace=True)
                bootstrap_values.append(_single_group_dti(group.iloc[sampled_positions]))
            bootstrap_level = "cell"
            n_bootstrap_units = int(group.shape[0])
        values = np.asarray(bootstrap_values, dtype=float)
        row = {column: value for column, value in zip(group_cols, keys)}
        row.update(
            {
                "DTI_ci_lower": float(np.quantile(values, lower_q)),
                "DTI_ci_upper": float(np.quantile(values, upper_q)),
                "DTI_bootstrap_mean": float(np.mean(values)),
                "n_bootstrap": int(n_bootstrap),
                "bootstrap_level": bootstrap_level,
                "n_bootstrap_units": n_bootstrap_units,
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)
