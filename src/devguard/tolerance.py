"""Developmental Tolerance Index summaries."""

from __future__ import annotations

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
