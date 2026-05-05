"""Organoid or donor heterogeneity control summaries."""

from __future__ import annotations

import pandas as pd


def summarize_false_positive_by_sample(
    frame: pd.DataFrame,
    *,
    sample_col: str = "sample_id",
    class_col: str = "normality_class",
) -> pd.DataFrame:
    """Estimate control false-positive rates by heldout sample labels."""

    if sample_col not in frame.columns:
        raise ValueError(f"Missing sample column: {sample_col}")
    grouped = frame.groupby(sample_col, dropna=False)
    rows = []
    for sample_id, sample in grouped:
        n_cells = sample.shape[0]
        n_abnormal = int((sample[class_col] == "abnormal_off_normal").sum()) if class_col in sample.columns else 0
        rows.append(
            {
                sample_col: sample_id,
                "n_cells": n_cells,
                "n_abnormal_off_normal": n_abnormal,
                "false_positive_rate": n_abnormal / n_cells if n_cells else 0.0,
            }
        )
    return pd.DataFrame(rows)
