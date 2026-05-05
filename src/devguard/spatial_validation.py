"""Lightweight spatial/tomo-seq validation summaries."""

from __future__ import annotations

import pandas as pd


def axis_enrichment(
    frame: pd.DataFrame,
    *,
    axis_col: str,
    class_col: str = "normality_class",
) -> pd.DataFrame:
    if axis_col not in frame.columns:
        raise ValueError(f"Missing axis column: {axis_col}")
    counts = frame.groupby([axis_col, class_col], dropna=False).size().reset_index(name="n_cells")
    totals = counts.groupby(axis_col, dropna=False)["n_cells"].transform("sum")
    counts["fraction_within_axis_bin"] = counts["n_cells"] / totals
    return counts
