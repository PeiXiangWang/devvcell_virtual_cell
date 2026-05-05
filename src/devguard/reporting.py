"""Reporting helpers for pipeline outputs."""

from __future__ import annotations

import pandas as pd


def sample_count_table(obs: pd.DataFrame) -> pd.DataFrame:
    columns = ["dataset_id", "condition", "time_point", "lineage", "sample_id"]
    available = [column for column in columns if column in obs.columns]
    return obs.groupby(available, dropna=False).size().reset_index(name="n_cells")


def assert_no_rdeg_inputs(paths: list[str]) -> None:
    forbidden = ["scLine_pro.h5ad", "rdeg", "response_recovery", "devvcell_lite"]
    lowered = " ".join(str(path).lower() for path in paths)
    hits = [item for item in forbidden if item.lower() in lowered]
    if hits:
        raise ValueError(f"DevGuard mainline inputs must not include legacy/RDEG assets: {hits}")
