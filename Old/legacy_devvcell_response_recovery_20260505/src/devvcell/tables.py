"""Table helpers for DevVCell command-line analyses."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from devvcell.io import resolve_project_path


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def read_table(path_like: str | Path) -> pd.DataFrame:
    """Read CSV/TSV/Parquet by suffix."""

    path = resolve_project_path(path_like)
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def write_table(frame: pd.DataFrame, path_like: str | Path) -> Path:
    """Write a table, falling back from parquet to CSV if no parquet engine exists."""

    path = resolve_project_path(path_like)
    ensure_parent(path)
    if path.suffix.lower() == ".parquet":
        try:
            frame.to_parquet(path, index=False)
            return path
        except ImportError:
            fallback = path.with_suffix(".csv")
            frame.to_csv(fallback, index=False)
            return fallback
    if path.suffix.lower() in {".tsv", ".txt"}:
        frame.to_csv(path, index=False, sep="\t")
    else:
        frame.to_csv(path, index=False)
    return path


def minmax(values: pd.Series) -> pd.Series:
    series = pd.to_numeric(values, errors="coerce").astype(float)
    lo = series.min(skipna=True)
    hi = series.max(skipna=True)
    if not np.isfinite(lo) or not np.isfinite(hi) or np.isclose(lo, hi):
        return pd.Series(np.zeros(len(series)), index=series.index)
    return (series - lo) / (hi - lo)


def latent_columns(frame: pd.DataFrame, prefixes: Iterable[str] = ("latent_",)) -> list[str]:
    columns: list[str] = []
    for column in frame.columns:
        for prefix in prefixes:
            if column.startswith(prefix):
                columns.append(column)
                break
    return sorted(columns, key=lambda item: (len(item), item))


def numeric_matrix(frame: pd.DataFrame, columns: list[str]) -> np.ndarray:
    if not columns:
        raise ValueError("No numeric columns were supplied.")
    return frame[columns].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float)
