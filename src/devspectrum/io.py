"""I/O helpers for DevSpectrum."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from devguard.io import ensure_dir, load_json, read_h5ad, write_dataframe, write_json, write_manifest


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, low_memory=False)


def resolve_output(root: str | Path, *parts: str) -> Path:
    return ensure_dir(Path(root).joinpath(*parts))


def load_config(path: str | Path) -> dict[str, Any]:
    return load_json(path)


__all__ = [
    "ensure_dir",
    "load_config",
    "read_h5ad",
    "read_table",
    "resolve_output",
    "write_dataframe",
    "write_json",
    "write_manifest",
]
