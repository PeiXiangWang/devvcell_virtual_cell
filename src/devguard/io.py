"""Input/output helpers for DevGuard."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


def ensure_dir(path: str | Path) -> Path:
    """Create a directory and return it as a Path."""

    resolved = Path(path)
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def load_json(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(payload: dict[str, Any], path: str | Path) -> Path:
    output = Path(path)
    ensure_dir(output.parent)
    with output.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return output


def read_h5ad(path: str | Path):
    import anndata as ad

    return ad.read_h5ad(path)


def write_h5ad(adata, path: str | Path) -> Path:
    output = Path(path)
    ensure_dir(output.parent)
    adata.write_h5ad(output)
    return output


def write_dataframe(frame: pd.DataFrame, path: str | Path, index: bool = False) -> Path:
    output = Path(path)
    ensure_dir(output.parent)
    suffix = output.suffix.lower()
    if suffix == ".parquet":
        frame.to_parquet(output, index=index)
    else:
        frame.to_csv(output, index=index)
    return output


def dataset_registry_frame(config: dict[str, Any]) -> pd.DataFrame:
    return pd.DataFrame(config.get("datasets", []))


def write_manifest(
    path: str | Path,
    *,
    name: str,
    inputs: list[str] | None = None,
    outputs: list[str] | None = None,
    parameters: dict[str, Any] | None = None,
    metrics: dict[str, Any] | None = None,
) -> Path:
    payload = {
        "name": name,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": inputs or [],
        "outputs": outputs or [],
        "parameters": parameters or {},
        "metrics": metrics or {},
    }
    return write_json(payload, path)
