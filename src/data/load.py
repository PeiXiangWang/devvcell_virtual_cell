from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd

from .schema import inspect_h5ad, inspect_mtx_bundle, markdown_audit, summarize_records


DATA_EXTENSIONS = {
    ".h5ad",
    ".loom",
    ".zarr",
    ".mtx",
    ".csv",
    ".tsv",
    ".parquet",
    ".rds",
    ".rda",
}


def discover_data_files(roots: list[str] | None = None) -> list[Path]:
    roots = roots or ["data", "datasets", "input", "raw", "processed", "notebooks"]
    files: list[Path] = []
    for root in roots:
        base = Path(root)
        if not base.exists():
            continue
        for path in base.rglob("*"):
            if path.is_file():
                suffixes = "".join(path.suffixes[-2:]).lower()
                if path.suffix.lower() in DATA_EXTENSIONS or suffixes in {".rds.gz", ".rda.gz", ".txt.gz"}:
                    files.append(path)
    return sorted(set(files))


def audit_data_files(files: list[Path], deep_h5ad_size_limit: int = 1_100_000_000) -> tuple[list[dict[str, Any]], dict[str, Any], str]:
    records: list[dict[str, Any]] = []
    for path in files:
        low = path.name.lower()
        if low.endswith(".h5ad"):
            if path.stat().st_size > deep_h5ad_size_limit and path.name != "scLine_pro.h5ad":
                records.append(
                    {
                        "path": str(path),
                        "format": "h5ad",
                        "readable": False,
                        "size_bytes": int(path.stat().st_size),
                        "note": "deep inspection skipped for >1.1GB file during quick audit; file remains catalogued for full audit",
                    }
                )
            elif path.name == "scLine_pro.h5ad":
                records.append(
                    {
                        "path": str(path),
                        "format": "h5ad",
                        "readable": False,
                        "size_bytes": int(path.stat().st_size),
                        "note": "large 12GB atlas catalogued; use data/processed/cell_level_subset_v1.h5ad for default deep quick audit",
                        "time_fields": ["author_day", "development_stage"],
                        "cell_type_fields": ["author_major_cell_cluster", "author_cell_type", "cell_type"],
                        "condition_fields": ["donor_id", "author_experimental_id", "sex"],
                    }
                )
            else:
                records.append(inspect_h5ad(path))
        elif low.endswith(".mtx"):
            records.append(inspect_mtx_bundle(path))
        else:
            records.append({"path": str(path), "format": path.suffix.lower().lstrip("."), "readable": False, "note": "catalogued but not deeply inspected"})
    summary = summarize_records(records)
    return records, summary, markdown_audit(records, summary)


def load_h5ad(path: str | Path) -> ad.AnnData:
    return ad.read_h5ad(path)


def first_existing_column(frame: pd.DataFrame, candidates: list[str]) -> str | None:
    lower = {str(c).lower(): str(c) for c in frame.columns}
    for candidate in candidates:
        if candidate in frame.columns:
            return candidate
        hit = lower.get(str(candidate).lower())
        if hit is not None:
            return hit
    for candidate in candidates:
        token = str(candidate).lower()
        for col in frame.columns:
            if token in str(col).lower():
                return str(col)
    return None
