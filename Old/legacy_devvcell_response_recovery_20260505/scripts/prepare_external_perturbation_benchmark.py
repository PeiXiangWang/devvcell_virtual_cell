"""Download and inspect a small public scPerturb H5AD benchmark.

This script establishes a reproducible external perturbation-data entry point
for future DevVCell-vs-CellOT/scGen-style benchmarking. The default file is
Datlinger/Bock 2021 because it contains enough controls, guide pairs and TCR
stimulation contexts for a lightweight guide-transfer benchmark.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from devvcell.io import resolve_project_path, write_json  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record-id", default="10044268", help="Zenodo record ID for scPerturb H5AD files.")
    parser.add_argument(
        "--file-name",
        default="DatlingerBock2021.h5ad",
        help="H5AD file within the Zenodo record.",
    )
    parser.add_argument("--output-dir", default="data/external/scperturb")
    parser.add_argument("--force", action="store_true", help="Redownload even if the target file exists.")
    return parser.parse_args()


def fetch_record(record_id: str) -> dict:
    url = f"https://zenodo.org/api/records/{record_id}"
    with urllib.request.urlopen(url, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def find_file(record: dict, file_name: str) -> dict:
    for item in record.get("files", []):
        if item.get("key") == file_name:
            return item
    available = [item.get("key") for item in record.get("files", [])]
    raise FileNotFoundError(f"{file_name!r} not found in Zenodo record. Available examples: {available[:10]}")


def download_file(file_record: dict, target_path: Path, force: bool) -> None:
    if target_path.exists() and not force:
        return
    target_path.parent.mkdir(parents=True, exist_ok=True)
    url = file_record["links"]["self"]
    tmp_path = target_path.with_suffix(target_path.suffix + ".part")
    with urllib.request.urlopen(url, timeout=300) as response, tmp_path.open("wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
    tmp_path.replace(target_path)


def load_cached_source(manifest_path: Path) -> dict:
    if not manifest_path.exists():
        return {}
    with manifest_path.open("r", encoding="utf-8") as handle:
        cached = json.load(handle)
    source = cached.get("source", {})
    return source if isinstance(source, dict) else {}


def column_candidates(columns: list[str]) -> dict[str, list[str]]:
    lowered = {col: col.lower() for col in columns}
    return {
        "perturbation": [
            col
            for col, lower in lowered.items()
            if any(token in lower for token in ["perturb", "condition", "target", "guide", "grna", "sgrna", "gene"])
        ],
        "batch_or_sample": [
            col
            for col, lower in lowered.items()
            if any(token in lower for token in ["batch", "sample", "library", "gem", "well"])
        ],
        "cell_annotation": [
            col
            for col, lower in lowered.items()
            if any(token in lower for token in ["cell", "type", "cluster", "leiden", "louvain"])
        ],
    }


def inspect_h5ad(path: Path) -> dict:
    adata = ad.read_h5ad(path, backed="r")
    obs_columns = list(map(str, adata.obs.columns))
    var_columns = list(map(str, adata.var.columns))
    candidates = column_candidates(obs_columns)
    preview: dict[str, list[str]] = {}
    for col in candidates["perturbation"][:8]:
        series = adata.obs[col].astype(str)
        preview[col] = series.value_counts(dropna=False).head(10).index.astype(str).tolist()
    summary = {
        "path": str(path.relative_to(PROJECT_ROOT)).replace("\\", "/"),
        "size_bytes": int(path.stat().st_size),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "obs_columns": obs_columns,
        "var_columns": var_columns,
        "candidate_columns": candidates,
        "perturbation_value_preview": preview,
    }
    adata.file.close()
    return summary


def main() -> None:
    args = parse_args()
    out_dir = resolve_project_path(args.output_dir)
    target_path = out_dir / args.file_name
    manifest_path = out_dir / "scperturb_benchmark_manifest.json"
    fetch_warning = None
    file_record: dict[str, object] = {}
    try:
        record = fetch_record(args.record_id)
        file_record = find_file(record, args.file_name)
        download_file(file_record, target_path, args.force)
    except (urllib.error.URLError, TimeoutError, FileNotFoundError) as exc:
        if not target_path.exists():
            raise
        fetch_warning = f"{type(exc).__name__}: {exc}"
        cached_source = load_cached_source(manifest_path)
        file_record = {
            "checksum": cached_source.get("selected_file_md5"),
            "size": cached_source.get("selected_file_size_bytes", target_path.stat().st_size),
        }
    h5ad_summary = inspect_h5ad(target_path)

    source = {
        "name": "scPerturb Single-Cell Perturbation Data",
        "zenodo_record_id": args.record_id,
        "zenodo_url": f"https://zenodo.org/records/{args.record_id}",
        "selected_file": args.file_name,
        "selected_file_md5": file_record.get("checksum"),
        "selected_file_size_bytes": int(file_record.get("size") or target_path.stat().st_size),
        "associated_publication": "Peidli et al., Nature Methods 2024, doi:10.1038/s41592-023-02144-y",
    }
    if fetch_warning is not None:
        source["zenodo_fetch_warning"] = fetch_warning

    manifest = {
        "analysis": "external_scperturb_benchmark_ingest",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source": source,
        "h5ad": h5ad_summary,
        "next_steps": [
            "Choose perturbation/control columns after manual schema inspection.",
            "Map perturbation labels to gene symbols and shared genes with DevVCell.",
            "Implement control-to-perturbed response benchmark against scGen/CellOT-style baselines.",
            "Report perturbation-specific effects rather than global batch or composition differences.",
        ],
    }
    write_json(manifest_path, manifest)
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
