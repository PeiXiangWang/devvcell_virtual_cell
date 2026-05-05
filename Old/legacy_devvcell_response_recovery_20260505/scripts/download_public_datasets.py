"""Download registered public DevVCell datasets with size guards and checksums."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
import urllib.request
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from devvcell.io import load_json, resolve_project_path  # noqa: E402
from devvcell.tables import write_table  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", default="config/external_datasets.json")
    parser.add_argument("--dataset", action="append", help="Dataset id to download. Can be passed multiple times.")
    parser.add_argument("--metadata-only", action="store_true", help="Download only files <= max_mb.")
    parser.add_argument("--max-mb", type=float, default=50.0, help="Default safety limit for metadata downloads.")
    parser.add_argument("--destination", default=None)
    return parser.parse_args()


def parse_size_mb(size: str | None) -> float | None:
    if not size:
        return None
    match = re.match(r"\s*([0-9.]+)\s*([KMG]b)\s*$", str(size), flags=re.IGNORECASE)
    if not match:
        return None
    value = float(match.group(1))
    unit = match.group(2).lower()
    if unit == "kb":
        return value / 1024.0
    if unit == "gb":
        return value * 1024.0
    return value


def iter_files(registry: dict, selected: set[str] | None) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for section in ["primary_candidates", "normal_manifold_validation", "gonad_niche_validation"]:
        for item in registry.get(section, []):
            dataset_id = str(item.get("dataset_id", ""))
            if selected and dataset_id not in selected:
                continue
            for file_item in item.get("processed_files", []):
                if file_item.get("url"):
                    rows.append(
                        {
                            "section": section,
                            "dataset_id": dataset_id,
                            "title": item.get("title"),
                            "file_name": file_item.get("name"),
                            "file_size": file_item.get("size"),
                            "file_size_mb": parse_size_mb(file_item.get("size")),
                            "url": file_item.get("url"),
                        }
                    )
    return rows


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def download(url: str, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as response, output.open("wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)


def main() -> None:
    args = parse_args()
    registry = load_json(args.registry)
    selected = set(args.dataset) if args.dataset else None
    destination = resolve_project_path(args.destination or registry["download_policy"]["default_destination"])
    rows = iter_files(registry, selected)
    manifest_rows: list[dict[str, object]] = []

    for row in rows:
        size_mb = row.get("file_size_mb")
        if args.metadata_only and size_mb is not None and float(size_mb) > args.max_mb:
            manifest_rows.append({**row, "status": "skipped_size_guard", "local_path": "", "sha256": ""})
            continue
        output = destination / str(row["dataset_id"]) / str(row["file_name"])
        if not output.exists():
            download(str(row["url"]), output)
        manifest_rows.append(
            {
                **row,
                "status": "downloaded" if output.exists() else "missing",
                "local_path": str(output),
                "bytes": int(output.stat().st_size) if output.exists() else 0,
                "sha256": sha256(output) if output.exists() else "",
                "gzip_valid": gzip_valid(output) if output.suffix == ".gz" and output.exists() else "",
            }
        )

    manifest = pd.DataFrame(manifest_rows)
    path = write_table(manifest, destination / "download_manifest.csv")
    print(json.dumps({"manifest": str(path), "n_files": len(manifest)}, indent=2))


def gzip_valid(path: Path) -> bool:
    try:
        with gzip.open(path, "rb") as handle:
            handle.read(1)
        return True
    except OSError:
        return False


if __name__ == "__main__":
    main()
