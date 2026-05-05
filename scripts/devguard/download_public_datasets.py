from __future__ import annotations

import argparse
import hashlib
import time
import urllib.request
from pathlib import Path

import _bootstrap  # noqa: F401
import pandas as pd

from devguard.io import dataset_registry_frame, ensure_dir, load_json, write_dataframe, write_manifest


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(url: str, output: Path, timeout: int = 60) -> str:
    ensure_dir(output.parent)
    tmp = output.with_suffix(output.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "DevGuard/0.1"})
    with urllib.request.urlopen(request, timeout=timeout) as response, tmp.open("wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
    tmp.replace(output)
    return sha256_file(output)


def download_registered_files(
    config_path: str | Path,
    *,
    datasets: set[str] | None = None,
    allow_large: bool = False,
    max_mb: float = 100.0,
) -> pd.DataFrame:
    config = load_json(config_path)
    registry = dataset_registry_frame(config)
    rows = []
    for _, dataset in registry.iterrows():
        dataset_id = dataset["dataset_id"]
        if datasets and dataset_id not in datasets:
            continue
        for public_file in dataset.get("public_files", []) or []:
            local_path = Path(public_file["local_path"])
            size_mb = float(public_file.get("size_mb", 0))
            status = "pending"
            sha256 = ""
            error = ""
            if local_path.exists():
                status = "present"
                sha256 = sha256_file(local_path)
            elif size_mb > max_mb and not allow_large:
                status = "skipped_size_guard"
            else:
                try:
                    started = time.time()
                    sha256 = download_file(public_file["url"], local_path)
                    status = "downloaded"
                    error = f"elapsed_seconds={time.time() - started:.1f}"
                except Exception as exc:  # pragma: no cover - network dependent
                    status = "failed"
                    error = str(exc)
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "file_name": public_file["name"],
                    "url": public_file["url"],
                    "local_path": str(local_path),
                    "size_mb": size_mb,
                    "status": status,
                    "sha256": sha256,
                    "note": error,
                }
            )
    result = pd.DataFrame(rows)
    output = Path("results/devguard/dataset_metadata/public_download_manifest.csv")
    write_dataframe(result, output)
    write_manifest(
        "results/devguard/dataset_metadata/public_download_manifest.json",
        name="download_public_datasets",
        inputs=[str(config_path)],
        outputs=[str(output)],
        parameters={"allow_large": allow_large, "max_mb": max_mb, "datasets": sorted(datasets or [])},
        metrics={"n_files": int(result.shape[0]), "n_downloaded": int((result["status"] == "downloaded").sum()) if not result.empty else 0},
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Download public DevGuard dataset files registered in datasets_mouse.json.")
    parser.add_argument("--config", default="config/devguard/datasets_mouse.json")
    parser.add_argument("--datasets", nargs="*", help="Optional dataset IDs to download.")
    parser.add_argument("--allow-large", action="store_true", help="Allow files larger than --max-mb.")
    parser.add_argument("--max-mb", type=float, default=100.0)
    args = parser.parse_args()
    download_registered_files(args.config, datasets=set(args.datasets or []), allow_large=args.allow_large, max_mb=args.max_mb)


if __name__ == "__main__":
    main()
