"""Validate DevVCell external perturbation and response-recovery outputs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from devvcell.io import load_json, resolve_project_path  # noqa: E402
from devvcell.tables import read_table, write_table  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--response-config", default="config/response_recovery.json")
    parser.add_argument("--transfer-config", default="config/perturbation_transfer.json")
    return parser.parse_args()


def existing(path_like: str) -> Path | None:
    path = resolve_project_path(path_like)
    if path.exists():
        return path
    if path.suffix == ".parquet" and path.with_suffix(".csv").exists():
        return path.with_suffix(".csv")
    return None


def main() -> None:
    args = parse_args()
    response_cfg = load_json(args.response_config)
    transfer_cfg = load_json(args.transfer_config)
    rows: list[dict[str, object]] = []

    checks = {
        "external_response_dictionary": transfer_cfg["output"]["external_response_dictionary"],
        "transferred_response_by_stage_celltype": transfer_cfg["output"]["transferred_response_by_stage_celltype"],
        "response_recovery_classes": response_cfg["output"]["response_recovery_classes"],
        "stage_vulnerability_response_recovery": response_cfg["output"]["stage_vulnerability_response_recovery"],
        "window_enrichment_statistics": response_cfg["output"].get("window_enrichment_statistics", ""),
        "minimal_rescue_control_matrix": response_cfg["output"].get("minimal_rescue_control_matrix", ""),
    }
    for name, path_like in checks.items():
        path = existing(path_like)
        if path is None:
            rows.append({"check": name, "status": "missing", "rows": 0, "details": str(resolve_project_path(path_like))})
            continue
        frame = read_table(path)
        numeric = frame.select_dtypes(include=[np.number])
        finite_fraction = float(np.isfinite(numeric.to_numpy()).mean()) if numeric.size else 1.0
        pass_threshold = 0.90 if name == "window_enrichment_statistics" else 0.99
        rows.append(
            {
                "check": name,
                "status": "pass" if len(frame) > 0 and finite_fraction >= pass_threshold else "warn",
                "rows": int(len(frame)),
                "numeric_finite_fraction": finite_fraction,
                "details": str(path),
            }
        )

    classes_path = existing(response_cfg["output"]["response_recovery_classes"])
    if classes_path is not None:
        classes = read_table(classes_path)
        if "response_recovery_class" in classes.columns:
            for klass, count in classes["response_recovery_class"].value_counts().items():
                rows.append({"check": f"class_count:{klass}", "status": "info", "rows": int(count), "details": "response_recovery_classes"})

    output = resolve_project_path(response_cfg["output"]["tables"]) / "external_validation_report.csv"
    path = write_table(pd.DataFrame(rows), output)
    print(f"Wrote validation report: {path}")


if __name__ == "__main__":
    main()
