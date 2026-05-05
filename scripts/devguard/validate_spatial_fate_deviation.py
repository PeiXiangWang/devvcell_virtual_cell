from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import pandas as pd

from devguard.io import ensure_dir, write_dataframe, write_manifest
from devguard.spatial_validation import axis_enrichment


def validate_spatial_fate_deviation(
    input_csv: str | Path,
    output_dir: str | Path = "results/devguard/spatial_validation",
    axis_col: str = "tomo_axis_bin",
) -> Path:
    output = ensure_dir(output_dir)
    frame = pd.read_csv(input_csv)
    output_path = output / "fate_deviation_spatial_axis.csv"
    if axis_col in frame.columns:
        result = axis_enrichment(frame, axis_col=axis_col)
        status = "completed"
    else:
        result = pd.DataFrame(
            [
                {
                    "status": "skipped",
                    "reason": f"Input table does not contain spatial/tomo axis column {axis_col}.",
                }
            ]
        )
        status = "skipped_missing_axis"
    write_dataframe(result, output_path)
    write_manifest(
        output / "spatial_validation_manifest.json",
        name="validate_spatial_fate_deviation",
        inputs=[str(input_csv)],
        outputs=[str(output_path)],
        parameters={"axis_col": axis_col},
        metrics={"status": status},
    )
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate fate deviation against spatial/tomo axis labels.")
    parser.add_argument("--input", default="results/devguard/perturbation_classification/cell_normality_classes.csv")
    parser.add_argument("--output-dir", default="results/devguard/spatial_validation")
    parser.add_argument("--axis-col", default="tomo_axis_bin")
    args = parser.parse_args()
    validate_spatial_fate_deviation(args.input, args.output_dir, args.axis_col)


if __name__ == "__main__":
    main()
