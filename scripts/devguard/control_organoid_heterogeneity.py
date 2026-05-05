from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import pandas as pd

from devguard.heterogeneity import summarize_false_positive_by_sample
from devguard.io import ensure_dir, write_dataframe, write_manifest


def control_organoid_heterogeneity(
    input_csv: str | Path = "results/devguard/perturbation_classification/cell_normality_classes.csv",
    output_dir: str | Path = "results/devguard/heterogeneity_control",
) -> Path:
    output = ensure_dir(output_dir)
    frame = pd.read_csv(input_csv)
    summary = summarize_false_positive_by_sample(frame)
    summary_path = output / "organoid_false_positive_rate.csv"
    write_dataframe(summary, summary_path)
    write_manifest(
        output / "organoid_heterogeneity_manifest.json",
        name="control_organoid_heterogeneity",
        inputs=[str(input_csv)],
        outputs=[str(summary_path)],
        parameters={"note": "MVP summary. Full leave-one-organoid-out requires GSE212050 barcode metadata."},
        metrics={"n_samples": int(summary.shape[0])},
    )
    return summary_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize DevGuard organoid heterogeneity controls.")
    parser.add_argument("--input", default="results/devguard/perturbation_classification/cell_normality_classes.csv")
    parser.add_argument("--output-dir", default="results/devguard/heterogeneity_control")
    args = parser.parse_args()
    control_organoid_heterogeneity(args.input, args.output_dir)


if __name__ == "__main__":
    main()
