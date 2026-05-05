from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import pandas as pd

from devguard.io import ensure_dir, write_dataframe, write_manifest
from devguard.tolerance import bootstrap_dti_ci, compute_developmental_tolerance_index, vulnerable_windows


def compute_dti(
    input_csv: str | Path = "results/devguard/perturbation_classification/cell_normality_classes.csv",
    output_dir: str | Path = "results/devguard/tolerance_index",
    n_bootstrap: int = 200,
) -> Path:
    output = ensure_dir(output_dir)
    frame = pd.read_csv(input_csv)
    dti = compute_developmental_tolerance_index(frame)
    ci = bootstrap_dti_ci(frame, n_bootstrap=n_bootstrap)
    merge_cols = ["time_point", "lineage", "perturbation_name"]
    dti = dti.merge(ci, on=merge_cols, how="left")
    dti_path = output / "developmental_tolerance_index.csv"
    vulnerable_path = output / "vulnerable_windows.csv"
    write_dataframe(dti, dti_path)
    write_dataframe(vulnerable_windows(dti), vulnerable_path)
    write_manifest(
        output / "developmental_tolerance_index_manifest.json",
        name="compute_developmental_tolerance_index",
        inputs=[str(input_csv)],
        outputs=[str(dti_path), str(vulnerable_path)],
        parameters={"n_bootstrap": n_bootstrap},
        metrics={"n_windows": int(dti.shape[0]), "min_dti": float(dti["DTI"].min())},
    )
    return dti_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute DevGuard Developmental Tolerance Index.")
    parser.add_argument("--input", default="results/devguard/perturbation_classification/cell_normality_classes.csv")
    parser.add_argument("--output-dir", default="results/devguard/tolerance_index")
    parser.add_argument("--n-bootstrap", type=int, default=200)
    args = parser.parse_args()
    compute_dti(args.input, args.output_dir, args.n_bootstrap)


if __name__ == "__main__":
    main()
