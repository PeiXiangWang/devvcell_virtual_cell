from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401

from devspectrum.io import ensure_dir, read_table, write_dataframe, write_manifest
from devspectrum.rescue import rank_rescue_candidates


def compute_spectral_rescue_candidates(
    residuals_csv: str | Path = "results/devspectrum/chimera_projection/chimera_endpoint_spectral_residuals.csv",
    output_dir: str | Path = "results/devspectrum/rescue_candidates",
    *,
    cohort: str = "Tal1_chimera",
) -> Path:
    output = ensure_dir(output_dir)
    residuals = read_table(residuals_csv)
    candidates, report = rank_rescue_candidates(residuals, cohort=cohort)
    candidate_path = output / "spectral_rescue_candidates.csv"
    report_path = output / "tal1_spectral_rescue_report.md"
    write_dataframe(candidates, candidate_path)
    report_path.write_text(report, encoding="utf-8")
    write_manifest(
        output / "spectral_rescue_candidates_manifest.json",
        name="compute_spectral_rescue_candidates",
        inputs=[str(residuals_csv)],
        outputs=[str(candidate_path), str(report_path)],
        parameters={"cohort": cohort},
        metrics={"n_candidates": int(candidates.shape[0])},
    )
    return candidate_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute in silico spectral rescue candidates.")
    parser.add_argument("--residuals", default="results/devspectrum/chimera_projection/chimera_endpoint_spectral_residuals.csv")
    parser.add_argument("--output-dir", default="results/devspectrum/rescue_candidates")
    parser.add_argument("--cohort", default="Tal1_chimera")
    args = parser.parse_args()
    compute_spectral_rescue_candidates(args.residuals, args.output_dir, cohort=args.cohort)


if __name__ == "__main__":
    main()
