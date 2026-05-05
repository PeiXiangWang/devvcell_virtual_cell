from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import pandas as pd

from devspectrum.devguard_link import correlate_dti_with_spectral_residuals, failure_mode_spectral_signature, write_link_summary
from devspectrum.io import ensure_dir, load_config, read_table, write_dataframe, write_manifest


def link_spectrum_to_devguard(config_path: str | Path) -> Path:
    config = load_config(config_path)
    output = ensure_dir(config.get("devguard_link_output_dir", "results/devspectrum/devguard_link"))
    residual_csv = config.get("residuals_csv", "results/devspectrum/chimera_projection/chimera_endpoint_spectral_residuals.csv")
    summary_csv = config.get("residual_summary_csv", "results/devspectrum/chimera_projection/tal1_t_spectral_residual_summary.csv")
    residuals = read_table(residual_csv)
    residual_summary = read_table(summary_csv)
    signature = failure_mode_spectral_signature(residuals)
    dti_tables = {}
    for cohort, path in config.get("devguard_outputs", {}).get("dti", {}).items():
        dti = read_table(path)
        dti_tables["Tal1_chimera" if cohort.lower() == "tal1" else "T_chimera"] = dti
    correlation = correlate_dti_with_spectral_residuals(residual_summary, dti_tables)
    summary_md = write_link_summary(correlation, signature)
    correlation_path = output / "devguard_spectral_correlation.csv"
    signature_path = output / "failure_mode_spectral_signature.csv"
    summary_path = output / "devguard_spectrum_link_summary.md"
    write_dataframe(correlation, correlation_path)
    write_dataframe(signature, signature_path)
    summary_path.write_text(summary_md, encoding="utf-8")
    write_manifest(
        output / "devguard_link_manifest.json",
        name="link_spectrum_to_devguard",
        inputs=[str(residual_csv), str(summary_csv), *config.get("devguard_outputs", {}).get("dti", {}).values()],
        outputs=[str(correlation_path), str(signature_path), str(summary_path)],
        parameters=config,
        metrics={"n_correlation_rows": int(correlation.shape[0]), "n_signature_rows": int(signature.shape[0])},
    )
    return correlation_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Link DevSpectrum residuals to DevGuard DTI and failure classes.")
    parser.add_argument("--config", default="config/devspectrum/devspectrum_chimera_projection.json")
    args = parser.parse_args()
    link_spectrum_to_devguard(args.config)


if __name__ == "__main__":
    main()
