from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import pandas as pd

from devspectrum.io import ensure_dir, load_config, write_dataframe, write_manifest
from devspectrum.perturbation_compare import endpoint_residuals, perturbation_fingerprint, score_classified_endpoint


def project_chimera_endpoint_spectra(config_path: str | Path) -> Path:
    config = load_config(config_path)
    output = ensure_dir(config.get("output_dir", "results/devspectrum/chimera_projection"))
    module_registry = config.get("modules", {}).get("module_registry", "config/devspectrum/module_registry_mouse.json")
    control_scores, coverage = score_classified_endpoint(
        config["control_input_h5ad"],
        config["devguard_outputs"]["control_classes"],
        cohort="Control",
        module_registry=module_registry,
    )
    score_frames = []
    for cohort, h5ad_path in config.get("chimera_inputs", {}).items():
        class_path = config["devguard_outputs"][f"{cohort}_classes"]
        label = "Tal1_chimera" if cohort.lower() == "tal1" else ("T_chimera" if cohort.lower() == "t" else cohort)
        scores, _ = score_classified_endpoint(h5ad_path, class_path, cohort=label, module_registry=module_registry)
        score_frames.append(scores)
    perturbation_scores = pd.concat(score_frames, axis=0, ignore_index=True, sort=False)
    residuals, summary = endpoint_residuals(
        control_scores,
        perturbation_scores,
        min_cells_per_group=int(config.get("min_cells_per_group", 10)),
    )
    fingerprint = perturbation_fingerprint(residuals)
    residual_path = output / "chimera_endpoint_spectral_residuals.csv"
    summary_path = output / "tal1_t_spectral_residual_summary.csv"
    fingerprint_path = output.parent / "perturbation_spectra" / "perturbation_spectral_fingerprint.csv"
    coverage_path = output / "chimera_endpoint_module_gene_coverage.csv"
    write_dataframe(residuals, residual_path)
    write_dataframe(summary, summary_path)
    write_dataframe(fingerprint, fingerprint_path)
    write_dataframe(coverage, coverage_path)
    write_manifest(
        output / "chimera_projection_manifest.json",
        name="project_chimera_endpoint_spectra",
        inputs=[config["control_input_h5ad"], *config.get("chimera_inputs", {}).values()],
        outputs=[str(residual_path), str(summary_path), str(fingerprint_path), str(coverage_path)],
        parameters=config,
        metrics={"n_residual_rows": int(residuals.shape[0]), "n_fingerprint_rows": int(fingerprint.shape[0])},
    )
    return residual_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Project Tal1/T chimera endpoints into DevSpectrum module-residual space.")
    parser.add_argument("--config", default="config/devspectrum/devspectrum_chimera_projection.json")
    args = parser.parse_args()
    project_chimera_endpoint_spectra(args.config)


if __name__ == "__main__":
    main()
