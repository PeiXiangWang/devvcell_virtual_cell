from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401

from devspectrum.io import ensure_dir, load_config, read_table, write_manifest
from devspectrum.plotting import plot_energy_heatmap, plot_fingerprint, plot_reconstruction_summary


def build_devspectrum_figures(config_path: str | Path = "config/devspectrum/figure_plan.json") -> list[Path]:
    config = load_config(config_path)
    output = ensure_dir(config.get("output_dir", "results/devspectrum/figures"))
    spectral = read_table(config.get("spectral_features", "results/devspectrum/spectral_fits/spectral_features.csv"))
    reconstruction = read_table(config.get("reconstruction_summary", "results/devspectrum/missing_stage_reconstruction/reconstruction_summary.csv"))
    fingerprint = read_table(config.get("fingerprint", "results/devspectrum/perturbation_spectra/perturbation_spectral_fingerprint.csv"))
    outputs = [
        plot_energy_heatmap(spectral, "low_frequency_energy", output / "figure2_low_frequency_energy.png", title="Low-frequency developmental energy"),
        plot_energy_heatmap(spectral, "high_frequency_energy", output / "figure2_high_frequency_energy.png", title="High-frequency developmental energy"),
        plot_energy_heatmap(spectral, "spectral_entropy", output / "figure2_spectral_entropy.png", title="Spectral entropy"),
        plot_reconstruction_summary(reconstruction, output / "figure3_missing_stage_reconstruction.png"),
        plot_fingerprint(fingerprint, output / "figure4_perturbation_spectral_fingerprint.png"),
    ]
    write_manifest(
        output / "figure_manifest.json",
        name="build_devspectrum_figures",
        inputs=[str(config_path)],
        outputs=[str(path) for path in outputs],
        parameters=config,
        metrics={"n_figures": len(outputs)},
    )
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description="Build DevSpectrum MVP figures.")
    parser.add_argument("--config", default="config/devspectrum/figure_plan.json")
    args = parser.parse_args()
    build_devspectrum_figures(args.config)


if __name__ == "__main__":
    main()
