from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401

from build_devspectrum_figures import build_devspectrum_figures
from build_stage_lineage_module_timeseries import build_timeseries_from_config
from compute_spectral_rescue_candidates import compute_spectral_rescue_candidates
from fit_spectral_basis import fit_spectral_basis
from link_spectrum_to_devguard import link_spectrum_to_devguard
from project_chimera_endpoint_spectra import project_chimera_endpoint_spectra
from reconstruct_missing_stages import reconstruct_missing_stages


def run_pipeline(mode: str = "quick") -> None:
    control_config = "config/devspectrum/devspectrum_gse212050_control.json"
    basis_config = "config/devspectrum/spectral_basis.json"
    reconstruction_config = "config/devspectrum/missing_stage_reconstruction.json"
    chimera_config = "config/devspectrum/devspectrum_chimera_projection.json"
    quick = mode == "quick"
    timeseries = build_timeseries_from_config(control_config, quick=quick)
    fit_spectral_basis(timeseries, basis_config)
    reconstruct_missing_stages(reconstruction_config)
    if mode in {"main", "full"}:
        project_chimera_endpoint_spectra(chimera_config)
        link_spectrum_to_devguard(chimera_config)
        compute_spectral_rescue_candidates()
        build_devspectrum_figures()


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the DevSpectrum MVP pipeline.")
    parser.add_argument("--mode", choices=["quick", "main", "full"], default="quick")
    args = parser.parse_args()
    run_pipeline(args.mode)


if __name__ == "__main__":
    main()
