from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import pandas as pd

from devspectrum.io import ensure_dir, read_table, write_dataframe, write_manifest


def compare_control_perturbation_spectra(
    spectral_features: str | Path = "results/devspectrum/spectral_fits/spectral_features.csv",
    output_dir: str | Path = "results/devspectrum/perturbation_spectra",
) -> Path:
    output = ensure_dir(output_dir)
    features = read_table(spectral_features)
    control = features[features["condition"].astype(str).str.lower().eq("control")]
    rows = []
    for _, row in features.iterrows():
        if str(row.get("condition", "")).lower() == "control":
            continue
        match = control[
            control["lineage"].astype(str).eq(str(row["lineage"]))
            & control["module_name"].astype(str).eq(str(row["module_name"]))
            & control["basis_method"].astype(str).eq(str(row["basis_method"]))
        ]
        if match.empty:
            continue
        baseline = match.iloc[0]
        rows.append(
            {
                "condition": row["condition"],
                "lineage": row["lineage"],
                "module_name": row["module_name"],
                "basis_method": row["basis_method"],
                "low_frequency_delta": row.get("low_frequency_energy", 0) - baseline.get("low_frequency_energy", 0),
                "high_frequency_delta": row.get("high_frequency_energy", 0) - baseline.get("high_frequency_energy", 0),
                "entropy_delta": row.get("spectral_entropy", 0) - baseline.get("spectral_entropy", 0),
            }
        )
    out = pd.DataFrame(rows)
    path = output / "control_perturbation_spectral_delta.csv"
    write_dataframe(out, path)
    write_manifest(
        output / "control_perturbation_spectral_delta_manifest.json",
        name="compare_control_perturbation_spectra",
        inputs=[str(spectral_features)],
        outputs=[str(path)],
        parameters={},
        metrics={"n_rows": int(out.shape[0])},
    )
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare control and perturbation spectra when multi-timepoint perturbation data are available.")
    parser.add_argument("--spectral-features", default="results/devspectrum/spectral_fits/spectral_features.csv")
    parser.add_argument("--output-dir", default="results/devspectrum/perturbation_spectra")
    args = parser.parse_args()
    compare_control_perturbation_spectra(args.spectral_features, args.output_dir)


if __name__ == "__main__":
    main()
