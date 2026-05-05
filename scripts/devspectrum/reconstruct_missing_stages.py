from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401

from devspectrum.io import ensure_dir, load_config, read_table, write_dataframe, write_manifest
from devspectrum.reconstruction import leave_one_timepoint_predictions, summarize_reconstruction


def reconstruct_missing_stages(config_path: str | Path) -> Path:
    config = load_config(config_path)
    output = ensure_dir(config.get("output_dir", "results/devspectrum/missing_stage_reconstruction"))
    timeseries_csv = config.get("timeseries", "results/devspectrum/timeseries/stage_lineage_module_timeseries.csv")
    basis_config = config.get("basis_config", {})
    timeseries = read_table(timeseries_csv)
    methods = config.get("baseline_methods", ["mean", "linear", "spline"]) + config.get("spectral_methods", ["dct", "wavelet", "dct_wavelet"])
    predictions = leave_one_timepoint_predictions(timeseries, methods=methods, basis_config=basis_config)
    shuffle_predictions = leave_one_timepoint_predictions(
        timeseries,
        methods=methods,
        basis_config=basis_config,
        shuffle_time=True,
        seed=int(config.get("bootstrap", {}).get("seed", 42)),
    )
    predictions = predictions.assign(control_type="observed_time")
    shuffle_predictions = shuffle_predictions.assign(control_type="shuffle_time")
    import pandas as pd

    combined = pd.concat([predictions, shuffle_predictions], axis=0, ignore_index=True)
    metrics, summary = summarize_reconstruction(combined)
    predictions_path = output / "reconstruction_predictions.csv"
    metrics_path = output / "reconstruction_metrics.csv"
    summary_path = output / "reconstruction_summary.csv"
    write_dataframe(combined, predictions_path)
    write_dataframe(metrics, metrics_path)
    write_dataframe(summary, summary_path)
    write_manifest(
        output / "missing_stage_reconstruction_manifest.json",
        name="reconstruct_missing_stages",
        inputs=[str(timeseries_csv), str(config_path)],
        outputs=[str(predictions_path), str(metrics_path), str(summary_path)],
        parameters=config,
        metrics={"n_predictions": int(combined.shape[0]), "n_metric_rows": int(metrics.shape[0])},
    )
    return summary_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Run DevSpectrum missing-stage reconstruction benchmark.")
    parser.add_argument("--config", default="config/devspectrum/missing_stage_reconstruction.json")
    args = parser.parse_args()
    reconstruct_missing_stages(args.config)


if __name__ == "__main__":
    main()
