from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import joblib
import pandas as pd

from devspectrum.basis import basis_feature_row
from devspectrum.io import ensure_dir, load_config, read_table, write_dataframe, write_manifest


def fit_spectral_basis(
    timeseries_csv: str | Path,
    config_path: str | Path,
    *,
    output_dir: str | Path = "results/devspectrum/spectral_fits",
) -> Path:
    config = load_config(config_path)
    timeseries = read_table(timeseries_csv)
    output = ensure_dir(output_dir)
    methods = config.get("basis_methods", ["dct", "wavelet", "spline"])
    group_cols = ["dataset_id", "condition", "lineage", "module_name", "feature_type"]
    aggregate = (
        timeseries.groupby(group_cols + ["time_numeric", "time_point"], dropna=False, observed=True)
        .agg(feature_value=("feature_value", "mean"), n_cells=("n_cells", "sum"))
        .reset_index()
    )
    rows = []
    quality_rows = []
    models = {}
    for keys, group in aggregate.groupby(group_cols, dropna=False, observed=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        group = group.sort_values("time_numeric")
        if group["time_numeric"].nunique() < 3:
            continue
        times = group["time_numeric"].to_numpy(dtype=float)
        values = group["feature_value"].to_numpy(dtype=float)
        key_dict = {column: value for column, value in zip(group_cols, keys)}
        for method in methods:
            row, fit = basis_feature_row(times, values, method, config.get(method, {}))
            row.update(key_dict)
            row["n_time_points"] = int(group["time_numeric"].nunique())
            row["n_cells"] = int(group["n_cells"].sum())
            rows.append(row)
            quality_rows.append(
                {
                    **key_dict,
                    "basis_method": method,
                    "reconstruction_error": row.get("reconstruction_error", 0.0),
                    "n_time_points": int(group["time_numeric"].nunique()),
                }
            )
            models[tuple([*keys, method])] = fit
    features = pd.DataFrame(rows)
    coeff_cols = [column for column in features.columns if column.startswith("c")]
    coefficients = features[group_cols + ["basis_method", *coeff_cols]].copy() if not features.empty else pd.DataFrame()
    quality = pd.DataFrame(quality_rows)
    feature_path = output / "spectral_features.csv"
    coefficient_path = output / "spectral_coefficients.csv"
    quality_path = output / "spectral_fit_quality.csv"
    model_path = output / "control_spectral_model.joblib"
    write_dataframe(features, feature_path)
    write_dataframe(coefficients, coefficient_path)
    write_dataframe(quality, quality_path)
    joblib.dump({"models": models, "config": config, "timeseries": str(timeseries_csv)}, model_path)
    write_manifest(
        output / "spectral_fits_manifest.json",
        name="fit_spectral_basis",
        inputs=[str(timeseries_csv), str(config_path)],
        outputs=[str(feature_path), str(coefficient_path), str(quality_path), str(model_path)],
        parameters=config,
        metrics={"n_feature_rows": int(features.shape[0]), "n_models": int(len(models))},
    )
    return feature_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Fit DevSpectrum spectral basis models.")
    parser.add_argument("--timeseries", default="results/devspectrum/timeseries/stage_lineage_module_timeseries.csv")
    parser.add_argument("--config", default="config/devspectrum/spectral_basis.json")
    parser.add_argument("--output-dir", default="results/devspectrum/spectral_fits")
    args = parser.parse_args()
    fit_spectral_basis(args.timeseries, args.config, output_dir=args.output_dir)


if __name__ == "__main__":
    main()
