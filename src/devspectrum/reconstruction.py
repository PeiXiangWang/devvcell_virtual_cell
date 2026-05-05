"""Missing-stage reconstruction benchmarks."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from devspectrum.basis import fit_basis


def _safe_corr(left: np.ndarray, right: np.ndarray, fn) -> float:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if left.size < 2 or np.nanstd(left) == 0 or np.nanstd(right) == 0:
        return float("nan")
    return float(fn(left, right).statistic)


def _mean_predict(train_times: np.ndarray, train_values: np.ndarray, query_times: np.ndarray) -> np.ndarray:
    return np.repeat(float(np.mean(train_values)), len(query_times))


def _linear_predict(train_times: np.ndarray, train_values: np.ndarray, query_times: np.ndarray) -> np.ndarray:
    order = np.argsort(train_times)
    return np.interp(query_times, train_times[order], train_values[order])


def leave_one_timepoint_predictions(
    timeseries: pd.DataFrame,
    *,
    methods: list[str] | tuple[str, ...] = ("mean", "linear", "spline", "dct", "wavelet", "dct_wavelet"),
    basis_config: dict | None = None,
    group_cols: list[str] | None = None,
    shuffle_time: bool = False,
    seed: int = 42,
) -> pd.DataFrame:
    """Predict each held-out time point for each lineage-module trajectory."""

    basis_config = basis_config or {}
    group_cols = group_cols or ["dataset_id", "condition", "lineage", "module_name", "feature_type"]
    data = timeseries.copy()
    if "n_samples" not in data.columns:
        data["n_samples"] = 1
    data = (
        data.groupby(group_cols + ["time_numeric", "time_point"], dropna=False, observed=True)
        .agg(feature_value=("feature_value", "mean"), n_cells=("n_cells", "sum"), n_samples=("n_samples", "max"))
        .reset_index()
    )
    rng = np.random.default_rng(seed)
    rows = []
    for keys, group in data.groupby(group_cols, dropna=False, observed=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        group = group.sort_values("time_numeric")
        times = group["time_numeric"].to_numpy(dtype=float)
        values = group["feature_value"].to_numpy(dtype=float)
        if np.unique(times).size < 3:
            continue
        fit_times = rng.permutation(times) if shuffle_time else times
        for holdout_idx, holdout_time in enumerate(times):
            train_mask = np.arange(times.size) != holdout_idx
            train_times = fit_times[train_mask]
            train_values = values[train_mask]
            query = np.asarray([fit_times[holdout_idx]], dtype=float)
            observed = float(values[holdout_idx])
            for method in methods:
                if method == "mean":
                    predicted = _mean_predict(train_times, train_values, query)
                elif method == "linear":
                    predicted = _linear_predict(train_times, train_values, query)
                elif method == "dct_wavelet":
                    dct_fit = fit_basis(train_times, train_values, "dct", basis_config.get("dct", {}))
                    wavelet_fit = fit_basis(train_times, train_values, "wavelet", basis_config.get("wavelet", {}))
                    predicted = (dct_fit.predict(query) + wavelet_fit.predict(query)) / 2.0
                else:
                    config = basis_config.get(method, {})
                    predicted = fit_basis(train_times, train_values, method, config).predict(query)
                row = {column: value for column, value in zip(group_cols, keys)}
                row.update(
                    {
                        "holdout_time_numeric": float(holdout_time),
                        "holdout_time_point": str(group.iloc[holdout_idx]["time_point"]),
                        "method": method,
                        "observed": observed,
                        "predicted": float(predicted[0]),
                        "residual": float(observed - predicted[0]),
                        "squared_error": float((observed - predicted[0]) ** 2),
                        "absolute_error": float(abs(observed - predicted[0])),
                        "shuffle_time": bool(shuffle_time),
                    }
                )
                rows.append(row)
    return pd.DataFrame(rows)


def summarize_reconstruction(predictions: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if predictions.empty:
        return pd.DataFrame(), pd.DataFrame()
    group_cols = ["method", "shuffle_time", "lineage", "module_name", "feature_type"]
    metrics = (
        predictions.groupby(group_cols, dropna=False, observed=True)
        .agg(
            mse=("squared_error", "mean"),
            mae=("absolute_error", "mean"),
            n_predictions=("observed", "size"),
            bias=("residual", "mean"),
        )
        .reset_index()
    )
    corr_rows = []
    for keys, group in predictions.groupby(group_cols, dropna=False, observed=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = {column: value for column, value in zip(group_cols, keys)}
        observed = group["observed"].to_numpy(dtype=float)
        predicted = group["predicted"].to_numpy(dtype=float)
        row["pearson"] = _safe_corr(observed, predicted, pearsonr)
        row["spearman"] = _safe_corr(observed, predicted, spearmanr)
        corr_rows.append(row)
    metrics = metrics.merge(pd.DataFrame(corr_rows), on=group_cols, how="left")
    overall = (
        metrics.groupby(["method", "shuffle_time"], dropna=False, observed=True)
        .agg(
            mean_mse=("mse", "mean"),
            median_mse=("mse", "median"),
            mean_mae=("mae", "mean"),
            median_mae=("mae", "median"),
            mean_pearson=("pearson", "mean"),
            mean_spearman=("spearman", "mean"),
            n_groups=("module_name", "size"),
        )
        .reset_index()
        .sort_values(["shuffle_time", "mean_mse"])
    )
    return metrics.sort_values(["method", "lineage", "module_name"]).reset_index(drop=True), overall.reset_index(drop=True)
