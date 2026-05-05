"""Plotting helpers for DevSpectrum MVP figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from devguard.io import ensure_dir


def _save(fig, path: str | Path) -> Path:
    output = Path(path)
    ensure_dir(output.parent)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)
    return output


def plot_energy_heatmap(features: pd.DataFrame, value_col: str, output: str | Path, *, title: str) -> Path:
    data = features[features["basis_method"].eq("dct")].copy()
    if data.empty:
        data = features.copy()
    pivot = data.pivot_table(index="lineage", columns="module_name", values=value_col, aggfunc="mean", fill_value=0)
    fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(pivot.columns)), max(4, 0.6 * len(pivot.index))))
    im = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=60, ha="right")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label=value_col)
    return _save(fig, output)


def plot_reconstruction_summary(summary: pd.DataFrame, output: str | Path) -> Path:
    data = summary[~summary["shuffle_time"].astype(bool)].sort_values("mean_mse")
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(data["method"].astype(str), data["mean_mse"])
    ax.set_ylabel("Mean MSE")
    ax.set_title("Missing-stage reconstruction")
    ax.tick_params(axis="x", rotation=30)
    return _save(fig, output)


def plot_fingerprint(fingerprint: pd.DataFrame, output: str | Path) -> Path:
    metrics = ["global_spectral_distance", "mean_absolute_residual", "high_frequency_burst_score"]
    data = fingerprint.set_index("perturbation_name")[metrics]
    fig, ax = plt.subplots(figsize=(6, 4))
    im = ax.imshow(data.to_numpy(), aspect="auto", cmap="magma")
    ax.set_yticks(range(len(data.index)))
    ax.set_yticklabels(data.index)
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(metrics, rotation=30, ha="right")
    ax.set_title("Perturbation spectral fingerprint")
    fig.colorbar(im, ax=ax)
    return _save(fig, output)
