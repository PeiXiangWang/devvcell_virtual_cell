"""Figure draft builders for DevGuard."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from devguard.io import ensure_dir
from devguard.tolerance import CLASS_COLUMNS


def plot_method_schematic(output_path: str | Path) -> Path:
    output = Path(output_path)
    ensure_dir(output.parent)
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.axis("off")
    boxes = [
        (0.05, 0.58, "Control mouse\nstage-lineage cells"),
        (0.29, 0.58, "Train embedding\nand reference"),
        (0.53, 0.58, "Calibrate\nconformal scores"),
        (0.77, 0.58, "Classify\nperturbed cells"),
        (0.53, 0.15, "within normal | delay | acceleration\nfate deviation | abnormal"),
    ]
    for x, y, text in boxes:
        ax.add_patch(plt.Rectangle((x, y), 0.18, 0.22, fill=False, linewidth=1.8))
        ax.text(x + 0.09, y + 0.11, text, ha="center", va="center", fontsize=10)
    for start, end in [(0.23, 0.29), (0.47, 0.53), (0.71, 0.77)]:
        ax.annotate("", xy=(end, 0.69), xytext=(start, 0.69), arrowprops={"arrowstyle": "->", "lw": 1.8})
    ax.annotate("", xy=(0.62, 0.37), xytext=(0.86, 0.58), arrowprops={"arrowstyle": "->", "lw": 1.8})
    ax.set_title("DevGuard conformal developmental normality", fontsize=14, pad=12)
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)
    return output


def plot_class_summary(summary: pd.DataFrame, output_path: str | Path) -> Path:
    output = Path(output_path)
    ensure_dir(output.parent)
    if summary.empty:
        raise ValueError("Class summary is empty.")
    label_cols = [column for column in ["perturbation_name", "time_point", "lineage"] if column in summary.columns]
    work = summary.copy()
    work["label"] = work[label_cols].astype(str).agg(" | ".join, axis=1)
    pivot = work.pivot_table(index="label", columns="normality_class", values="fraction", fill_value=0)
    for column in CLASS_COLUMNS:
        if column not in pivot.columns:
            pivot[column] = 0
    pivot = pivot[CLASS_COLUMNS]
    fig, ax = plt.subplots(figsize=(max(8, 0.7 * len(pivot)), 4.8))
    bottom = np.zeros(len(pivot))
    colors = ["#2a9d8f", "#e9c46a", "#f4a261", "#457b9d", "#e76f51"]
    for column, color in zip(CLASS_COLUMNS, colors):
        ax.bar(np.arange(len(pivot)), pivot[column].to_numpy(), bottom=bottom, label=column, color=color)
        bottom += pivot[column].to_numpy()
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of perturbed cells")
    ax.set_xticks(np.arange(len(pivot)))
    ax.set_xticklabels(pivot.index, rotation=45, ha="right")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)
    return output


def plot_dti_heatmap(dti: pd.DataFrame, output_path: str | Path) -> Path:
    output = Path(output_path)
    ensure_dir(output.parent)
    if dti.empty:
        raise ValueError("DTI table is empty.")
    row_col = "lineage" if "lineage" in dti.columns else dti.columns[0]
    col_col = "time_point" if "time_point" in dti.columns else dti.columns[1]
    pivot = dti.pivot_table(index=row_col, columns=col_col, values="DTI", aggfunc="mean")
    fig, ax = plt.subplots(figsize=(7, 4.8))
    image = ax.imshow(pivot.to_numpy(), vmin=-1, vmax=1, cmap="RdYlGn", aspect="auto")
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns.astype(str))
    ax.set_yticks(np.arange(pivot.shape[0]))
    ax.set_yticklabels(pivot.index.astype(str))
    ax.set_xlabel(col_col)
    ax.set_ylabel(row_col)
    fig.colorbar(image, ax=ax, label="DTI")
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)
    return output
