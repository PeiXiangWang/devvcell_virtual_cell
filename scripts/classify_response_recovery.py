"""Classify perturbed virtual-cell states as recovery, delay, deflection or off-manifold."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from devvcell.io import load_json, resolve_project_path  # noqa: E402
from devvcell.response_recovery import classify_from_latent_tables, summarize_response_recovery  # noqa: E402
from devvcell.tables import read_table, write_table  # noqa: E402


CLASS_COLORS = {
    "reversible_response": "#2f7d5c",
    "developmental_delay": "#c98b2c",
    "fate_deflection": "#7d4fa3",
    "off_manifold_collapse": "#b84343",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--response-config", default="config/response_recovery.json")
    parser.add_argument("--transfer-config", default="config/perturbation_transfer.json")
    return parser.parse_args()


def existing_table(path_like: str) -> Path:
    path = resolve_project_path(path_like)
    if path.exists():
        return path
    if path.suffix == ".parquet" and path.with_suffix(".csv").exists():
        return path.with_suffix(".csv")
    raise FileNotFoundError(path)


def plot_landscape(summary: pd.DataFrame, figure_path: Path) -> None:
    if summary.empty:
        return
    pivot = summary.pivot_table(
        index="stage_num",
        columns="response_recovery_class",
        values="class_fraction",
        fill_value=0.0,
        aggfunc="sum",
    )
    for column in CLASS_COLORS:
        if column not in pivot.columns:
            pivot[column] = 0.0
    pivot = pivot[list(CLASS_COLORS)]

    fig, ax = plt.subplots(figsize=(10, 5))
    bottom = None
    for column, color in CLASS_COLORS.items():
        values = pivot[column].to_numpy()
        ax.bar(pivot.index.astype(str), values, bottom=bottom, color=color, label=column)
        bottom = values if bottom is None else bottom + values
    ax.set_xlabel("Theiler stage")
    ax.set_ylabel("Class fraction")
    ax.set_title("DevVCell response-recovery landscape")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_path, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    response_cfg = load_json(args.response_config)
    transfer_cfg = load_json(args.transfer_config)

    transferred = read_table(existing_table(transfer_cfg["output"]["transferred_response_by_stage_celltype"]))
    centroids = read_table(existing_table(response_cfg["output"]["stage_celltype_centroids"]))
    classes = classify_from_latent_tables(transferred, centroids, response_cfg)
    summary = summarize_response_recovery(classes)

    classes_path = write_table(classes, response_cfg["output"]["response_recovery_classes"])
    summary_path = write_table(summary, response_cfg["output"]["stage_vulnerability_response_recovery"])
    plot_landscape(summary, resolve_project_path(response_cfg["output"]["figures"]) / "response_recovery_landscape.png")
    print(f"Wrote response-recovery classes: {classes_path}")
    print(f"Wrote stage summary: {summary_path}")


if __name__ == "__main__":
    main()
