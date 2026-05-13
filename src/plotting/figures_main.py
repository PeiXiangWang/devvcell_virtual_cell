from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.utils.config import ensure_dir


def _copy_panel(src: str, dst: str) -> None:
    path = Path(src)
    if path.exists():
        import shutil

        shutil.copyfile(path, dst)


def concept_figure(out_root: Path) -> None:
    ensure_dir(out_root / "main")
    fig, ax = plt.subplots(figsize=(10, 4.8), dpi=180)
    ax.axis("off")
    boxes = [
        (0.06, 0.62, "Time-series\nsingle-cell data"),
        (0.28, 0.62, "moscot/WOT\nOT teacher"),
        (0.50, 0.62, "Finite-cell\nagent simulator"),
        (0.72, 0.62, "Held-out\nreconstruction"),
        (0.50, 0.20, "local rules\nbirth/death\ndiffusion\nCCI\nphenomenological memory"),
        (0.72, 0.20, "in silico\ncontrol-layer\nperturbations"),
    ]
    for x, y, text in boxes:
        ax.add_patch(plt.Rectangle((x, y), 0.17, 0.18, facecolor="#eef4f8", edgecolor="#315b7d", linewidth=1.2))
        ax.text(x + 0.085, y + 0.09, text, ha="center", va="center", fontsize=9)
    arrows = [((0.23, 0.71), (0.28, 0.71)), ((0.45, 0.71), (0.50, 0.71)), ((0.67, 0.71), (0.72, 0.71)), ((0.585, 0.62), (0.585, 0.38)), ((0.67, 0.29), (0.72, 0.29))]
    for start, end in arrows:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", color="#315b7d", lw=1.4))
    ax.text(0.05, 0.93, "SwarmLineage-OT", fontsize=18, weight="bold")
    ax.text(0.05, 0.88, "OT-inferred pseudo-lineage supervises an executable swarm virtual-cell population", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_root / "main/figure1_concept.png")
    plt.close(fig)


def metrics_figure(out_root: Path, table_root: Path) -> None:
    ensure_dir(out_root / "main")
    path = table_root / "final_metrics.csv"
    if not path.exists():
        return
    data = pd.read_csv(path)
    order = data.groupby("model")["sinkhorn"].mean().sort_values().index.tolist()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), dpi=180)
    for ax, metric in zip(axes, ["sinkhorn", "celltype_composition_rmse"]):
        means = data.groupby("model")[metric].mean().reindex(order)
        std = data.groupby("model")[metric].std().reindex(order)
        ax.bar(range(len(order)), means, yerr=std, color="#6d9dc5", edgecolor="#25445c", linewidth=0.4)
        for i, model in enumerate(order):
            vals = data.loc[data["model"] == model, metric]
            ax.scatter(np.full(vals.shape, i), vals, color="#1b2935", s=8, zorder=3)
        ax.set_xticks(range(len(order)), order, rotation=65, ha="right", fontsize=6)
        ax.set_ylabel(metric)
        ax.set_title(metric)
    fig.tight_layout()
    fig.savefig(out_root / "main/figure3_heldout_reconstruction.png")
    plt.close(fig)


def mechanism_figure(out_root: Path, table_root: Path) -> None:
    ensure_dir(out_root / "main")
    growth = out_root / "ot_growth_map.png"
    if growth.exists():
        _copy_panel(str(growth), str(out_root / "main/figure4_growth_diffusion.png"))
    lr_path = table_root / "lr_knockout_predictions.csv"
    if lr_path.exists():
        data = pd.read_csv(lr_path)
        fig, ax = plt.subplots(figsize=(6.5, 4), dpi=180)
        labels = [f"{r.ligand}-{r.receptor}" for r in data.itertuples(index=False)]
        ax.bar(range(len(data)), data["predicted_entropy_change"], color="#8fbf9f", edgecolor="#244", linewidth=0.4)
        ax.set_xticks(range(len(data)), labels, rotation=45, ha="right", fontsize=7)
        ax.set_ylabel("predicted entropy change")
        ax.set_title("Exploratory LR control-layer perturbation")
        fig.tight_layout()
        fig.savefig(out_root / "main/figure5_cci_perturbation.png")
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    out_root = Path("figures/quick_fixture") if args.quick_fixture else Path("figures")
    table_root = Path("tables/quick_fixture") if args.quick_fixture else Path("tables")
    ensure_dir(out_root / "main")
    concept_figure(out_root)
    _copy_panel(str(out_root / "ot_lineage_graph.png"), str(out_root / "main/figure2_ot_teacher_lineage_graph.png"))
    _copy_panel(str(out_root / "ot_fate_umap.png"), str(out_root / "main/figure2_ot_fate_umap.png"))
    metrics_figure(out_root, table_root)
    mechanism_figure(out_root, table_root)
    _copy_panel(str(out_root / "ablation_barplots.png"), str(out_root / "main/figure6_ablation_summary.png"))
    print({"figures_main": str((out_root / "main").resolve())})


if __name__ == "__main__":
    main()
