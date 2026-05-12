from __future__ import annotations

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


def concept_figure() -> None:
    ensure_dir("figures/main")
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
    fig.savefig("figures/main/figure1_concept.png")
    plt.close(fig)


def metrics_figure() -> None:
    ensure_dir("figures/main")
    path = Path("tables/final_metrics.csv")
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
    fig.savefig("figures/main/figure3_heldout_reconstruction.png")
    plt.close(fig)


def mechanism_figure() -> None:
    ensure_dir("figures/main")
    if Path("figures/ot_growth_map.png").exists():
        _copy_panel("figures/ot_growth_map.png", "figures/main/figure4_growth_diffusion.png")
    if Path("tables/lr_knockout_predictions.csv").exists():
        data = pd.read_csv("tables/lr_knockout_predictions.csv")
        fig, ax = plt.subplots(figsize=(6.5, 4), dpi=180)
        labels = [f"{r.ligand}-{r.receptor}" for r in data.itertuples(index=False)]
        ax.bar(range(len(data)), data["predicted_entropy_change"], color="#8fbf9f", edgecolor="#244", linewidth=0.4)
        ax.set_xticks(range(len(data)), labels, rotation=45, ha="right", fontsize=7)
        ax.set_ylabel("predicted entropy change")
        ax.set_title("Exploratory LR control-layer perturbation")
        fig.tight_layout()
        fig.savefig("figures/main/figure5_cci_perturbation.png")
        plt.close(fig)


def main() -> None:
    ensure_dir("figures/main")
    concept_figure()
    _copy_panel("figures/ot_lineage_graph.png", "figures/main/figure2_ot_teacher_lineage_graph.png")
    _copy_panel("figures/ot_fate_umap.png", "figures/main/figure2_ot_fate_umap.png")
    metrics_figure()
    mechanism_figure()
    _copy_panel("figures/ablation_barplots.png", "figures/main/figure6_ablation_summary.png")
    print({"figures_main": str(Path("figures/main").resolve())})


if __name__ == "__main__":
    main()
