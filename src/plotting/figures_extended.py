from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from src.utils.config import ensure_dir


def main() -> None:
    ensure_dir("figures/extended")
    for src, dst in [
        ("figures/data_overview_umap.png", "figures/extended/extended_data_overview_umap.png"),
        ("figures/time_celltype_counts.png", "figures/extended/extended_time_celltype_counts.png"),
        ("figures/ablation_radar_or_heatmap.png", "figures/extended/extended_ablation_heatmap.png"),
    ]:
        if Path(src).exists():
            import shutil

            shutil.copyfile(src, dst)
    if Path("tables/final_metrics.csv").exists():
        data = pd.read_csv("tables/final_metrics.csv")
        fig, ax = plt.subplots(figsize=(9, 4.5), dpi=180)
        data.boxplot(column="mmd_rbf", by="model", ax=ax, rot=70, grid=False)
        ax.set_title("Extended: MMD by model")
        fig.suptitle("")
        fig.tight_layout()
        fig.savefig("figures/extended/extended_mmd_boxplot.png")
        plt.close(fig)
    print({"figures_extended": str(Path("figures/extended").resolve())})


if __name__ == "__main__":
    main()
