from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from src.utils.config import ensure_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    out_root = Path("figures/quick_fixture") if args.quick_fixture else Path("figures")
    table_root = Path("tables/quick_fixture") if args.quick_fixture else Path("tables")
    ensure_dir(out_root / "extended")
    for src, dst in [
        (out_root / "data_overview_umap.png", out_root / "extended/extended_data_overview_umap.png"),
        (out_root / "time_celltype_counts.png", out_root / "extended/extended_time_celltype_counts.png"),
        (out_root / "ablation_radar_or_heatmap.png", out_root / "extended/extended_ablation_heatmap.png"),
    ]:
        if Path(src).exists():
            import shutil

            shutil.copyfile(src, dst)
    if (table_root / "final_metrics.csv").exists():
        data = pd.read_csv(table_root / "final_metrics.csv")
        fig, ax = plt.subplots(figsize=(9, 4.5), dpi=180)
        data.boxplot(column="mmd_rbf", by="model", ax=ax, rot=70, grid=False)
        ax.set_title("Extended: MMD by model")
        fig.suptitle("")
        fig.tight_layout()
        fig.savefig(out_root / "extended/extended_mmd_boxplot.png")
        plt.close(fig)
    print({"figures_extended": str((out_root / "extended").resolve())})


if __name__ == "__main__":
    main()
