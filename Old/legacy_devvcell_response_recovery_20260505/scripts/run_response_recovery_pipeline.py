"""Run the DevVCell response-recovery pipeline.

Modes:
  quick: RDEG-derived proxy classes for regression and manuscript planning.
  main: embryo manifold + external response dictionary + transfer + classes.
  full: main mode plus niche context. External validation datasets are tracked
        in config/external_datasets.json and should be downloaded into data/.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from devvcell.io import load_json, resolve_project_path  # noqa: E402
from devvcell.response_recovery import classify_from_latent_tables, summarize_response_recovery  # noqa: E402
from devvcell.tables import minmax, write_table  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["quick", "main", "full"], default="quick")
    parser.add_argument("--response-config", default="config/response_recovery.json")
    parser.add_argument("--transfer-config", default="config/perturbation_transfer.json")
    parser.add_argument("--niche-config", default="config/niche_context.json")
    parser.add_argument("--allow-missing-external", action="store_true")
    return parser.parse_args()


def run_script(script: str, *args: str) -> None:
    command = [sys.executable, str(PROJECT_ROOT / "scripts" / script), *args]
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)


def ensure_lite_outputs(cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    tables_dir = resolve_project_path(cfg["input"]["devvcell_lite_tables_dir"])
    stage_path = tables_dir / "stage_vulnerability.csv"
    priority_path = tables_dir / "perturbation_priority.csv"
    if not stage_path.exists() or not priority_path.exists():
        run_script("devvcell_lite.py")
    return pd.read_csv(stage_path), pd.read_csv(priority_path)


def proxy_centroids(stage_df: pd.DataFrame) -> pd.DataFrame:
    stage = stage_df.copy()
    for column in ["vulnerability_score", "baseline_dev_step_norm", "mean_outlier_ratio", "next_step_mse"]:
        if column not in stage.columns:
            stage[column] = 0.0
        stage[f"{column}_scaled"] = minmax(stage[column])

    rows: list[dict[str, object]] = []
    cell_offsets = {
        "neural_response_axis": 0.0,
        "mesoderm_alternative_axis": 1.0,
        "pgc_gonad_axis": 2.0,
    }
    for _, row in stage.iterrows():
        for cell_type, offset in cell_offsets.items():
            rows.append(
                {
                    "devvcell_system": "rdeg_proxy",
                    "stage": row["stage"],
                    "stage_num": int(row["stage_num"]),
                    "cell_type": cell_type,
                    "n_cells": int(row.get("total_cells", 1)),
                    "latent_01": float(minmax(stage["stage_num"]).loc[row.name]),
                    "latent_02": float(row["vulnerability_score_scaled"]),
                    "latent_03": float(row["baseline_dev_step_norm_scaled"]),
                    "latent_04": float(row["mean_outlier_ratio_scaled"]),
                    "latent_05": offset,
                    "latent_06": float(row["next_step_mse_scaled"]),
                }
            )
    return pd.DataFrame(rows)


def proxy_transferred(stage_df: pd.DataFrame, priority: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    proxy_cfg = cfg["quick_proxy"]
    top = priority.head(int(proxy_cfg["top_perturbations"])).copy()
    stage = stage_df.copy()
    stage["vulnerability_scaled"] = minmax(stage["vulnerability_score"])
    stage["outlier_scaled"] = minmax(stage.get("mean_outlier_ratio", pd.Series(np.zeros(len(stage)))))

    rows: list[dict[str, object]] = []
    for _, srow in stage.iterrows():
        in_window = cfg["stage_window_of_interest"]["min_stage"] <= int(srow["stage_num"]) <= cfg["stage_window_of_interest"]["max_stage"]
        for _, prow in top.iterrows():
            priority_score = float(prow.get("devvcell_priority_score", 0.0))
            response = float(proxy_cfg["response_scale"]) * (float(srow["vulnerability_scaled"]) + priority_score)
            delay = -float(proxy_cfg["delay_weight"]) * float(prow.get("developmental_delay_score", 0.0)) * (1.3 if in_window else 0.7)
            fate = float(proxy_cfg["fate_weight"]) * float(prow.get("fate_displacement_proxy", 0.0))
            off = float(proxy_cfg["off_manifold_weight"]) * float(srow["outlier_scaled"]) * (1.0 + priority_score)
            rescue = float(prow.get("best_rescue_fraction", 0.0))
            recovery = -float(proxy_cfg["rescue_weight"]) * rescue
            rows.append(
                {
                    "perturbation": prow.get("tf_name", "unknown"),
                    "external_cell_type": "RDEG-derived",
                    "stage": srow["stage"],
                    "stage_num": int(srow["stage_num"]),
                    "cell_type": "neural_response_axis",
                    "devvcell_system": "rdeg_proxy",
                    "transfer_confidence": 1.0,
                    "response_norm": response,
                    "source_response_norm": response,
                    "transfer_method": "rdeg_proxy_quick",
                    "response_latent_01": delay,
                    "response_latent_02": response,
                    "response_latent_03": recovery,
                    "response_latent_04": off,
                    "response_latent_05": fate,
                    "response_latent_06": off + response * 0.25,
                }
            )
    return pd.DataFrame(rows)


def plot_quick_landscape(summary: pd.DataFrame, cfg: dict) -> None:
    if summary.empty:
        return
    classes = [
        "reversible_response",
        "developmental_delay",
        "fate_deflection",
        "off_manifold_collapse",
    ]
    colors = ["#2f7d5c", "#c98b2c", "#7d4fa3", "#b84343"]
    pivot = summary.pivot_table(index="stage_num", columns="response_recovery_class", values="class_fraction", fill_value=0.0)
    for klass in classes:
        if klass not in pivot:
            pivot[klass] = 0.0
    fig, ax = plt.subplots(figsize=(10, 5))
    bottom = np.zeros(len(pivot))
    for klass, color in zip(classes, colors):
        values = pivot[klass].to_numpy()
        ax.bar(pivot.index.astype(str), values, bottom=bottom, color=color, label=klass)
        bottom += values
    ax.set_xlabel("Theiler stage")
    ax.set_ylabel("Class fraction")
    ax.set_title("Response-recovery quick landscape")
    ax.legend(frameon=False, ncol=2, fontsize=8)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    figure_dir = resolve_project_path(cfg["output"]["figures"])
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / "response_recovery_landscape_quick.png", dpi=220)
    plt.close(fig)


def run_quick(cfg: dict) -> None:
    stage_df, priority = ensure_lite_outputs(cfg)
    centroids = proxy_centroids(stage_df)
    transferred = proxy_transferred(stage_df, priority, cfg)
    classes = classify_from_latent_tables(transferred, centroids, cfg)
    summary = summarize_response_recovery(classes)

    paths = cfg["output"]
    write_table(centroids, paths["stage_celltype_centroids"])
    write_table(transferred, "results/perturbation_transfer/tables/transferred_response_by_stage_celltype.csv")
    classes_path = write_table(classes, paths["response_recovery_classes"])
    summary_path = write_table(summary, paths["stage_vulnerability_response_recovery"])
    plot_quick_landscape(summary, cfg)
    print(f"Quick response-recovery classes: {classes_path}")
    print(f"Quick stage summary: {summary_path}")


def run_main(args: argparse.Namespace) -> None:
    run_script("build_developmental_manifold.py", "--config", args.response_config)
    export_args = ["--config", args.transfer_config]
    if args.allow_missing_external:
        export_args.append("--allow-missing")
    run_script("export_external_perturbation_response.py", *export_args)
    transfer_cfg = load_json(args.transfer_config)
    response_dict = resolve_project_path(transfer_cfg["output"]["external_response_dictionary"])
    response_dict_csv = response_dict.with_suffix(".csv")
    if args.allow_missing_external and (not response_dict.exists()) and (not response_dict_csv.exists()):
        print("External response dictionary is not available yet; stopped after manifold build and acquisition manifest.")
        return
    run_script("transfer_perturbation_response_ot.py", "--config", args.transfer_config)
    run_script("classify_response_recovery.py", "--response-config", args.response_config, "--transfer-config", args.transfer_config)
    run_script("analyze_response_recovery_statistics.py", "--config", args.response_config)
    run_script("compute_minimal_rescue_control.py", "--response-config", args.response_config, "--transfer-config", args.transfer_config)


def main() -> None:
    args = parse_args()
    cfg = load_json(args.response_config)
    if args.mode == "quick":
        run_quick(cfg)
        return
    run_main(args)
    if args.mode == "full":
        run_script("build_niche_context.py", "--config", args.niche_config)


if __name__ == "__main__":
    main()
