"""Run multi-seed and ablation checks for the DevVCell research package."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from copy import deepcopy
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from devvcell.io import load_json, resolve_project_path, write_json  # noqa: E402


plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial Unicode MS", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/ablation_suite.json")
    parser.add_argument(
        "--skip-transition",
        action="store_true",
        help="Skip transition ablations and only aggregate existing outputs if present.",
    )
    parser.add_argument(
        "--skip-stimulus",
        action="store_true",
        help="Skip stimulus ablations and only aggregate existing outputs if present.",
    )
    return parser.parse_args()


def run_command(command: list[str]) -> None:
    print("Running:", " ".join(command))
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)


def write_temp_config(path: Path, config: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(config, handle, ensure_ascii=False, indent=2)


def transition_configs(base: dict, suite: dict, results_dir: Path) -> list[tuple[str, Path, dict]]:
    configs: list[tuple[str, Path, dict]] = []
    tcfg = suite["transition"]
    for seed in tcfg["seeds"]:
        for method in tcfg["pairing_methods"]:
            run_name = f"transition_seed{seed}_{method}"
            cfg = deepcopy(base)
            cfg["seed"] = int(seed)
            cfg["results_dir"] = str(results_dir / "runs" / run_name)
            cfg["model"]["latent_dim"] = int(tcfg["latent_dim"])
            cfg["model"]["max_pairs_per_transition"] = int(tcfg["max_pairs_per_transition"])
            cfg["model"]["pairing"]["method"] = method
            cfg["model"]["pairing"]["sinkhorn_epsilon"] = float(tcfg["sinkhorn_epsilon"])
            cfg["model"]["pairing"]["sinkhorn_iterations"] = int(tcfg["sinkhorn_iterations"])
            cfg["model"]["mlp"]["enabled"] = bool(tcfg["mlp_enabled"])
            if "context_residual_mlp" in tcfg:
                cfg["model"]["context_residual_mlp"] = deepcopy(tcfg["context_residual_mlp"])
            cfg_path = results_dir / "configs" / f"{run_name}.json"
            configs.append((run_name, cfg_path, cfg))
    return configs


def run_transition_ablations(base: dict, suite: dict, results_dir: Path, skip: bool) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for run_name, cfg_path, cfg in transition_configs(base, suite, results_dir):
        write_temp_config(cfg_path, cfg)
        run_dir = resolve_project_path(cfg["results_dir"])
        summary_path = run_dir / "training_summary.json"
        metrics_path = run_dir / "tables" / "cell_level_transition_metrics.csv"
        if not skip:
            run_command([sys.executable, "scripts/train_cell_transition_baseline.py", "--config", str(cfg_path)])
        if not summary_path.exists() or not metrics_path.exists():
            raise FileNotFoundError(f"Missing transition ablation output for {run_name}")
        metrics = pd.read_csv(metrics_path)
        for row in metrics.itertuples(index=False):
            rows.append(
                {
                    "run_name": run_name,
                    "seed": int(cfg["seed"]),
                    "pairing_method": cfg["model"]["pairing"]["method"],
                    "latent_dim": int(cfg["model"]["latent_dim"]),
                    "max_pairs_per_transition": int(cfg["model"]["max_pairs_per_transition"]),
                    "model": row.model,
                    "system": row.system,
                    "src_stage": int(row.src_stage),
                    "tgt_stage": int(row.tgt_stage),
                    "n_pairs": int(row.n_pairs),
                    "pair_latent_mse": float(row.pair_latent_mse),
                    "centroid_latent_mse": float(row.centroid_latent_mse),
                    "rbf_mmd": float(row.rbf_mmd),
                }
            )
    table = pd.DataFrame(rows)
    table.to_csv(results_dir / "tables" / "transition_ablation_metrics.csv", index=False)
    summary = (
        table.groupby(["pairing_method", "model"], as_index=False)
        .agg(
            mean_pair_latent_mse=("pair_latent_mse", "mean"),
            sd_pair_latent_mse=("pair_latent_mse", "std"),
            mean_centroid_latent_mse=("centroid_latent_mse", "mean"),
            mean_rbf_mmd=("rbf_mmd", "mean"),
            n_evaluations=("pair_latent_mse", "size"),
        )
        .sort_values(["pairing_method", "mean_pair_latent_mse"])
    )
    summary.to_csv(results_dir / "tables" / "transition_ablation_summary.csv", index=False)
    return summary


def stimulus_configs(base: dict, suite: dict, results_dir: Path) -> list[tuple[str, Path, dict, Path]]:
    configs: list[tuple[str, Path, dict, Path]] = []
    scfg = suite["stimulus"]
    for variant in scfg["variants"]:
        cfg = deepcopy(base)
        cfg["stimulus"]["top_n_tfs"] = int(scfg["top_n_tfs"])
        cfg["stimulus"]["max_cells_per_stage_system"] = int(scfg["max_cells_per_stage_system"])
        cfg["stimulus"]["global_grn_weight"] = float(variant["global_grn_weight"])
        cfg["stimulus"]["system_edge_weight"] = float(variant["system_edge_weight"])
        run_name = str(variant["name"])
        out_dir = results_dir / "runs" / f"stimulus_{run_name}"
        cfg_path = results_dir / "configs" / f"stimulus_{run_name}.json"
        configs.append((run_name, cfg_path, cfg, out_dir))
    return configs


def run_stimulus_ablations(base: dict, suite: dict, results_dir: Path, skip: bool) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for run_name, cfg_path, cfg, out_dir in stimulus_configs(base, suite, results_dir):
        write_temp_config(cfg_path, cfg)
        response_path = out_dir / "tables" / "cell_level_tf_grn_stimulus_response.csv"
        tf_summary_path = out_dir / "tables" / "cell_level_tf_grn_stimulus_summary.csv"
        if not skip:
            run_command(
                [
                    sys.executable,
                    "scripts/run_stimulus_response_head.py",
                    "--config",
                    str(cfg_path),
                    "--model-results-dir",
                    base["results_dir"],
                    "--output-dir",
                    str(out_dir),
                ]
            )
        if not response_path.exists() or not tf_summary_path.exists():
            raise FileNotFoundError(f"Missing stimulus ablation output for {run_name}")
        tf_summary = pd.read_csv(tf_summary_path)
        for row in tf_summary.itertuples(index=False):
            rows.append(
                {
                    "variant": run_name,
                    "global_grn_weight": float(cfg["stimulus"]["global_grn_weight"]),
                    "system_edge_weight": float(cfg["stimulus"]["system_edge_weight"]),
                    "tf_name": row.tf_name,
                    "mean_stimulus_response_norm": float(row.mean_stimulus_response_norm),
                    "max_stimulus_response_norm": float(row.max_stimulus_response_norm),
                    "mean_fate_displacement": float(row.mean_fate_displacement),
                    "mean_recovery_probability": float(row.mean_recovery_probability),
                    "n_system_stage_tests": int(row.n_system_stage_tests),
                }
            )
    table = pd.DataFrame(rows)
    table.to_csv(results_dir / "tables" / "stimulus_ablation_tf_summary.csv", index=False)
    summary = (
        table.groupby("variant", as_index=False)
        .agg(
            mean_response=("mean_stimulus_response_norm", "mean"),
            max_response=("max_stimulus_response_norm", "max"),
            mean_recovery_probability=("mean_recovery_probability", "mean"),
            n_tfs=("tf_name", "nunique"),
        )
        .sort_values("mean_response", ascending=False)
    )
    summary.to_csv(results_dir / "tables" / "stimulus_ablation_summary.csv", index=False)
    return table


def plot_transition_summary(summary: pd.DataFrame, figures_dir: Path) -> None:
    model_labels = {
        "ridge": "岭回归",
        "identity": "恒等映射",
        "mean_shift": "均值平移",
    }
    pairing_labels = {
        "nearest": "最近邻",
        "sinkhorn": "Sinkhorn OT",
    }
    labeled = summary.copy()
    labeled["model_label"] = labeled["model"].map(model_labels).fillna(labeled["model"])
    labeled["pairing_label"] = labeled["pairing_method"].map(pairing_labels).fillna(labeled["pairing_method"])
    pivot = labeled.pivot(index="model_label", columns="pairing_label", values="mean_pair_latent_mse").fillna(0.0)
    pivot = pivot.loc[pivot.mean(axis=1).sort_values().index]
    fig, ax = plt.subplots(figsize=(7.4, 4.5))
    pivot.plot(kind="bar", ax=ax, color=["#4e79a7", "#59a14f", "#b07aa1", "#f28e2b"][: len(pivot.columns)])
    ax.set_ylabel("平均 heldout latent MSE")
    ax.set_xlabel("transition 模型")
    ax.set_title("多 seed transition 消融")
    ax.legend(title="配对方法")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(figures_dir / "transition_ablation_summary.png", dpi=220)
    plt.close(fig)


def plot_stimulus_summary(table: pd.DataFrame, figures_dir: Path) -> None:
    variant_labels = {
        "full_grn": "完整 GRN",
        "no_global_grn": "仅系统边",
        "no_system_edges": "仅全局 GRN",
    }
    labeled = table.copy()
    labeled["variant_label"] = labeled["variant"].map(variant_labels).fillna(labeled["variant"])
    pivot = labeled.pivot(index="tf_name", columns="variant_label", values="mean_stimulus_response_norm").fillna(0.0)
    order = pivot.mean(axis=1).sort_values(ascending=False).index
    pivot = pivot.loc[order]
    fig, ax = plt.subplots(figsize=(8.2, max(4.2, 0.32 * len(pivot))))
    im = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="cividis")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=25, ha="right")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_title("TF/GRN stimulus head 消融")
    ax.set_xlabel("消融设置")
    ax.set_ylabel("TF")
    fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, label="平均响应强度")
    fig.tight_layout()
    fig.savefig(figures_dir / "stimulus_ablation_heatmap.png", dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    suite = load_json(args.config)
    base = load_json(suite["base_config"])
    results_dir = resolve_project_path(suite["results_dir"])
    tables_dir = results_dir / "tables"
    figures_dir = results_dir / "figures"
    for path in [tables_dir, figures_dir, results_dir / "configs", results_dir / "runs"]:
        path.mkdir(parents=True, exist_ok=True)

    transition_summary = run_transition_ablations(base, suite, results_dir, skip=args.skip_transition)
    stimulus_table = run_stimulus_ablations(base, suite, results_dir, skip=args.skip_stimulus)
    plot_transition_summary(transition_summary, figures_dir)
    plot_stimulus_summary(stimulus_table, figures_dir)

    summary = {
        "analysis": "devvcell_ablation_suite_v1",
        "transition_runs": int(len(suite["transition"]["seeds"]) * len(suite["transition"]["pairing_methods"])),
        "stimulus_variants": [v["name"] for v in suite["stimulus"]["variants"]],
        "outputs": {
            "transition_metrics": "results/ablation_v1/tables/transition_ablation_metrics.csv",
            "transition_summary": "results/ablation_v1/tables/transition_ablation_summary.csv",
            "stimulus_tf_summary": "results/ablation_v1/tables/stimulus_ablation_tf_summary.csv",
            "stimulus_summary": "results/ablation_v1/tables/stimulus_ablation_summary.csv",
            "transition_figure": "results/ablation_v1/figures/transition_ablation_summary.png",
            "stimulus_figure": "results/ablation_v1/figures/stimulus_ablation_heatmap.png",
        },
    }
    write_json(results_dir / "ablation_summary.json", summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
