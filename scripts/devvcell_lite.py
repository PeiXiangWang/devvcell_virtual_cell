"""DevVCell-lite prototype analysis.

This script reuses lightweight RDEG-derived tables and converts them into
prototype metrics for a developmental virtual cell project.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = PROJECT_ROOT / "config" / "devvcell_lite.json"


def load_config() -> dict:
    with CONFIG_PATH.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def stage_number(label: str) -> int:
    text = str(label)
    digits = "".join(ch if ch.isdigit() else " " for ch in text).split()
    if not digits:
        raise ValueError(f"Cannot parse stage number from {label!r}")
    return int(digits[-1])


def canonical_stage(label: str) -> str:
    return f"Theiler stage {stage_number(label)}"


def minmax(values: pd.Series) -> pd.Series:
    series = pd.to_numeric(values, errors="coerce").astype(float)
    lo = series.min(skipna=True)
    hi = series.max(skipna=True)
    if not np.isfinite(lo) or not np.isfinite(hi) or math.isclose(lo, hi):
        return pd.Series(np.zeros(len(series)), index=series.index)
    return (series - lo) / (hi - lo)


def weighted_sum(frame: pd.DataFrame, weights: dict[str, float]) -> pd.Series:
    score = pd.Series(np.zeros(len(frame)), index=frame.index, dtype=float)
    for column, weight in weights.items():
        if column in frame.columns:
            score = score + frame[column].fillna(0.0).astype(float) * float(weight)
    return score


def stage_sort_key(labels: Iterable[str]) -> list[int]:
    return [stage_number(label) for label in labels]


def load_tables(data_dir: Path) -> dict[str, pd.DataFrame]:
    return {
        "stage_module": pd.read_csv(data_dir / "stage_module_means.csv"),
        "temporal": pd.read_csv(data_dir / "temporal_sensitivity.csv"),
        "tf_knockout": pd.read_csv(data_dir / "tf_knockout_results.csv"),
        "rescue": pd.read_csv(data_dir / "rescue_experiments.csv"),
        "rollout_error": pd.read_csv(data_dir / "rollout_error_matrix.csv"),
        "ot_pair": pd.read_csv(data_dir / "ot_pair_metrics.csv"),
        "next_step": pd.read_csv(data_dir / "next_step_pair_metrics.csv"),
        "nodes": pd.read_csv(data_dir / "nodes_summary.csv"),
    }


def compute_stage_module_shift(stage_module: pd.DataFrame) -> pd.DataFrame:
    module_cols = [c for c in stage_module.columns if c.startswith("module_")]
    ordered = stage_module.copy()
    ordered["stage"] = ordered["stage"].map(canonical_stage)
    ordered["stage_num"] = ordered["stage"].map(stage_number)
    ordered = ordered.sort_values("stage_num").reset_index(drop=True)

    distances = []
    for idx, row in ordered.iterrows():
        if idx == len(ordered) - 1:
            distances.append(np.nan)
            continue
        current = row[module_cols].to_numpy(dtype=float)
        nxt = ordered.loc[idx + 1, module_cols].to_numpy(dtype=float)
        distances.append(float(np.linalg.norm(nxt - current)))
    return pd.DataFrame(
        {
            "stage": ordered["stage"],
            "stage_num": ordered["stage_num"],
            "baseline_dev_step_norm": distances,
        }
    )


def compute_stage_vulnerability(tables: dict[str, pd.DataFrame], config: dict) -> pd.DataFrame:
    temporal = tables["temporal"].copy()
    temporal["stage"] = temporal["stage"].map(canonical_stage)
    temporal_summary = (
        temporal.groupby("stage", as_index=False)
        .agg(
            mean_sensitivity=("sensitivity", "mean"),
            max_sensitivity=("sensitivity", "max"),
            n_sensitive_tfs=("tf_name", "nunique"),
        )
    )

    rollout = tables["rollout_error"].copy()
    rollout = rollout.rename(columns={rollout.columns[0]: "stage"})
    rollout["stage"] = rollout["stage"].map(canonical_stage)
    step_cols = [c for c in rollout.columns if c.startswith("step_")]
    rollout["mean_rollout_error"] = rollout[step_cols].mean(axis=1)
    rollout["max_rollout_error"] = rollout[step_cols].max(axis=1)

    nodes = tables["nodes"].copy()
    nodes["stage"] = nodes["stage"].map(canonical_stage)
    node_summary = (
        nodes.groupby("stage", as_index=False)
        .agg(
            n_nodes=("node_id", "nunique"),
            total_cells=("cell_count", "sum"),
            mean_outlier_ratio=("outlier_ratio", "mean"),
        )
    )

    module_shift = compute_stage_module_shift(tables["stage_module"])

    next_step = tables["next_step"].copy()
    next_step["stage"] = next_step["src_stage"].map(canonical_stage)
    next_step["next_step_mse"] = next_step["mse_with_reg"].astype(float)
    next_step["ot_reg_gain"] = (
        next_step["mse_linear"].astype(float) - next_step["mse_with_reg"].astype(float)
    ) / next_step["mse_linear"].replace(0, np.nan).astype(float)
    next_step = next_step[["stage", "next_step_mse", "ot_reg_gain"]]

    stage_df = temporal_summary.merge(rollout[["stage", "mean_rollout_error", "max_rollout_error"]], on="stage", how="outer")
    stage_df = stage_df.merge(node_summary, on="stage", how="outer")
    stage_df = stage_df.merge(module_shift, on="stage", how="outer")
    stage_df = stage_df.merge(next_step, on="stage", how="left")
    stage_df["stage_num"] = stage_df["stage"].map(stage_number)

    fill_cols = [
        "mean_sensitivity",
        "max_sensitivity",
        "mean_rollout_error",
        "max_rollout_error",
        "mean_outlier_ratio",
        "baseline_dev_step_norm",
        "next_step_mse",
        "ot_reg_gain",
    ]
    for column in fill_cols:
        stage_df[column] = pd.to_numeric(stage_df[column], errors="coerce")
        stage_df[column] = stage_df[column].fillna(stage_df[column].median(skipna=True))

    for column in [
        "mean_sensitivity",
        "mean_rollout_error",
        "baseline_dev_step_norm",
        "mean_outlier_ratio",
        "next_step_mse",
    ]:
        stage_df[f"{column}_norm"] = minmax(stage_df[column])

    weights = config["vulnerability_weights"]
    stage_df["vulnerability_score"] = weighted_sum(stage_df, weights)

    window = config["stage_window_of_interest"]
    stage_df["window"] = np.where(
        stage_df["stage_num"].between(window["min_stage"], window["max_stage"]),
        window["name"],
        "outside",
    )

    output_cols = [
        "stage",
        "stage_num",
        "window",
        "vulnerability_score",
        "mean_sensitivity",
        "mean_rollout_error",
        "baseline_dev_step_norm",
        "mean_outlier_ratio",
        "next_step_mse",
        "ot_reg_gain",
        "n_sensitive_tfs",
        "n_nodes",
        "total_cells",
    ]
    return stage_df[output_cols].sort_values("stage_num").reset_index(drop=True)


def compute_perturbation_priority(tables: dict[str, pd.DataFrame], config: dict) -> pd.DataFrame:
    tf = tables["tf_knockout"].copy()
    temporal = tables["temporal"].copy()
    temporal["stage"] = temporal["stage"].map(canonical_stage)
    temporal["stage_num"] = temporal["stage"].map(stage_number)

    window = config["stage_window_of_interest"]
    temporal["in_window"] = temporal["stage_num"].between(window["min_stage"], window["max_stage"])

    idx = temporal.groupby("tf_name")["sensitivity"].idxmax()
    peak = temporal.loc[idx, ["tf_name", "stage", "sensitivity"]].rename(
        columns={"stage": "peak_stage", "sensitivity": "peak_sensitivity"}
    )

    temporal_summary = (
        temporal.groupby("tf_name", as_index=False)
        .agg(
            mean_sensitivity=("sensitivity", "mean"),
            max_sensitivity=("sensitivity", "max"),
        )
        .merge(peak, on="tf_name", how="left")
    )
    window_summary = (
        temporal.loc[temporal["in_window"]]
        .groupby("tf_name", as_index=False)
        .agg(window_sensitivity=("sensitivity", "mean"))
    )

    rescue = tables["rescue"].copy()
    rescue_best = (
        rescue.sort_values("rescue_fraction", ascending=False)
        .groupby("knocked_out_tf", as_index=False)
        .first()
        .rename(
            columns={
                "knocked_out_tf": "tf_name",
                "rescuer_tf": "best_rescuer_tf",
                "rescue_fraction": "best_rescue_fraction",
            }
        )
    )

    priority = tf.merge(temporal_summary, on="tf_name", how="left")
    priority = priority.merge(window_summary, on="tf_name", how="left")
    priority = priority.merge(rescue_best, on="tf_name", how="left")
    priority["best_rescue_fraction"] = priority["best_rescue_fraction"].fillna(0.0)
    priority["best_rescuer_tf"] = priority["best_rescuer_tf"].fillna("not_observed")
    priority["window_sensitivity"] = priority["window_sensitivity"].fillna(priority["mean_sensitivity"])
    priority["feedback_cost"] = 1.0 - priority["best_rescue_fraction"].clip(lower=0.0, upper=1.0)
    priority["log_mass_shift"] = np.log10(priority["mass_shift_mean"].abs() + 1.0)

    norm_map = {
        "global_effect_score": "global_effect_norm",
        "developmental_delay_score": "developmental_delay_norm",
        "window_sensitivity": "window_sensitivity_norm",
        "feedback_cost": "feedback_cost_norm",
        "log_mass_shift": "log_mass_shift_norm",
    }
    for src, dst in norm_map.items():
        priority[dst] = minmax(priority[src])

    priority["devvcell_priority_score"] = weighted_sum(priority, config["perturbation_weights"])
    priority["response_amplitude_proxy"] = priority["global_effect_norm"]
    priority["fate_displacement_proxy"] = 0.5 * priority["developmental_delay_norm"] + 0.5 * priority["log_mass_shift_norm"]
    priority["recovery_probability_proxy"] = np.exp(-2.5 * priority["feedback_cost_norm"])

    output_cols = [
        "tf_name",
        "devvcell_priority_score",
        "response_amplitude_proxy",
        "fate_displacement_proxy",
        "feedback_cost",
        "recovery_probability_proxy",
        "global_effect_score",
        "developmental_delay_score",
        "mass_shift_mean",
        "window_sensitivity",
        "peak_stage",
        "peak_sensitivity",
        "best_rescuer_tf",
        "best_rescue_fraction",
        "n_affected_genes",
    ]
    return priority[output_cols].sort_values("devvcell_priority_score", ascending=False).reset_index(drop=True)


def compute_virtual_cell_rollout_proxy(
    stage_vulnerability: pd.DataFrame,
    perturbation_priority: pd.DataFrame,
    tables: dict[str, pd.DataFrame],
    config: dict,
) -> pd.DataFrame:
    top_tfs = perturbation_priority.head(5).copy()
    stage_module = compute_stage_module_shift(tables["stage_module"])
    base = stage_vulnerability.merge(stage_module[["stage", "baseline_dev_step_norm"]], on="stage", how="left", suffixes=("", "_raw"))
    rows = []

    for _, stage_row in base.iterrows():
        for _, tf_row in top_tfs.iterrows():
            window_boost = float(stage_row["vulnerability_score"])
            response = float(stage_row["baseline_dev_step_norm"]) * (1.0 + float(tf_row["response_amplitude_proxy"]) + window_boost)
            recovery_cost = response * float(tf_row["feedback_cost"])
            recovery_prob = math.exp(-min(5.0, recovery_cost * 12.0))
            rows.append(
                {
                    "stage": stage_row["stage"],
                    "stage_num": int(stage_row["stage_num"]),
                    "tf_name": tf_row["tf_name"],
                    "virtual_cell_response_proxy": response,
                    "feedback_cost_proxy": recovery_cost,
                    "recovery_probability_proxy": recovery_prob,
                    "baseline_dev_step_norm": stage_row["baseline_dev_step_norm"],
                    "stage_vulnerability_score": stage_row["vulnerability_score"],
                    "tf_priority_score": tf_row["devvcell_priority_score"],
                }
            )

    return pd.DataFrame(rows).sort_values(
        ["virtual_cell_response_proxy", "feedback_cost_proxy"], ascending=False
    )


def plot_stage_vulnerability(stage_df: pd.DataFrame, figures_dir: Path) -> None:
    colors = ["#b84343" if w != "outside" else "#4b7ba7" for w in stage_df["window"]]
    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.bar(stage_df["stage_num"].astype(str), stage_df["vulnerability_score"], color=colors)
    ax.set_xlabel("Theiler stage")
    ax.set_ylabel("DevVCell vulnerability score")
    ax.set_title("Prototype developmental vulnerability by stage")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(figures_dir / "stage_vulnerability.png", dpi=200)
    plt.close(fig)


def plot_perturbation_priority(priority: pd.DataFrame, figures_dir: Path, top_n: int) -> None:
    top = priority.head(top_n).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.barh(top["tf_name"], top["devvcell_priority_score"], color="#3f7f5f")
    ax.set_xlabel("DevVCell perturbation priority")
    ax.set_title(f"Top {top_n} TF perturbation candidates")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(figures_dir / f"perturbation_priority_top{top_n}.png", dpi=200)
    plt.close(fig)


def plot_recovery_scatter(priority: pd.DataFrame, figures_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.scatter(
        priority["response_amplitude_proxy"],
        priority["feedback_cost"],
        s=30 + priority["fate_displacement_proxy"] * 120,
        c=priority["devvcell_priority_score"],
        cmap="viridis",
        alpha=0.85,
        edgecolor="white",
        linewidth=0.5,
    )
    for _, row in priority.head(8).iterrows():
        ax.text(row["response_amplitude_proxy"], row["feedback_cost"], row["tf_name"], fontsize=8)
    ax.set_xlabel("Response amplitude proxy")
    ax.set_ylabel("Feedback cost proxy")
    ax.set_title("Perturbation response and recovery cost")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(figures_dir / "recovery_cost_vs_response.png", dpi=200)
    plt.close(fig)


def plot_module_step_distance(stage_df: pd.DataFrame, figures_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.plot(stage_df["stage_num"], stage_df["baseline_dev_step_norm"], marker="o", color="#6f4e7c")
    ax.set_xlabel("Theiler stage")
    ax.set_ylabel("Module step distance")
    ax.set_title("Normal-development module shift proxy")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(figures_dir / "module_step_distance.png", dpi=200)
    plt.close(fig)


def write_summary(
    config: dict,
    stage_vulnerability: pd.DataFrame,
    perturbation_priority: pd.DataFrame,
    rollout_proxy: pd.DataFrame,
    tables_dir: Path,
) -> None:
    window_name = config["stage_window_of_interest"]["name"]
    window_mean = float(stage_vulnerability.loc[stage_vulnerability["window"] == window_name, "vulnerability_score"].mean())
    outside_mean = float(stage_vulnerability.loc[stage_vulnerability["window"] == "outside", "vulnerability_score"].mean())

    summary = {
        "project": config["project"],
        "source_run": config["source_run"],
        "n_stages": int(stage_vulnerability["stage"].nunique()),
        "n_tf_candidates": int(perturbation_priority["tf_name"].nunique()),
        "top_vulnerability_stage": stage_vulnerability.sort_values("vulnerability_score", ascending=False).iloc[0].to_dict(),
        "window_vulnerability_mean": window_mean,
        "outside_vulnerability_mean": outside_mean,
        "window_vs_outside_ratio": window_mean / outside_mean if outside_mean else None,
        "top_perturbation": perturbation_priority.iloc[0].to_dict(),
        "top_rollout_proxy": rollout_proxy.iloc[0].to_dict(),
        "outputs": {
            "stage_vulnerability": "results/tables/stage_vulnerability.csv",
            "perturbation_priority": "results/tables/perturbation_priority.csv",
            "virtual_cell_rollout_proxy": "results/tables/virtual_cell_rollout_proxy.csv",
        },
    }
    with (tables_dir / "execution_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)


def main() -> None:
    config = load_config()
    data_dir = PROJECT_ROOT / config["data_dir"]
    results_dir = PROJECT_ROOT / config["results_dir"]
    tables_dir = results_dir / "tables"
    figures_dir = results_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    tables = load_tables(data_dir)
    stage_vulnerability = compute_stage_vulnerability(tables, config)
    perturbation_priority = compute_perturbation_priority(tables, config)
    rollout_proxy = compute_virtual_cell_rollout_proxy(stage_vulnerability, perturbation_priority, tables, config)

    stage_vulnerability.to_csv(tables_dir / "stage_vulnerability.csv", index=False)
    perturbation_priority.to_csv(tables_dir / "perturbation_priority.csv", index=False)
    rollout_proxy.to_csv(tables_dir / "virtual_cell_rollout_proxy.csv", index=False)

    top_n = int(config.get("top_n_perturbations", 12))
    plot_stage_vulnerability(stage_vulnerability, figures_dir)
    plot_perturbation_priority(perturbation_priority, figures_dir, top_n=top_n)
    plot_recovery_scatter(perturbation_priority, figures_dir)
    plot_module_step_distance(stage_vulnerability, figures_dir)
    write_summary(config, stage_vulnerability, perturbation_priority, rollout_proxy, tables_dir)

    print("DevVCell-lite prototype completed.")
    print(f"Top stage: {stage_vulnerability.sort_values('vulnerability_score', ascending=False).iloc[0]['stage']}")
    print(f"Top perturbation: {perturbation_priority.iloc[0]['tf_name']}")
    print(f"Results: {results_dir}")


if __name__ == "__main__":
    main()
