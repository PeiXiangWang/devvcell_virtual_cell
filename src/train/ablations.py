from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from src.utils.config import ensure_dir, load_config, write_text
from src.utils.stats import benjamini_hochberg, paired_test


PRIMARY_METRICS = ["sinkhorn", "mmd_rbf", "energy", "celltype_composition_rmse"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/train.yaml")
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    cfg = load_config(args.config)
    model_cfg = load_config(cfg.get("model_config", "configs/model.yaml"))
    if args.quick_fixture:
        model_cfg = dict(model_cfg)
        model_cfg["metrics_path"] = "tables/quick_fixture/final_metrics.csv"
        cfg = dict(cfg)
        cfg["ablation_metrics_path"] = "tables/quick_fixture/ablation_metrics.csv"
        cfg["ablation_stats_path"] = "tables/quick_fixture/ablation_statistical_tests.csv"
        cfg["negative_results_report"] = "reports/quick_fixture/negative_results.md"
        cfg["ablation_report"] = "reports/quick_fixture/ablation_interpretation.md"
    ensure_dir("figures")
    ensure_dir("tables")
    ensure_dir("reports")
    metrics = pd.read_csv(model_cfg.get("metrics_path", "tables/final_metrics.csv"))
    out_path = cfg.get("ablation_metrics_path", "tables/ablation_metrics.csv")
    metrics.to_csv(out_path, index=False)
    tests = []
    baseline_candidates = ["M0_linear_label_interpolation", "M1_intrinsic_neural", "M2_ot_teacher_force"]
    best_baseline = metrics[metrics["model"].isin(baseline_candidates)].groupby("model")["sinkhorn"].mean().sort_values().index[0]
    challenger = "M9_full_memory"
    for metric in PRIMARY_METRICS:
        tests.append(paired_test(metrics, metric, best_baseline, challenger))
    stat = pd.DataFrame(tests)
    stat["q_value"] = benjamini_hochberg(stat["p_value"].tolist())
    stat.to_csv(cfg.get("ablation_stats_path", "tables/ablation_statistical_tests.csv"), index=False)

    summary = metrics.groupby("model")[PRIMARY_METRICS].agg(["mean", "std"]).reset_index()
    ranked = metrics.groupby("model")["sinkhorn"].mean().sort_values()
    fig, ax = plt.subplots(figsize=(10, 5), dpi=160)
    order = ranked.index.tolist()
    means = metrics.groupby("model")["sinkhorn"].mean().reindex(order)
    std = metrics.groupby("model")["sinkhorn"].std().reindex(order)
    ax.bar(range(len(order)), means, yerr=std, color="#5b8db8", edgecolor="#234", linewidth=0.4)
    ax.set_xticks(range(len(order)), order, rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("Sinkhorn distance lower is better")
    ax.set_title("Held-out reconstruction ablation")
    fig.tight_layout()
    fig.savefig("figures/ablation_barplots.png")
    plt.close(fig)

    heat = metrics.groupby("model")[PRIMARY_METRICS].mean().reindex(order)
    heat_norm = (heat - heat.min()) / (heat.max() - heat.min()).replace(0, 1)
    fig, ax = plt.subplots(figsize=(7, 6), dpi=160)
    im = ax.imshow(heat_norm.to_numpy(), aspect="auto", cmap="magma_r")
    ax.set_yticks(range(len(order)), order, fontsize=7)
    ax.set_xticks(range(len(PRIMARY_METRICS)), PRIMARY_METRICS, rotation=45, ha="right")
    fig.colorbar(im, ax=ax, fraction=0.035, label="normalized error")
    ax.set_title("Ablation metric heatmap")
    fig.tight_layout()
    fig.savefig("figures/ablation_radar_or_heatmap.png")
    plt.close(fig)

    negative = metrics[metrics["model"].isin(["M10_shuffled_time_ot", "M11_random_lr_labels"])]
    neg_lines = [
        "# Negative Results",
        "",
        "Negative controls are retained and reported. Lower metric values are better for distribution distances.",
        "",
        negative.groupby("model")[PRIMARY_METRICS].mean().to_markdown(),
        "",
    ]
    write_text(cfg.get("negative_results_report", "reports/negative_results.md"), "\n".join(neg_lines))
    report = [
        "# Ablation Interpretation",
        "",
        "This report is generated directly from `tables/ablation_metrics.csv`; it is not a claim of publishable superiority.",
        "",
        "## Mean Metrics",
        "",
        summary.to_markdown(index=False),
        "",
        "## Primary Paired Tests",
        "",
        stat.to_markdown(index=False),
        "",
        "## Current Interpretation",
        "",
        f"Strongest non-reference baseline for paired diagnostic: `{best_baseline}`. `M0b_ot_interpolation` is retained as the OT teacher/reference interpolation, not as a competitor that the finite-agent model must beat. Modules are evaluated by teacher fidelity plus whether they provide stable mechanistic diagnostics and emergent-law signals.",
        "",
    ]
    write_text(cfg.get("ablation_report", "reports/ablation_interpretation.md"), "\n".join(report))
    print({"ablation_metrics": out_path, "stats": cfg.get("ablation_stats_path")})


if __name__ == "__main__":
    main()
