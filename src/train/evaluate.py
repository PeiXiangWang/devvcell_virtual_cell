from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from src.utils.config import ensure_dir, load_config, write_text


CORE = ["sinkhorn", "mmd_rbf", "energy", "celltype_composition_rmse"]


def _best_model(metrics: pd.DataFrame) -> str:
    score = metrics.groupby("model")[CORE].mean()
    ranks = score.rank(axis=0, ascending=True).mean(axis=1)
    return str(ranks.sort_values().index[0])


def _strongest_baseline(metrics: pd.DataFrame) -> str:
    candidates = ["M0_linear_label_interpolation", "M0b_ot_interpolation", "M1_intrinsic_neural", "M2_ot_teacher_force"]
    sub = metrics[metrics["model"].isin(candidates)]
    return _best_model(sub)


def _diagnose(metrics: pd.DataFrame, full: str, baseline: str) -> list[str]:
    mean = metrics.groupby("model")[CORE].mean()
    lines = []
    if full not in mean.index or baseline not in mean.index:
        return ["- Full model or strongest baseline missing from metrics."]
    diff = mean.loc[full] - mean.loc[baseline]
    worse = diff[diff > 0].sort_values(ascending=False)
    if worse.empty:
        lines.append("- Full model improves all core distance/composition metrics relative to the strongest baseline.")
    else:
        lines.append(f"- Full model is worse than `{baseline}` on: " + ", ".join(f"{k} (+{v:.4g})" for k, v in worse.items()))
        if "M7_ot_swarm_birth_death_diffusion" in mean.index and mean.loc["M7_ot_swarm_birth_death_diffusion", "sinkhorn"] > mean.loc[baseline, "sinkhorn"]:
            lines.append("- Diagnostic: swarm/birth/diffusion coefficients or noise are likely too strong for the current teacher; inspect event counts and sigma calibration.")
        if "M8_ot_swarm_birth_death_diffusion_cci" in mean.index and mean.loc["M8_ot_swarm_birth_death_diffusion_cci", "sinkhorn"] > mean.loc["M7_ot_swarm_birth_death_diffusion", "sinkhorn"]:
            lines.append("- Diagnostic: CCI graph is not improving reconstruction; LR edges may be sparse, mis-specified or not relevant to this dataset.")
        if "M9_full_memory" in mean.index and mean.loc["M9_full_memory", "sinkhorn"] > mean.loc["M8_ot_swarm_birth_death_diffusion_cci", "sinkhorn"]:
            lines.append("- Diagnostic: memory field is over-steering or poorly calibrated; no_memory vs memory ablation does not support retaining memory.")
    return lines


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/train.yaml")
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    cfg = load_config(args.config)
    model_cfg = load_config(cfg.get("model_config", "configs/model.yaml"))
    if args.quick_fixture:
        cfg = dict(cfg)
        model_cfg = dict(model_cfg)
        model_cfg["metrics_path"] = "tables/quick_fixture/final_metrics.csv"
        cfg["baseline_execution_matrix_path"] = "reports/quick_fixture/baseline_execution_matrix.csv"
        cfg["module_contribution_report"] = "reports/quick_fixture/module_contribution_audit.md"
        cfg["scientific_gap_report"] = "reports/quick_fixture/scientific_gap_audit.md"
    ensure_dir("reports/quick_fixture" if args.quick_fixture else "reports")
    metrics_path = model_cfg.get("metrics_path", "tables/final_metrics.csv")
    metrics = pd.read_csv(metrics_path)
    best = _best_model(metrics)
    strongest = _strongest_baseline(metrics)
    full = "M9_full_memory"
    mean = metrics.groupby("model")[["sinkhorn", "mmd_rbf", "energy", "knn_two_sample_accuracy", "celltype_composition_rmse"]].agg(["mean", "std"])
    ensure_dir("tables/quick_fixture" if args.quick_fixture else "tables")
    mean.to_csv("tables/quick_fixture/final_model_metric_summary.csv" if args.quick_fixture else "tables/final_model_metric_summary.csv")
    comparison = metrics[metrics["model"].isin([strongest, full])].groupby("model")[CORE].mean()
    supported = False
    wins = 0
    if strongest in comparison.index and full in comparison.index:
        wins = int((comparison.loc[full] < comparison.loc[strongest]).sum())
        supported = bool(wins >= 2)
    diagnostics = _diagnose(metrics, full, strongest)
    baseline_path = cfg.get("baseline_execution_matrix_path", "reports/baseline_execution_matrix.csv")
    baseline_matrix = pd.read_csv(baseline_path) if Path(baseline_path).exists() else pd.DataFrame()

    module_report = [
        "# Module Contribution Audit",
        "",
        f"- strongest baseline: `{strongest}`",
        f"- full model: `{full}`",
        f"- full-model core metric wins over strongest baseline: {wins}/4",
        "",
        "## Mean Metrics",
        "",
        metrics.groupby("model")[CORE].mean().to_markdown(),
        "",
        "## Diagnostics",
        "",
        *diagnostics,
        "",
    ]
    write_text(cfg.get("module_contribution_report", "reports/module_contribution_audit.md"), "\n".join(module_report))

    scientific_gap = [
        "# Scientific Gap Audit",
        "",
        f"- Best mean-rank model: `{best}`.",
        f"- Strongest baseline: `{strongest}`.",
        f"- Full model passes the predefined gate of beating the strongest baseline on at least two core metrics: {supported}.",
        "- High-level claims are allowed only with native_moscot/native_wot or externally validated teacher. The current fallback teacher remains toy_sinkhorn_fallback unless reports say otherwise.",
        "",
        "## Baseline Execution Matrix",
        "",
        baseline_matrix.to_markdown(index=False) if not baseline_matrix.empty else "No baseline execution matrix found.",
        "",
        "## Required Next Steps Before Strong Claims",
        "",
        "- Run native moscot TemporalProblem and extract native transport plans, or validate the fallback teacher externally.",
        "- Add real external held-out dataset validation and lineage/perturbation validation.",
        "- Tune trainable swarm, event and memory coefficients only on training times, then re-run the strict holdout gate.",
        "- Keep negative controls in the report; shuffled time/LR controls must degrade performance.",
        "",
    ]
    write_text(cfg.get("scientific_gap_report", "reports/scientific_gap_audit.md"), "\n".join(scientific_gap))

    editorial = [
        "# Editorial Assessment",
        "",
        f"- Best mean-rank model in current run: `{best}`.",
        f"- Full model beats strongest baseline `{strongest}` on at least two core metrics: {supported}.",
        "- Current evidence level: computational research prototype.",
        "",
        "## Readiness",
        "",
        "Not sufficient for Nature/Nature Methods/Nature Biotechnology unless the gate is passed with native or externally validated teacher and strong baselines.",
        "",
    ]
    write_text("reports/quick_fixture/editorial_assessment.md" if args.quick_fixture else "reports/editorial_assessment.md", "\n".join(editorial))

    retained = [
        "# Final Retained Results and Methods",
        "",
        "This document reports only retained methods and results from the current run. It does not present fallback teacher results as native moscot/WOT.",
        "",
        "## Data Split",
        "",
        "The pipeline supports `strict_time_holdout` and `teacher_edge_holdout`. Under strict holdout, held-out cells are excluded from HVG/SVD fitting, teacher construction, model training and hyperparameter selection; they are transformed only for evaluation.",
        "",
        "## Teacher",
        "",
        "Teacher backend is recorded in `processed/**/ot_teacher_summary.json` and AnnData `uns['swarmlineage_ot_teacher']['backend']`. toy_sinkhorn_fallback is a toy fallback and cannot support high-level moscot/WOT claims.",
        "",
        "## Model",
        "",
        "M1 intrinsic dynamics is trained without OT velocity targets. M2 and later variants can use OT teacher velocity. Birth/death uses stochastic event simulation; memory is a fate-specific kNN field; CCI uses sender-receiver LR graph signals.",
        "",
        "## Results",
        "",
        f"- best mean-rank model: `{best}`",
        f"- strongest baseline: `{strongest}`",
        f"- full model gate passed: {supported}",
        "",
        metrics.groupby("model")[CORE].mean().to_markdown(),
        "",
        "## Interpretation",
        "",
        "\n".join(diagnostics),
        "",
        "## Limitations",
        "",
        "- No result is written as a success if the full model does not beat the strongest baseline.",
        "- Native CellRank2/TrajectoryNet/MIOFlow/TIGON are only counted as compared methods if marked `executed=True` in the baseline matrix.",
        "- Perturbation and CCI results are exploratory until validated by matched perturbation or spatial data.",
        "",
    ]
    write_text("manuscript/final_retained_results_and_methods.md", "\n".join(retained))
    print({"best_model": best, "strongest_baseline": strongest, "full_model_gate_supported": supported, "wins": wins})


if __name__ == "__main__":
    main()
