from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from src.utils.config import load_config, write_text


def _best_model(metrics: pd.DataFrame) -> str:
    score = metrics.groupby("model")[["sinkhorn", "mmd_rbf", "energy", "celltype_composition_rmse"]].mean()
    ranks = score.rank(axis=0, ascending=True).mean(axis=1)
    return str(ranks.sort_values().index[0])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/train.yaml")
    args = parser.parse_args()
    cfg = load_config(args.config)
    model_cfg = load_config(cfg.get("model_config", "configs/model.yaml"))
    metrics_path = model_cfg.get("metrics_path", "tables/final_metrics.csv")
    metrics = pd.read_csv(metrics_path)
    best = _best_model(metrics)
    mean = metrics.groupby("model")[["sinkhorn", "mmd_rbf", "energy", "knn_two_sample_accuracy", "celltype_composition_rmse"]].agg(["mean", "std"])
    mean.to_csv("tables/final_model_metric_summary.csv")
    baseline = "M0_ot_interpolation"
    full = "M9_full_pheromone"
    comparison = metrics[metrics["model"].isin([baseline, full])].groupby("model")[["sinkhorn", "mmd_rbf", "energy", "celltype_composition_rmse"]].mean()
    supported = False
    if baseline in comparison.index and full in comparison.index:
        supported = bool((comparison.loc[full] < comparison.loc[baseline]).sum() >= 2)
    editorial = [
        "# Editorial Assessment",
        "",
        f"- Best mean-rank model in current quick run: `{best}`.",
        f"- Full model beats OT interpolation on at least two core metrics: {supported}.",
        "- Current evidence level: computational prototype on a subsampled real mouse developmental dataset.",
        "",
        "## Nature/Nature Methods/Nature Biotechnology Readiness",
        "",
    ]
    if supported:
        editorial.append("The quick run supports further development, but it is still not sufficient for a Nature-level submission without native moscot/WOT runs, external validation, lineage or perturbation validation, and stronger baseline execution.")
    else:
        editorial.append("Not sufficient for Nature-level submission. The current evidence does not pass the predefined superiority gate for the full model, or the comparison remains too lightweight.")
    editorial += [
        "",
        "## Largest Shortfalls",
        "",
        "- Native moscot, WOT, TIGON, TrajectoryNet/MIOFlow and CellRank2 baselines are not all executed end-to-end in this quick prototype.",
        "- The teacher is pseudo-lineage from stage snapshots, not lineage tracing.",
        "- Perturbation validation is exploratory unless a matched perturbation time-series is added.",
        "- The agent simulator is intentionally minimal and requires scalability and hyperparameter sensitivity before a methods-journal claim.",
        "",
    ]
    write_text("reports/editorial_assessment.md", "\n".join(editorial))
    wet = [
        "# Next Wet-Lab Validation",
        "",
        "Minimal validation should test model-predicted shifts in fate, growth and diffusion rather than only expression.",
        "",
        "| priority | perturbation | expected readout | assay | decision criterion |",
        "|---:|---|---|---|---|",
        "| 1 | FGF/FGFR modulation in gastruloid or embryo-derived culture | mesoderm/neural fate ratio, MKI67/TOP2A, branch entropy | perturb-and-profile scRNA-seq time course | direction matches SwarmLineage-OT and negative LR control fails |",
        "| 2 | CXCL12/CXCR4 niche-axis perturbation | migration/dispersion and lineage bias | spatial transcriptomics or barcoded live endpoint | spatial/latent dispersion changes in predicted direction |",
        "| 3 | WNT/FZD perturbation | primitive streak/mesoderm branch probability | short time-course scRNA-seq | branch probability and proliferation marker shift |",
        "| 4 | density titration | density-dependent birth/death hazard | controlled aggregate size series | hazard changes with density after cell-cycle adjustment |",
        "| 5 | lineage barcoding follow-up | ancestor-descendant calibration | CellTag/CRISPR barcode with snapshots | OT teacher and simulator improve over interpolation |",
        "",
    ]
    write_text("reports/next_wetlab_validation.md", "\n".join(wet))
    retained = [
        "# Final Retained Results and Methods",
        "",
        "This document keeps only the current retained prototype results and methods. Installation failures and exploratory noise are kept in `logs/` and `reports/negative_results.md` rather than used as claims.",
        "",
        "## Data and Preprocessing",
        "",
        "Default input: `data/processed/cell_level_subset_v1.h5ad`, subsampled and processed into `processed/swarmlineage_input.h5ad`. The dataset contains ordered developmental stages and lineage/cell-type labels. Stage is used as ordered developmental time.",
        "",
        "## OT Teacher",
        "",
        "Adjacent-stage entropic OT couplings are computed through the `run_moscot` entry point. In the current quick prototype, native moscot availability is recorded but an auditable POT/SciPy fallback generates the couplings. Outputs are `processed/ot_teacher.h5ad`, `processed/ot_couplings/*.npz`, and `processed/ot_fate_probabilities.parquet`.",
        "",
        "## SwarmLineage-OT Model",
        "",
        "A minimal PyTorch velocity model is trained to match OT barycentric descendant vectors. Agent simulations add optional birth/death resampling, adaptive diffusion, local swarm rules, CCI modulation and phenomenological memory.",
        "",
        "## Main Quantitative Results",
        "",
        f"Best mean-rank model: `{best}`.",
        "",
        mean.to_markdown(),
        "",
        "## Gate Status",
        "",
        f"Full model passes the quick superiority gate over OT interpolation on at least two core metrics: {supported}.",
        "",
        "## Limitations",
        "",
        "- Current results are computational evidence only.",
        "- Teacher lineage is OT-inferred pseudo-lineage, not true lineage tracing.",
        "- Native moscot/WOT/TIGON/TrajectoryNet/MIOFlow/CellRank2 baselines remain required for high-impact claims.",
        "- Perturbation predictions are exploratory until matched perturbation time-series or wet-lab validation is performed.",
        "",
        "## Reproducibility",
        "",
        "Run `bash reproducibility/run_all.sh` on Unix/WSL systems, or `powershell -ExecutionPolicy Bypass -File .\\reproducibility\\run_all.ps1` on this Windows environment. The PowerShell run was verified locally. Manifest is written to `reproducibility/manifest.json`.",
        "",
    ]
    write_text("manuscript/final_retained_results_and_methods.md", "\n".join(retained))
    manuscript = [
        "# SwarmLineage-OT: lineage-supervised swarm virtual cells from optimal transport pseudo-lineages",
        "",
        "## Abstract",
        "",
        "Optimal transport can infer pseudo-lineage maps from destructive single-cell snapshots, but it does not by itself produce an executable virtual cell population. We introduce SwarmLineage-OT, a prototype lineage-supervised swarm virtual-cell model in which finite cellular agents use local rules, density-dependent birth-death, adaptive diffusion and cell-cell communication to generate developmental trajectories constrained by OT-inferred couplings.",
        "",
        "## Results",
        "",
        "We audited local single-cell resources, built a real-data stage-based OT teacher, trained a minimal finite-agent simulator, and evaluated held-out stage reconstruction across eleven ablations and negative controls. The current quick-run result is a research prototype, not a final Nature-level result.",
        "",
        "## Discussion",
        "",
        "The central contribution is the conversion of OT pseudo-lineage into executable finite-agent supervision. The main limitation is that validation remains computational and relies on fallback OT couplings in the quick run.",
        "",
    ]
    write_text("manuscript/manuscript.md", "\n".join(manuscript))
    methods = [
        "# Methods",
        "",
        "See `manuscript/theory.md` and `manuscript/methods_theory.tex` for the mathematical formulation. The implemented pipeline uses AnnData preprocessing, PCA latent states, adjacent-stage entropic OT, PyTorch velocity fitting, and agent-level ablations.",
        "",
        "Statistical summaries use paired seed-level comparisons, bootstrap confidence intervals and Benjamini-Hochberg correction where applicable.",
        "",
    ]
    write_text("manuscript/methods.md", "\n".join(methods))
    supplementary = [
        "# Supplementary",
        "",
        "- Data audit: `reports/data_audit.md`",
        "- Literature positioning: `reports/literature_positioning.md`",
        "- OT teacher report: `reports/ot_teacher_report.md`",
        "- Ablation interpretation: `reports/ablation_interpretation.md`",
        "- Negative results: `reports/negative_results.md`",
        "- Perturbation evaluation: `reports/perturbation_evaluation.md`",
        "",
    ]
    write_text("manuscript/supplementary.md", "\n".join(supplementary))
    print({"best_model": best, "full_model_gate_supported": supported})


if __name__ == "__main__":
    main()
