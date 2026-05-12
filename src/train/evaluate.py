from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from src.discovery import birth_death_law, branch_nucleation, cci_branch_bias, diffusion_law, memory_hysteresis, phase_diagram
from src.utils.config import ensure_dir, load_config, write_text


CORE = ["sinkhorn", "mmd_rbf", "energy", "celltype_composition_rmse"]
OT_REFERENCE = "M0b_ot_interpolation"
PRIMARY_AGENT = "M9_full_memory"


def _best_model(metrics: pd.DataFrame) -> str:
    score = metrics.groupby("model")[CORE].mean()
    ranks = score.rank(axis=0, ascending=True).mean(axis=1)
    return str(ranks.sort_values().index[0])


def _output_paths(quick_fixture: bool) -> tuple[Path, Path]:
    table_dir = Path("tables/quick_fixture") if quick_fixture else Path("tables")
    report_dir = Path("reports/quick_fixture") if quick_fixture else Path("reports")
    ensure_dir(table_dir)
    ensure_dir(report_dir)
    return table_dir, report_dir


def _configure(config_path: str, quick_fixture: bool) -> tuple[dict, dict]:
    cfg = load_config(config_path)
    model_cfg = load_config(cfg.get("model_config", "configs/model.yaml"))
    if quick_fixture:
        cfg = dict(cfg)
        model_cfg = dict(model_cfg)
        model_cfg["metrics_path"] = "tables/quick_fixture/final_metrics.csv"
        model_cfg["event_log_path"] = "tables/quick_fixture/birth_death_event_log.csv"
        model_cfg["teacher_path"] = "processed/quick_fixture/ot_teacher.h5ad"
        model_cfg["model_dir"] = "results/quick_fixture/models"
        cfg["baseline_execution_matrix_path"] = "reports/quick_fixture/baseline_execution_matrix.csv"
        cfg["module_contribution_report"] = "reports/quick_fixture/module_contribution_audit.md"
        cfg["scientific_gap_report"] = "reports/quick_fixture/scientific_gap_audit.md"
    return cfg, model_cfg


def _event_count_stability(event_path: str | Path, models: list[str]) -> dict[str, float]:
    path = Path(event_path)
    if not path.exists():
        return {model: np.nan for model in models}
    events = pd.read_csv(path)
    if events.empty or not {"seed", "variant"}.issubset(events.columns):
        return {model: np.nan for model in models}
    counts = events.groupby(["variant", "seed"]).size().rename("n_events").reset_index()
    out: dict[str, float] = {}
    all_seeds = sorted(events["seed"].unique())
    for model in models:
        model_counts = counts[counts["variant"] == model].set_index("seed")["n_events"].reindex(all_seeds, fill_value=0)
        mean = float(model_counts.mean())
        out[model] = 0.0 if mean == 0.0 else float(model_counts.std(ddof=0) / max(mean, 1e-8))
    return out


def _teacher_fidelity(metrics: pd.DataFrame, model_cfg: dict, table_dir: Path, report_dir: Path) -> tuple[pd.DataFrame, bool]:
    mean = metrics.groupby("model")[["sinkhorn", "mmd_rbf", "energy", "knn_two_sample_accuracy", "celltype_composition_rmse"]].mean()
    if OT_REFERENCE not in mean.index:
        raise ValueError(f"Missing required OT reference row `{OT_REFERENCE}` in metrics.")
    reference = mean.loc[OT_REFERENCE]
    stability = _event_count_stability(model_cfg.get("event_log_path", "tables/birth_death_event_log.csv"), list(mean.index))
    rows = []
    composition_tolerance = float(model_cfg.get("composition_rmse_tolerance", 0.08))
    for model, row in mean.iterrows():
        is_negative_control = str(model).startswith("M10") or str(model).startswith("M11")
        branch_fate_error = float(row["celltype_composition_rmse"])
        manifold_escape_rate = float(np.clip(2.0 * (row["knn_two_sample_accuracy"] - 0.5), 0.0, 1.0))
        relative_sinkhorn = float(row["sinkhorn"] / max(reference["sinkhorn"], 1e-12))
        relative_mmd = float(row["mmd_rbf"] / max(reference["mmd_rbf"], 1e-12))
        fate_calibration = float(np.clip(1.0 - branch_fate_error / max(composition_tolerance, 1e-12), 0.0, 1.0))
        event_stability = stability.get(str(model), np.nan)
        gate = bool(
            (not is_negative_control)
            and model == OT_REFERENCE
            or (
                (not is_negative_control)
                and
                relative_sinkhorn <= float(model_cfg.get("teacher_fidelity_sinkhorn_ratio_max", 1.75))
                and relative_mmd <= float(model_cfg.get("teacher_fidelity_mmd_ratio_max", 4.0))
                and branch_fate_error <= composition_tolerance
                and manifold_escape_rate <= float(model_cfg.get("manifold_escape_rate_max", 0.45))
                and (np.isnan(event_stability) or event_stability <= float(model_cfg.get("event_count_stability_max", 2.0)))
            )
        )
        rows.append(
            {
                "model": model,
                "reference_model": OT_REFERENCE,
                "control_role": "negative_control" if is_negative_control else "evaluated_model",
                "mean_sinkhorn": float(row["sinkhorn"]),
                "ot_reference_sinkhorn": float(reference["sinkhorn"]),
                "relative_sinkhorn_to_ot_reference": relative_sinkhorn,
                "mean_mmd_rbf": float(row["mmd_rbf"]),
                "ot_reference_mmd_rbf": float(reference["mmd_rbf"]),
                "relative_mmd_to_ot_reference": relative_mmd,
                "composition_rmse": branch_fate_error,
                "composition_rmse_tolerance": composition_tolerance,
                "composition_rmse_pass": bool(branch_fate_error <= composition_tolerance),
                "branch_fate_error": branch_fate_error,
                "manifold_escape_rate": manifold_escape_rate,
                "event_count_stability": event_stability,
                "fate_probability_calibration": fate_calibration,
                "teacher_fidelity_gate": gate,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(table_dir / "teacher_fidelity_metrics.csv", index=False)
    primary = out[out["model"] == PRIMARY_AGENT]
    if primary.empty:
        primary_gate = False
        primary_line = f"- Primary agent `{PRIMARY_AGENT}` was not present in the metrics table."
    else:
        primary_gate = bool(primary["teacher_fidelity_gate"].iloc[0])
        p = primary.iloc[0]
        primary_line = (
            f"- Primary agent `{PRIMARY_AGENT}`: relative_sinkhorn={p['relative_sinkhorn_to_ot_reference']:.3f}, "
            f"relative_mmd={p['relative_mmd_to_ot_reference']:.3f}, branch_fate_error={p['branch_fate_error']:.4f}, "
            f"manifold_escape_rate={p['manifold_escape_rate']:.3f}, gate={primary_gate}."
        )
    report = [
        "# Teacher Fidelity Audit",
        "",
        "`M0b_ot_interpolation` is treated as an oracle-like OT teacher/reference interpolation, not as a competitor that the finite-agent model must beat.",
        "",
        primary_line,
        "",
        "## Fidelity Metrics",
        "",
        out.to_markdown(index=False),
        "",
        "## Interpretation Rule",
        "",
        "- If the agent model trails the OT reference but remains within tolerance, the teacher fidelity gate can pass.",
        "- If the agent model deviates substantially from the reference, this is reported as insufficient teacher fidelity, not as failure to beat OT.",
        "- High-level biological claims still require native moscot/WOT or external validation of the teacher.",
        "",
    ]
    write_text(report_dir / "teacher_fidelity_audit.md", "\n".join(report))
    return out, primary_gate


def _run_discovery(config: str, quick_fixture: bool, table_dir: Path, report_dir: Path) -> tuple[pd.DataFrame, bool]:
    modules = [
        diffusion_law.run,
        birth_death_law.run,
        branch_nucleation.run,
        memory_hysteresis.run,
        cci_branch_bias.run,
        phase_diagram.run,
    ]
    rows = []
    for run in modules:
        try:
            result = run(config, quick_fixture)
            rows.append({"law": result.get("law"), "gate": bool(result.get("gate")), "table": result.get("table"), "status": "executed"})
        except Exception as exc:  # discovery failures are recorded, not hidden
            rows.append({"law": getattr(run, "__module__", "unknown").split(".")[-1], "gate": False, "table": "", "status": f"{type(exc).__name__}: {exc}"})
    summary = pd.DataFrame(rows)
    summary.to_csv(table_dir / "emergent_law_gate_summary.csv", index=False)
    n_pass = int(summary["gate"].sum()) if "gate" in summary else 0
    emergent_gate = bool(n_pass >= 3)
    report = [
        "# Emergent Law Summary",
        "",
        "The core scientific objective is to learn microscopic finite-agent rules that realize the OT-inferred pseudo-lineage and expose laws that OT interpolation does not directly return.",
        "",
        f"- passed discovery gates: {n_pass}/{summary.shape[0]}",
        f"- emergent_law_gate: {emergent_gate}",
        "",
        "## Discovery Gates",
        "",
        summary.to_markdown(index=False),
        "",
        "## Result Interpretation",
        "",
        "- These analyses are mechanistic probes of the trained simulator, not wet-lab validation.",
        "- Stable laws can support the research prototype even when the OT reference has the lowest reconstruction error.",
        "- If both teacher fidelity and all emergent-law gates fail, the current scientific hypothesis is not supported.",
        "",
    ]
    write_text(report_dir / "emergent_law_summary.md", "\n".join(report))
    return summary, emergent_gate


def _write_module_report(metrics: pd.DataFrame, fidelity: pd.DataFrame, laws: pd.DataFrame, mechanistic_gate: bool, cfg: dict, report_dir: Path) -> None:
    mean = metrics.groupby("model")[CORE].mean()
    primary_fidelity = fidelity[fidelity["model"] == PRIMARY_AGENT]
    if primary_fidelity.empty:
        primary_text = f"`{PRIMARY_AGENT}` missing."
    else:
        row = primary_fidelity.iloc[0]
        primary_text = (
            f"`{PRIMARY_AGENT}` fidelity: relative_sinkhorn={row['relative_sinkhorn_to_ot_reference']:.3f}, "
            f"relative_mmd={row['relative_mmd_to_ot_reference']:.3f}, branch_fate_error={row['branch_fate_error']:.4f}."
        )
    module_report = [
        "# Module Contribution Audit",
        "",
        "This audit no longer treats OT interpolation as a method that the agent model must beat. OT interpolation is the teacher/reference map.",
        "",
        f"- primary agent: `{PRIMARY_AGENT}`",
        f"- {primary_text}",
        f"- mechanistic_usefulness_gate: {mechanistic_gate}",
        "",
        "## Reconstruction Context",
        "",
        mean.to_markdown(),
        "",
        "## Emergent-Law Contributions",
        "",
        laws.to_markdown(index=False),
        "",
        "## Interpretation",
        "",
        "- `M0b_ot_interpolation` can remain the lowest-error reconstruction because it is the reference interpolation derived from the teacher.",
        "- Swarm, birth/death, diffusion, CCI and memory are evaluated by fidelity plus whether they expose robust order parameters, hazards, phase regimes or perturbation sensitivities.",
        "- Modules that degrade fidelity without yielding stable laws should remain exploratory or be disabled in retained biological claims.",
        "",
    ]
    write_text(cfg.get("module_contribution_report", str(report_dir / "module_contribution_audit.md")), "\n".join(module_report))


def _write_scientific_gap_report(
    metrics: pd.DataFrame,
    fidelity_gate: bool,
    emergent_gate: bool,
    mechanistic_gate: bool,
    cfg: dict,
    report_dir: Path,
) -> None:
    baseline_path = cfg.get("baseline_execution_matrix_path", "reports/baseline_execution_matrix.csv")
    baseline_matrix = pd.read_csv(baseline_path) if Path(baseline_path).exists() else pd.DataFrame()
    best = _best_model(metrics)
    report = [
        "# Scientific Gap Audit",
        "",
        "Scientific goal after re-scoping: OT gives the developmental map; SwarmLineage-OT learns microscopic rules that realize the map and reveal emergent developmental laws.",
        "",
        f"- best mean-rank reconstruction row: `{best}`",
        f"- OT reference row: `{OT_REFERENCE}`",
        f"- teacher_fidelity_gate: {fidelity_gate}",
        f"- emergent_law_gate: {emergent_gate}",
        f"- mechanistic_usefulness_gate: {mechanistic_gate}",
        "",
        "## Baseline/Reference Execution Matrix",
        "",
        baseline_matrix.to_markdown(index=False) if not baseline_matrix.empty else "No baseline execution matrix found.",
        "",
        "## Remaining Gaps",
        "",
        "- Native moscot TemporalProblem or externally validated teacher is still required for strong claims beyond a toy fallback.",
        "- Emergent laws must be checked across seeds, held-out times and at least one external developmental dataset.",
        "- CCI and memory laws are computational hypotheses unless supported by matched spatial, perturbation or lineage-tracing data.",
        "- The manuscript must not state that SwarmLineage-OT outperforms OT; the correct claim is teacher fidelity plus mechanistic discovery.",
        "",
    ]
    write_text(cfg.get("scientific_gap_report", str(report_dir / "scientific_gap_audit.md")), "\n".join(report))


def _write_retained_methods(
    metrics: pd.DataFrame,
    fidelity_gate: bool,
    emergent_gate: bool,
    mechanistic_gate: bool,
    table_dir: Path,
    quick_fixture: bool,
) -> None:
    fidelity_path = table_dir / "teacher_fidelity_metrics.csv"
    law_path = table_dir / "emergent_law_gate_summary.csv"
    retained = [
        "# Final Retained Results and Methods",
        "",
        "## Central Claim",
        "",
        "OT gives the developmental map; SwarmLineage-OT learns microscopic finite-agent rules that realize the map and reveal emergent developmental laws.",
        "",
        "The retained manuscript must not claim that SwarmLineage-OT outperforms OT interpolation. `M0b_ot_interpolation` is an oracle-like teacher/reference interpolation.",
        "",
        "## Evaluation Gates",
        "",
        f"- teacher_fidelity_gate: {fidelity_gate}",
        f"- emergent_law_gate: {emergent_gate}",
        f"- mechanistic_usefulness_gate: {mechanistic_gate}",
        "",
        "## Retained Metrics",
        "",
        f"- teacher fidelity metrics: `{fidelity_path.as_posix()}`",
        f"- emergent law gates: `{law_path.as_posix()}`",
        "",
        metrics.groupby("model")[CORE].mean().to_markdown(),
        "",
        "## Methods Retained",
        "",
        "- strict time holdout and teacher-edge holdout support leakage-resistant evaluation.",
        "- OT teacher construction records backend status; toy fallback is not presented as native moscot/WOT.",
        "- `SwarmLineageDynamics` represents trainable intrinsic, teacher, swarm, birth/death, adaptive diffusion, CCI and memory components.",
        "- Stochastic birth/death uses event simulation and writes event logs.",
        "- Discovery modules estimate diffusion, growth, branch nucleation, memory hysteresis, CCI branch bias and phase-regime laws.",
        "",
        "## Interpretation",
        "",
        "If `M0b_ot_interpolation` has the lowest reconstruction error, this is expected for a teacher/reference. The agent model is retained when it stays close enough to the teacher and yields stable mechanistic laws.",
        "",
        "## Limitations",
        "",
        "- Current emergent laws are computational hypotheses, not validated biological mechanisms.",
        "- Strong biological claims require native moscot/WOT or external teacher validation plus external or perturbation validation.",
        "- If teacher fidelity is poor and no emergent-law gates are stable, the current scientific hypothesis should be reported as unsupported.",
        "",
    ]
    out_path = "manuscript/final_retained_results_and_methods.quick_fixture.md" if quick_fixture else "manuscript/final_retained_results_and_methods.md"
    write_text(out_path, "\n".join(retained))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/train.yaml")
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    cfg, model_cfg = _configure(args.config, args.quick_fixture)
    table_dir, report_dir = _output_paths(args.quick_fixture)
    metrics_path = model_cfg.get("metrics_path", "tables/final_metrics.csv")
    metrics = pd.read_csv(metrics_path)
    mean = metrics.groupby("model")[["sinkhorn", "mmd_rbf", "energy", "knn_two_sample_accuracy", "celltype_composition_rmse"]].agg(["mean", "std"])
    mean.to_csv(table_dir / "final_model_metric_summary.csv")

    fidelity, teacher_fidelity_gate = _teacher_fidelity(metrics, model_cfg, table_dir, report_dir)
    laws, emergent_law_gate = _run_discovery(args.config, args.quick_fixture, table_dir, report_dir)
    primary = fidelity[fidelity["model"] == PRIMARY_AGENT]
    catastrophic_divergence = True
    if not primary.empty:
        p = primary.iloc[0]
        catastrophic_divergence = bool(
            p["relative_sinkhorn_to_ot_reference"] > 2.5
            or p["relative_mmd_to_ot_reference"] > 8.0
            or p["manifold_escape_rate"] > 0.75
        )
    mechanistic_usefulness_gate = bool((not catastrophic_divergence) and (teacher_fidelity_gate or emergent_law_gate) and int(laws["gate"].sum()) >= 2)

    _write_module_report(metrics, fidelity, laws, mechanistic_usefulness_gate, cfg, report_dir)
    _write_scientific_gap_report(metrics, teacher_fidelity_gate, emergent_law_gate, mechanistic_usefulness_gate, cfg, report_dir)
    _write_retained_methods(metrics, teacher_fidelity_gate, emergent_law_gate, mechanistic_usefulness_gate, table_dir, args.quick_fixture)

    editorial = [
        "# Editorial Assessment",
        "",
        "Current evidence level: computational research prototype.",
        "",
        f"- teacher_fidelity_gate: {teacher_fidelity_gate}",
        f"- emergent_law_gate: {emergent_law_gate}",
        f"- mechanistic_usefulness_gate: {mechanistic_usefulness_gate}",
        "",
        "Not sufficient for Nature/Nature Methods/Nature Biotechnology until native or externally validated teacher results, external validation and stronger biological evidence are added.",
        "",
    ]
    write_text(report_dir / "editorial_assessment.md", "\n".join(editorial))

    print(
        {
            "ot_reference": OT_REFERENCE,
            "primary_agent": PRIMARY_AGENT,
            "teacher_fidelity_gate": teacher_fidelity_gate,
            "emergent_law_gate": emergent_law_gate,
            "mechanistic_usefulness_gate": mechanistic_usefulness_gate,
        }
    )


if __name__ == "__main__":
    main()
