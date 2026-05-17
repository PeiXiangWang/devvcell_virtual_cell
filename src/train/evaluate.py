from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.discovery import birth_death_law, branch_nucleation, cci_branch_bias, diffusion_law, memory_hysteresis, phase_diagram
from src.discovery.common import TIER_ORDER, bh_q_values, configure_paths, output_dirs, tier_at_least
from src.utils.config import ensure_dir, load_config, write_text


CORE = ["sinkhorn", "mmd_rbf", "energy", "celltype_composition_rmse"]
OT_REFERENCE = "M0b_ot_interpolation"
PRIMARY_AGENT = "M5_ot_swarm"


def _best_model(metrics: pd.DataFrame) -> str:
    score = metrics.groupby("model")[CORE].mean()
    ranks = score.rank(axis=0, ascending=True).mean(axis=1)
    return str(ranks.sort_values().index[0])


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


def _metric_tier(row: pd.Series, thresholds: dict) -> tuple[str, list[str]]:
    failure_reasons: list[str] = []
    for tier in ("strong", "acceptable", "weak"):
        t = thresholds.get(tier, {})
        checks = {
            "relative_sinkhorn_to_ot_reference": row["relative_sinkhorn_to_ot_reference"] <= float(t.get("relative_sinkhorn_max", np.inf)),
            "relative_mmd_to_ot_reference": row["relative_mmd_to_ot_reference"] <= float(t.get("relative_mmd_max", np.inf)),
            "relative_energy_to_ot_reference": row["relative_energy_to_ot_reference"] <= float(t.get("relative_energy_max", np.inf)),
            "composition_rmse": row["composition_rmse"] <= float(t.get("composition_rmse_max", np.inf)),
            "manifold_escape_rate": row["manifold_escape_rate"] <= float(t.get("manifold_escape_rate_max", np.inf)),
        }
        if all(checks.values()):
            return tier, failure_reasons
        if tier == "weak":
            for key, passed in checks.items():
                if not passed:
                    failure_reasons.append(f"{key}_above_{tier}_threshold")
    return "fail", failure_reasons or ["above_weak_threshold"]


def _teacher_fidelity(metrics: pd.DataFrame, model_cfg: dict, discovery_cfg: dict, table_dir: Path, report_dir: Path, fig_dir: Path) -> tuple[pd.DataFrame, str]:
    mean = metrics.groupby("model")[["sinkhorn", "mmd_rbf", "energy", "knn_two_sample_accuracy", "celltype_composition_rmse"]].mean()
    if OT_REFERENCE not in mean.index:
        raise ValueError(f"Missing required OT reference row `{OT_REFERENCE}` in metrics.")
    reference = mean.loc[OT_REFERENCE]
    stability = _event_count_stability(model_cfg.get("event_log_path", "tables/birth_death_event_log.csv"), list(mean.index))
    thresholds = discovery_cfg.get("teacher_fidelity", {})
    rows = []
    for model, metric_row in mean.iterrows():
        is_negative_control = str(model).startswith("M10") or str(model).startswith("M11")
        branch_fate_error = float(metric_row["celltype_composition_rmse"])
        row = {
            "model": model,
            "reference_model": OT_REFERENCE,
            "control_role": "negative_control" if is_negative_control else "evaluated_model",
            "negative_control_expected_to_fail": bool(is_negative_control),
            "mean_sinkhorn": float(metric_row["sinkhorn"]),
            "ot_reference_sinkhorn": float(reference["sinkhorn"]),
            "relative_sinkhorn_to_ot_reference": float(metric_row["sinkhorn"] / max(reference["sinkhorn"], 1e-12)),
            "mean_mmd_rbf": float(metric_row["mmd_rbf"]),
            "ot_reference_mmd_rbf": float(reference["mmd_rbf"]),
            "relative_mmd_to_ot_reference": float(metric_row["mmd_rbf"] / max(reference["mmd_rbf"], 1e-12)),
            "mean_energy": float(metric_row["energy"]),
            "ot_reference_energy": float(reference["energy"]),
            "relative_energy_to_ot_reference": float(metric_row["energy"] / max(reference["energy"], 1e-12)),
            "composition_rmse": branch_fate_error,
            "branch_fate_error": branch_fate_error,
            "manifold_escape_rate": float(np.clip(2.0 * (metric_row["knn_two_sample_accuracy"] - 0.5), 0.0, 1.0)),
            "event_count_stability": stability.get(str(model), np.nan),
            "fate_probability_calibration": float(np.clip(1.0 - branch_fate_error / max(thresholds.get("weak", {}).get("composition_rmse_max", 0.05), 1e-12), 0.0, 1.0)),
        }
        numeric_tier, reasons = _metric_tier(pd.Series(row), thresholds)
        if is_negative_control:
            expected_met = (
                row["relative_sinkhorn_to_ot_reference"] > 1.05
                or row["relative_mmd_to_ot_reference"] > 1.20
                or row["branch_fate_error"] > float(thresholds.get("strong", {}).get("composition_rmse_max", 0.005))
            )
            row["negative_control_outcome"] = "fail_expectation_met" if expected_met else "fail_expectation_not_met"
            row["numeric_teacher_fidelity_tier"] = numeric_tier
            row["teacher_fidelity_tier"] = "fail"
            row["failure_reasons"] = "negative_control_expected_to_fail"
        else:
            row["negative_control_outcome"] = "not_applicable"
            row["numeric_teacher_fidelity_tier"] = numeric_tier
            row["teacher_fidelity_tier"] = numeric_tier
            row["failure_reasons"] = ";".join(reasons)
        row["gate_pass"] = row["teacher_fidelity_tier"] in {"acceptable", "strong"}
        row["strong_gate"] = row["teacher_fidelity_tier"] == "strong"
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(table_dir / "teacher_fidelity_metrics.csv", index=False)

    fig, ax = plt.subplots(figsize=(8, 4.3), dpi=160)
    colors = {"strong": "#2f7d32", "acceptable": "#6d9dc5", "weak": "#d9a441", "fail": "#b65d5d"}
    plot = out.sort_values("relative_sinkhorn_to_ot_reference")
    ax.bar(plot["model"], plot["relative_sinkhorn_to_ot_reference"], color=[colors.get(t, "#999") for t in plot["teacher_fidelity_tier"]])
    ax.axhline(1.0, color="black", lw=0.8)
    ax.set_ylabel("relative Sinkhorn to OT reference")
    ax.set_title("Teacher fidelity tiers")
    ax.tick_params(axis="x", rotation=70, labelsize=6)
    fig.tight_layout()
    fig.savefig(fig_dir / "teacher_fidelity_tiers.png")
    plt.close(fig)

    primary = out[out["model"] == PRIMARY_AGENT]
    primary_tier = str(primary["teacher_fidelity_tier"].iloc[0]) if not primary.empty else "fail"
    primary_line = "Primary agent missing." if primary.empty else (
        f"Primary `{PRIMARY_AGENT}` tier: {primary_tier}; "
        f"relative_sinkhorn={float(primary['relative_sinkhorn_to_ot_reference'].iloc[0]):.3f}; "
        f"relative_mmd={float(primary['relative_mmd_to_ot_reference'].iloc[0]):.3f}; "
        f"composition={float(primary['composition_rmse'].iloc[0]):.4f}."
    )
    write_text(
        report_dir / "teacher_fidelity_audit.md",
        "\n".join(
            [
                "# Teacher Fidelity Audit",
                "",
                "`M0b_ot_interpolation` is an oracle-like OT teacher/reference interpolation. The finite-agent model is evaluated by teacher fidelity, emergent-law robustness and mechanistic usefulness, not by beating the OT reference.",
                "",
                f"- {primary_line}",
                "- Negative controls are marked as expected-to-fail controls and cannot pass as normal evaluated models.",
                "",
                "## Metrics",
                "",
                out.to_markdown(index=False),
                "",
            ]
        ),
    )
    return out, primary_tier


def _run_discovery(config: str, quick_fixture: bool, table_dir: Path, report_dir: Path) -> tuple[pd.DataFrame, str]:
    modules = [diffusion_law.run, birth_death_law.run, branch_nucleation.run, memory_hysteresis.run, cci_branch_bias.run, phase_diagram.run]
    rows = []
    for run in modules:
        try:
            rows.append(run(config, quick_fixture))
        except Exception as exc:
            rows.append(
                {
                    "law": getattr(run, "__module__", "unknown").split(".")[-1],
                    "tier": "fail",
                    "gate_pass": False,
                    "strong_gate": False,
                    "effect_size": np.nan,
                    "effect_ci_low": np.nan,
                    "effect_ci_high": np.nan,
                    "permutation_p": 1.0,
                    "permutation_q": 1.0,
                    "negative_control_pass": False,
                    "seed_stability_pass": False,
                    "rollout_based": False,
                    "directly_supervised_or_encoded": False,
                    "interpretation_level": "unsupported",
                    "table": "",
                    "report": "",
                    "status": f"{type(exc).__name__}: {exc}",
                }
            )
    summary = pd.DataFrame(rows)
    if not summary.empty:
        summary["permutation_q"] = bh_q_values(summary["permutation_p"].fillna(1.0).tolist())
    summary.to_csv(table_dir / "emergent_law_gate_summary.csv", index=False)
    if (summary["tier"] == "strong").sum() >= 2:
        emergent_tier = "strong"
    elif summary["tier"].isin(["acceptable", "strong"]).sum() >= 3:
        emergent_tier = "acceptable"
    elif summary["tier"].isin(["weak", "acceptable", "strong"]).sum() >= 3:
        emergent_tier = "weak"
    else:
        emergent_tier = "fail"
    write_text(
        report_dir / "emergent_law_summary.md",
        "\n".join(
            [
                "# Emergent Law Summary",
                "",
                "Discovery modules are graded as fail, weak, acceptable or strong using effect size, seed stability, negative controls and rollout support.",
                "",
                f"- emergent_law_tier: {emergent_tier}",
                f"- laws at least acceptable: {int(summary['tier'].isin(['acceptable', 'strong']).sum())}/{summary.shape[0]}",
                f"- strong laws: {int((summary['tier'] == 'strong').sum())}/{summary.shape[0]}",
                "",
                summary.to_markdown(index=False),
                "",
            ]
        ),
    )
    return summary, emergent_tier


def _teacher_backend_status(model_cfg: dict, table_dir: Path, report_dir: Path) -> tuple[pd.DataFrame, bool]:
    candidates = [Path(model_cfg.get("teacher_path", "processed/ot_teacher.h5ad")).with_name("ot_teacher_summary.json"), Path("processed/ot_teacher_summary.json")]
    payload = {}
    for path in candidates:
        if path.exists():
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
                break
            except Exception:
                payload = {}
    backend = str(payload.get("teacher_backend", "unknown"))
    native_validated = backend in {"native_moscot", "native_wot"}
    teacher_path = Path(model_cfg.get("teacher_path", "processed/ot_teacher.h5ad"))
    run_summary_path = Path("processed/quick_fixture/ot_couplings/moscot_run_summary.json") if "quick_fixture" in str(teacher_path) else Path("processed/ot_couplings/moscot_run_summary.json")
    run_summary = {}
    if run_summary_path.exists():
        try:
            run_summary = json.loads(run_summary_path.read_text(encoding="utf-8"))
        except Exception:
            run_summary = {}
    native_pairs = run_summary.get("pairs", [])
    native_status = run_summary.get("native_moscot_status", {})
    table = pd.DataFrame(
        [
            {
                "teacher_backend": backend,
                "native_teacher_available": native_validated,
                "external_teacher_validation": False,
                "native_teacher_claims_allowed": native_validated,
                "strong_biological_claims_allowed": False,
                "nature_level_claim_allowed": False,
                "native_temporalproblem_pairs": len(native_pairs) if native_validated else 0,
                "native_requirements_file": "reproducibility/native_moscot_requirements.txt",
                "status": "native teacher pending" if not native_validated else "native teacher available; external validation still required",
            }
        ]
    )
    table.to_csv(table_dir / "teacher_backend_status.csv", index=False)
    write_text(
        report_dir / "native_teacher_status.md",
        "\n".join(
            [
                "# Native Teacher Status",
                "",
                f"- teacher_backend: {backend}",
                f"- native_teacher_available: {native_validated}",
                "- external_teacher_validation: False",
                f"- native TemporalProblem pairs extracted: {len(native_pairs) if native_validated else 0}",
                f"- native moscot status: {native_status}",
                "- clean native stack: `reproducibility/native_moscot_requirements.txt`",
                "- Native moscot/WOT teacher claims are allowed only when `native_teacher_available=True`.",
                "- Strong biological claims remain forbidden without external, lineage, perturbation or experimental validation.",
                "- Nature-level claims are forbidden for the current prototype.",
                "",
            ]
        ),
    )
    return table, native_validated


def _mechanistic_tier(teacher_tier: str, laws: pd.DataFrame, native_validated: bool, discovery_cfg: dict) -> str:
    acceptable = int(laws["tier"].isin(["acceptable", "strong"]).sum()) if not laws.empty else 0
    strong = int((laws["tier"] == "strong").sum()) if not laws.empty else 0
    req = discovery_cfg.get("mechanistic_usefulness", {})
    strong_req = req.get("strong_requires", {})
    acc_req = req.get("acceptable_requires", {})
    if (
        tier_at_least(teacher_tier, strong_req.get("teacher_fidelity_tier_at_least", "acceptable"))
        and acceptable >= int(strong_req.get("emergent_laws_at_least_acceptable", 4))
        and strong >= int(strong_req.get("emergent_laws_at_least_strong", 2))
        and (native_validated or not bool(strong_req.get("external_or_native_teacher_validation", True)))
    ):
        return "strong"
    if (
        tier_at_least(teacher_tier, acc_req.get("teacher_fidelity_tier_at_least", "acceptable"))
        and acceptable >= int(acc_req.get("emergent_laws_at_least_acceptable", 3))
        and strong >= int(acc_req.get("emergent_laws_at_least_strong", 1))
    ):
        return "acceptable"
    if tier_at_least(teacher_tier, "weak") and laws["tier"].isin(["weak", "acceptable", "strong"]).sum() >= 2:
        return "weak"
    return "fail"


def _clone_boundary_context() -> dict | None:
    """Return the latest clone-aware evidence boundary, preferring the outcome-preserving audit."""
    candidates = [
        ("outcome_preserving", Path("tables/clone_outcome_preserving_native_final_summary.csv")),
        ("clone_stratified", Path("tables/clone_stratified_native_final_summary.csv")),
    ]
    for source, path in candidates:
        if not path.exists():
            continue
        try:
            table = pd.read_csv(path)
        except Exception:
            continue
        if table.empty:
            continue
        row = table.iloc[0].to_dict()
        return {
            "source": source,
            "status": str(row.get("final_clone_aware_status", "unknown")),
            "tier": str(row.get("final_clone_aware_tier", "unknown")),
            "interpretation": str(row.get("interpretation", "No interpretation recorded.")),
        }
    return None


def _clone_boundary_report_lines(ctx: dict | None) -> list[str]:
    if not ctx:
        return [
            "## Clone-Aware Evidence Boundary",
            "",
            "No finalized clone-aware audit table was found. Clone-level fate-diversification support is not retained by default.",
            "",
        ]
    boundary = (
        "Clone-level fate-diversification prediction is not supported under current tested datasets and native sampling strategies."
        if ctx["tier"] in {"fail", "weak"}
        else "Clone-aware support is a computational candidate that still requires independent validation."
    )
    return [
        "## Clone-Aware Evidence Boundary",
        "",
        f"- latest_audit_source: `{ctx['source']}`",
        f"- final_clone_aware_status: `{ctx['status']}`",
        f"- final_clone_aware_tier: `{ctx['tier']}`",
        f"- interpretation: {ctx['interpretation']}",
        f"- boundary: {boundary}",
        "",
        "Retained: branch nucleation / transient condensation-before-divergence as a time-series order-parameter computational hypothesis, supported by internal native moscot and E1 MouseGastrulationData with M5_ot_swarm as the evidence-selected primary model.",
        "",
        "Not retained: clone-level fate-diversification prediction from condensation, topological-neighbour-specific mechanism, swarm-required causality, birth/death, memory, CCI, or diffusion as an independent discovery.",
        "",
        "Clone-aware analyses are stress tests for the time-series branch signature, not the retained main claim. They do not justify presenting condensation as a clone-level predictor unless the primary outcome is supported across datasets after covariate, matched and negative-control analyses.",
        "",
    ]


def _developmental_atlas_context() -> dict | None:
    path = Path("tables/developmental_branch_window_atlas_final_summary.csv")
    if not path.exists():
        return None
    try:
        table = pd.read_csv(path)
    except Exception:
        return None
    if table.empty:
        return None
    row = table.iloc[0].to_dict()
    summary_path = Path("tables/developmental_branch_window_summary.csv")
    registry_path = Path("tables/developmental_time_series_dataset_registry.csv")
    analyzed: list[str] = []
    downloaded: list[str] = []
    independent_native: list[str] = []
    weak_or_failed: list[str] = []
    if summary_path.exists():
        try:
            summary = pd.read_csv(summary_path)
            analyzed = summary.get("dataset_id", pd.Series(dtype=str)).astype(str).tolist()
            independent_native = summary[
                summary.get("independence_tier", pd.Series(dtype=str)).astype(str).str.contains("independent", na=False)
                & summary.get("teacher_backend", pd.Series(dtype=str)).astype(str).eq("native_moscot")
            ].get("dataset_id", pd.Series(dtype=str)).astype(str).tolist()
            weak_or_failed = summary[
                summary.get("external_support_tier", pd.Series(dtype=str)).astype(str).isin(["weak", "fail"])
            ].get("dataset_id", pd.Series(dtype=str)).astype(str).tolist()
        except Exception:
            analyzed = []
    if registry_path.exists():
        try:
            registry = pd.read_csv(registry_path)
            downloaded = registry[
                registry.get("download_success", pd.Series(dtype=bool)).astype(str).str.lower().eq("true")
                & registry.get("selected_for_analysis", pd.Series(dtype=bool)).astype(str).str.lower().eq("true")
            ].get("dataset_id", pd.Series(dtype=str)).astype(str).tolist()
        except Exception:
            downloaded = []
    return {
        "tier": str(row.get("developmental_branch_window_overall_tier", "unknown")),
        "interpretation": str(row.get("interpretation", "No interpretation recorded.")),
        "new_datasets_attempted": int(row.get("new_datasets_attempted", 0)),
        "new_datasets_analyzed": int(row.get("new_datasets_analyzed", 0)),
        "acceptable_external_datasets": int(row.get("acceptable_external_datasets", 0)),
        "v2_final_sprint_directions_attempted": int(row.get("v2_final_sprint_directions_attempted", 0)),
        "v2_final_sprint_usable_datasets": int(row.get("v2_final_sprint_usable_datasets", 0)),
        "v2_final_sprint_acceptable_datasets": int(row.get("v2_final_sprint_acceptable_datasets", 0)),
        "v2_spatial_validation_status": str(row.get("v2_spatial_validation_status", "")),
        "final_manuscript_line": str(row.get("final_manuscript_line", "")),
        "analyzed_datasets": analyzed,
        "downloaded_datasets": downloaded,
        "independent_native_analyzed": independent_native,
        "weak_or_failed_datasets": weak_or_failed,
    }


def _developmental_atlas_report_lines(ctx: dict | None) -> list[str]:
    if not ctx:
        return [
            "## Developmental Time-Series Atlas",
            "",
            "No finalized developmental branch-window atlas summary was found. The retained time-series claim remains based on internal native moscot and E1 support only.",
            "",
        ]
    v2_lines = []
    if int(ctx.get("v2_final_sprint_directions_attempted", 0)) > 0:
        v2_lines = [
            f"- v2_final_sprint_directions_attempted: {ctx['v2_final_sprint_directions_attempted']}",
            f"- v2_final_sprint_usable_datasets: {ctx['v2_final_sprint_usable_datasets']}",
            f"- v2_final_sprint_acceptable_datasets: {ctx['v2_final_sprint_acceptable_datasets']}",
            f"- v2_spatial_validation_status: {ctx['v2_spatial_validation_status'] or 'not_recorded'}",
            f"- final_manuscript_line: {ctx['final_manuscript_line'] or 'not_recorded'}",
        ]
    return [
        "## Developmental Time-Series Atlas",
        "",
        f"- atlas_tier: `{ctx['tier']}`",
        f"- datasets_attempted: {ctx['new_datasets_attempted']}",
        f"- datasets_analyzed: {ctx['new_datasets_analyzed']}",
        f"- acceptable_external_datasets: {ctx['acceptable_external_datasets']}",
        f"- analyzed_datasets: {', '.join(ctx.get('analyzed_datasets', [])) or 'none'}",
        f"- downloaded_new_dataset: {', '.join(ctx.get('downloaded_datasets', [])) or 'none'}",
        f"- independent_native_analyzed: {', '.join(ctx.get('independent_native_analyzed', [])) or 'none'}",
        f"- interpretation: {ctx['interpretation']}",
        *v2_lines,
        "",
        "The atlas is used to define the current external boundary of the branch-window order-parameter hypothesis. Weak or failed atlas rows must not be written as cross-dataset validation. A detected branch-like window without condensation-before-divergence, or with unclean controls, does not upgrade the retained claim.",
        "",
    ]


def _grn_context() -> dict | None:
    path = Path("tables/grn_evidence_matrix.csv")
    if not path.exists():
        return None
    try:
        table = pd.read_csv(path)
    except Exception:
        return None
    if table.empty:
        return None
    grn = table[table["claim"].astype(str).str.contains("GRN/regulon", case=False, na=False)]
    known = table[table["claim"].astype(str).str.contains("Known developmental", case=False, na=False)]
    return {
        "grn_tier": str(grn.iloc[0].get("tier", "unknown")) if not grn.empty else "unknown",
        "grn_allowed": str(grn.iloc[0].get("allowed_language", "No GRN interpretation recorded.")) if not grn.empty else "No GRN interpretation recorded.",
        "known_tier": str(known.iloc[0].get("tier", "unknown")) if not known.empty else "unknown",
    }


def _grn_report_lines(ctx: dict | None) -> list[str]:
    if not ctx:
        return [
            "## GRN / Regulon Evidence Boundary",
            "",
            "No finalized GRN evidence matrix was found. Regulatory mechanism support is not retained by default.",
            "",
        ]
    return [
        "## GRN / Regulon Evidence Boundary",
        "",
        f"- final_GRN_tier: `{ctx['grn_tier']}`",
        f"- known_TF_program_recovery_tier: `{ctx['known_tier']}`",
        f"- interpretation: {ctx['grn_allowed']}",
        "- boundary: GRN/regulon analysis is a computational audit and candidate-generation layer. It does not establish causal GRN control, validated TF perturbation, experimental validation, or a proven regulatory mechanism.",
        "",
    ]


def _write_reports(
    metrics: pd.DataFrame,
    fidelity: pd.DataFrame,
    teacher_tier: str,
    laws: pd.DataFrame,
    emergent_tier: str,
    mechanistic_tier: str,
    native_validated: bool,
    cfg: dict,
    table_dir: Path,
    report_dir: Path,
    quick_fixture: bool,
) -> None:
    mean = metrics.groupby("model")[CORE].mean()
    clone_ctx = None if quick_fixture else _clone_boundary_context()
    clone_lines = [] if quick_fixture else _clone_boundary_report_lines(clone_ctx)
    atlas_ctx = None if quick_fixture else _developmental_atlas_context()
    atlas_lines = [] if quick_fixture else _developmental_atlas_report_lines(atlas_ctx)
    grn_ctx = None if quick_fixture else _grn_context()
    grn_lines = [] if quick_fixture else _grn_report_lines(grn_ctx)
    mech = pd.DataFrame(
        [
            {
                "teacher_fidelity_tier": teacher_tier,
                "emergent_law_tier": emergent_tier,
                "mechanistic_usefulness_tier": mechanistic_tier,
                "mechanistic_gate_pass": mechanistic_tier in {"acceptable", "strong"},
                "strong_gate": mechanistic_tier == "strong",
                "laws_at_least_acceptable": int(laws["tier"].isin(["acceptable", "strong"]).sum()),
                "laws_strong": int((laws["tier"] == "strong").sum()),
                "native_or_external_teacher_validation": bool(native_validated),
            }
        ]
    )
    mech.to_csv(table_dir / "mechanistic_usefulness_summary.csv", index=False)
    retained = laws[laws["interpretation_level"].isin(["retained_computational_hypothesis", "rollout_supported_mechanistic_probe", "encoded_control_law_recovery"])]
    exploratory = laws[laws["interpretation_level"].isin(["exploratory_sensitivity", "demonstration_only"])]
    unsupported = laws[laws["interpretation_level"] == "unsupported"]

    write_text(
        cfg.get("module_contribution_report", str(report_dir / "module_contribution_audit.md")),
        "\n".join(
            [
                "# Module Contribution Audit",
                "",
                "OT interpolation is the teacher/reference map, not a competitor to beat.",
                "",
                f"- teacher_fidelity_tier: {teacher_tier}",
                f"- emergent_law_tier: {emergent_tier}",
                f"- mechanistic_usefulness_tier: {mechanistic_tier}",
                "",
                "## Reconstruction Context",
                "",
                mean.to_markdown(),
                "",
                "## Discovery Tiers",
                "",
                laws.to_markdown(index=False),
                "",
            ]
        ),
    )
    write_text(
        cfg.get("scientific_gap_report", str(report_dir / "scientific_gap_audit.md")),
        "\n".join(
            [
                "# Scientific Gap Audit",
                "",
                "OT gives the developmental map; SwarmLineage-OT converts it into finite-agent dynamics and audits branch-window order parameters in developmental time-series data.",
                "",
                f"- best mean-rank reconstruction row: `{_best_model(metrics)}`",
                f"- OT reference row: `{OT_REFERENCE}`",
                f"- teacher_fidelity_tier: {teacher_tier}",
                f"- emergent_law_tier: {emergent_tier}",
                f"- mechanistic_usefulness_tier: {mechanistic_tier}",
                f"- native_or_external_teacher_validation: {native_validated}",
                "",
                "## Remaining Gaps",
                "",
                "- Native moscot/WOT or external teacher validation status is reported in `teacher_backend_status.csv`.",
                "- CCI and memory are computational probes; experimental validation is absent.",
                "- No manuscript claim may frame SwarmLineage-OT as surpassing the OT reference.",
                "- High-impact submission claims remain unsupported without external lineage, perturbation or experimental validation.",
                "" if quick_fixture else "- Clone-level fate-diversification prediction is not retained; the current main line remains a time-series order-parameter hypothesis.",
                "",
            ]
            + ([] if quick_fixture else _clone_boundary_report_lines(clone_ctx))
            + ([] if quick_fixture else _developmental_atlas_report_lines(atlas_ctx))
            + ([] if quick_fixture else _grn_report_lines(grn_ctx))
        ),
    )
    final_lines = [
        "# Final Retained Results and Methods",
        "",
        "## Central Claim",
        "",
        "SwarmLineage-OT converts native OT-inferred developmental maps into finite-agent virtual-cell dynamics and currently retains a branch-window order-parameter hypothesis, transient condensation-before-divergence, rather than a clone-fate, CCI, memory, birth/death or topological-specific mechanism claim.",
        "",
        "`M0b_ot_interpolation` is an oracle-like OT teacher/reference interpolation. The finite-agent model is evaluated by teacher fidelity, emergent-law robustness and mechanistic usefulness, not by beating the OT reference.",
        "",
        "## Tier Summary",
        "",
        mech.to_markdown(index=False),
        "",
        "## Retained Computational Hypotheses",
        "",
        retained[["law", "tier", "interpretation_level", "rollout_based", "directly_supervised_or_encoded"]].to_markdown(index=False) if not retained.empty else "None retained at acceptable evidence level.",
        "",
        *clone_lines,
        *atlas_lines,
        *grn_lines,
        "## Exploratory / Demonstration Only",
        "",
        exploratory[["law", "tier", "interpretation_level", "rollout_based", "directly_supervised_or_encoded"]].to_markdown(index=False) if not exploratory.empty else "None.",
        "",
        "## Unsupported Claims",
        "",
        unsupported[["law", "tier", "status"]].to_markdown(index=False) if not unsupported.empty else "None.",
        "",
        "## Core Metrics",
        "",
        mean.to_markdown(),
        "",
        "## Limitations",
        "",
        "- Current results are computational hypotheses.",
        "- Some laws are encoded control-law recoveries and must not be written as independent biological discoveries.",
        "- Native moscot teacher extraction removes the toy-fallback blocker for teacher construction, but not the need for external validation.",
        "- No experimental validation or causal mechanism is claimed.",
        "",
    ]
    out_path = "manuscript/final_retained_results_and_methods.quick_fixture.md" if quick_fixture else "manuscript/final_retained_results_and_methods.md"
    write_text(out_path, "\n".join(final_lines))
    write_text(
        report_dir / "discovery_hardening_summary.md",
        "\n".join(["# Discovery Hardening Summary", "", mech.to_markdown(index=False), "", laws.to_markdown(index=False), ""]),
    )
    write_text(
        report_dir / "editorial_assessment.md",
        "\n".join(
            [
                "# Editorial Assessment",
                "",
                "Current evidence level: discovery-hardened computational research prototype.",
                "",
                f"- teacher_fidelity_tier: {teacher_tier}",
                f"- emergent_law_tier: {emergent_tier}",
                f"- mechanistic_usefulness_tier: {mechanistic_tier}",
                "- Not ready for high-impact submission without external teacher validation and biological validation.",
                "" if quick_fixture else "- Clone-aware fate-diversification support is not a retained main claim unless the latest clone audit reaches acceptable cross-dataset support.",
                "" if quick_fixture or not clone_ctx else f"- latest clone audit: {clone_ctx['status']} ({clone_ctx['tier']}); {clone_ctx['interpretation']}",
                "" if quick_fixture or not atlas_ctx else f"- developmental atlas: {atlas_ctx['tier']}; {atlas_ctx['interpretation']}",
                "" if quick_fixture or not atlas_ctx else f"- developmental atlas analyzed datasets: {', '.join(atlas_ctx.get('analyzed_datasets', [])) or 'none'}",
                "" if quick_fixture or not atlas_ctx else "- E5 zebrafish is independent and native-moscot analyzed, but did not reproduce condensation-before-divergence and controls were not clean; it is boundary evidence, not validation.",
                "" if quick_fixture or not grn_ctx else f"- GRN/regulon audit: {grn_ctx['grn_tier']}; {grn_ctx['grn_allowed']}",
                "",
            ]
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/train.yaml")
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    cfg, model_cfg, discovery_cfg = configure_paths(args.config, args.quick_fixture)
    table_dir, report_dir, fig_dir = output_dirs(cfg, discovery_cfg)
    metrics = pd.read_csv(model_cfg.get("metrics_path", "tables/final_metrics.csv"))
    ensure_dir(table_dir)
    metrics.groupby("model")[["sinkhorn", "mmd_rbf", "energy", "knn_two_sample_accuracy", "celltype_composition_rmse"]].agg(["mean", "std"]).to_csv(table_dir / "final_model_metric_summary.csv")

    fidelity, teacher_tier = _teacher_fidelity(metrics, model_cfg, discovery_cfg, table_dir, report_dir, fig_dir)
    laws, emergent_tier = _run_discovery(args.config, args.quick_fixture, table_dir, report_dir)
    _backend, native_validated = _teacher_backend_status(model_cfg, table_dir, report_dir)
    mechanistic_tier = _mechanistic_tier(teacher_tier, laws, native_validated, discovery_cfg)
    _write_reports(metrics, fidelity, teacher_tier, laws, emergent_tier, mechanistic_tier, native_validated, cfg, table_dir, report_dir, args.quick_fixture)
    print(
        {
            "ot_reference": OT_REFERENCE,
            "primary_agent": PRIMARY_AGENT,
            "teacher_fidelity_tier": teacher_tier,
            "emergent_law_tier": emergent_tier,
            "mechanistic_usefulness_tier": mechanistic_tier,
        }
    )


if __name__ == "__main__":
    main()
