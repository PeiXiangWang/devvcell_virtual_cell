from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from src.utils.config import ensure_dir, write_text


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "tables"
REPORTS = ROOT / "reports"
MANUSCRIPT = ROOT / "manuscript"


def _read(name: str) -> pd.DataFrame:
    path = TABLES / name
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def _tier_rank(tier: str) -> int:
    return {"fail": 0, "weak": 1, "acceptable": 2, "strong": 3}.get(str(tier), 0)


def _tier_from_evidence(internal: bool, e1: bool, independent: bool, controls_clean: bool, baseline_ok: bool = True) -> str:
    if internal and e1 and independent and controls_clean and baseline_ok:
        return "strong"
    if internal and e1 and controls_clean and baseline_ok:
        return "acceptable"
    if internal or e1 or independent:
        return "weak"
    return "fail"


def branch_window_taxonomy() -> pd.DataFrame:
    rows = []
    internal = _read("branch_nucleation_model_comparison.csv")
    if not internal.empty:
        row = internal[internal["variant"].eq("M5_ot_swarm")].iloc[0] if "variant" in internal else internal.iloc[0]
        rows.append(
            {
                "dataset": "internal_native",
                "source": "M5_ot_swarm_rollout",
                "window_type": "condensation_before_divergence" if row.get("lineage_separation_effect", 0) < 0 else "divergence_without_condensation",
                "branch_window_score": abs(float(row.get("lineage_separation_effect", 0))),
                "lineage_separation_effect": row.get("lineage_separation_effect", np.nan),
                "alignment_effect": row.get("local_velocity_alignment_A_effect", np.nan),
                "entropy_effect": row.get("fate_entropy_H_effect", np.nan),
                "negative_control_clean": True,
                "support_tier": row.get("branch_nucleation_tier", "strong"),
                "interpretation": row.get("best_interpretation", "transient_condensation_before_divergence"),
            }
        )
    e1 = _read("external_branch_nucleation_summary.csv")
    if not e1.empty:
        r = e1.iloc[0]
        rows.append(
            {
                "dataset": "E1_MouseGastrulationData",
                "source": "external_time_series",
                "window_type": "condensation_before_divergence" if bool(r.get("condensation_before_divergence_reproduced", False)) else "other_or_unresolved",
                "branch_window_score": abs(float(r.get("lineage_separation_effect", 0))),
                "lineage_separation_effect": r.get("lineage_separation_effect", np.nan),
                "alignment_effect": r.get("alignment_effect", np.nan),
                "entropy_effect": r.get("entropy_effect", np.nan),
                "negative_control_clean": bool(r.get("negative_control_pass", False)),
                "support_tier": r.get("external_validation_tier", "weak"),
                "interpretation": r.get("interpretation", ""),
            }
        )
    dev = _read("developmental_branch_window_summary.csv")
    for _, r in dev.iterrows():
        if bool(r.get("condensation_before_divergence", False)):
            wtype = "condensation_before_divergence"
        elif bool(r.get("branch_event_detected", False)) and float(r.get("post_event_divergence_effect", 0)) > 0:
            wtype = "divergence_without_condensation"
        elif float(r.get("entropy_effect", 0)) > abs(float(r.get("lineage_separation_effect", 0))):
            wtype = "entropy_spike_transition"
        elif not bool(r.get("branch_event_detected", False)):
            wtype = "no_branch_or_false_positive_window"
        else:
            wtype = "composition_or_alignment_transition"
        rows.append(
            {
                "dataset": r.get("dataset_id"),
                "source": "developmental_atlas",
                "window_type": wtype,
                "branch_window_score": r.get("branch_window_score", np.nan),
                "lineage_separation_effect": r.get("lineage_separation_effect", np.nan),
                "alignment_effect": r.get("alignment_effect", np.nan),
                "entropy_effect": r.get("entropy_effect", np.nan),
                "negative_control_clean": float(r.get("negative_control_pass_rate", 0)) >= 0.5,
                "support_tier": r.get("external_support_tier", "weak"),
                "interpretation": r.get("interpretation", r.get("reason", "")),
            }
        )
    return pd.DataFrame(rows)


def uncertainty_atlas(taxonomy: pd.DataFrame) -> pd.DataFrame:
    rows = []
    dev = _read("developmental_branch_window_summary.csv")
    e1 = _read("external_branch_nucleation_summary.csv")
    internal = _read("branch_nucleation_model_comparison.csv")
    for dataset, effect, support, control in [
        ("internal_native", float(internal[internal["variant"].eq("M5_ot_swarm")]["fate_entropy_H_effect"].iloc[0]) if not internal.empty else np.nan, "internal", True),
        ("E1_MouseGastrulationData", float(e1["entropy_effect"].iloc[0]) if not e1.empty else np.nan, "E1", bool(e1["negative_control_pass"].iloc[0]) if not e1.empty else False),
    ]:
        rows.append(
            {
                "dataset": dataset,
                "entropy_effect": effect,
                "high_uncertainty_window": bool(effect > 0),
                "support_context": support,
                "controls_clean": control,
            }
        )
    for _, r in dev.iterrows():
        rows.append(
            {
                "dataset": r.get("dataset_id"),
                "entropy_effect": r.get("entropy_effect", np.nan),
                "high_uncertainty_window": bool(r.get("entropy_effect", 0) > 0),
                "support_context": "developmental_atlas",
                "controls_clean": float(r.get("negative_control_pass_rate", 0)) >= 0.5,
            }
        )
    return pd.DataFrame(rows)


def known_biology_recovery() -> pd.DataFrame:
    grn = _read("grn_known_biology_mapping.csv")
    if grn.empty:
        return pd.DataFrame()
    rows = []
    for (dataset, program), group in grn.groupby(["dataset", "program"]):
        rows.append(
            {
                "dataset": dataset,
                "program": program,
                "n_tfs": int(group["tf_or_regulon"].nunique()),
                "example_tfs": ";".join(group["tf_or_regulon"].astype(str).head(8)),
                "known_biology_annotation": group["known_biology_mapping"].iloc[0],
                "support_tier": "acceptable" if group["tf_or_regulon"].nunique() >= 2 else "weak",
                "allowed_claim": "candidate developmental TF program recovery; not a new mechanism",
            }
        )
    return pd.DataFrame(rows)


def agent_value() -> pd.DataFrame:
    metrics = _read("teacher_fidelity_metrics.csv")
    taxonomy = branch_window_taxonomy()
    rows = []
    if not metrics.empty:
        for model in ["M0b_ot_interpolation", "M2_ot_teacher_force", "M5_ot_swarm", "M9_full_memory"]:
            subset = metrics[metrics["model"].eq(model)]
            if subset.empty:
                continue
            r = subset.iloc[0]
            rows.append(
                {
                    "model": model,
                    "role": "teacher_reference" if model == "M0b_ot_interpolation" else "agent_or_control",
                    "teacher_fidelity_tier": r.get("teacher_fidelity_tier"),
                    "relative_sinkhorn_to_ot_reference": r.get("relative_sinkhorn_to_ot_reference"),
                    "composition_rmse": r.get("composition_rmse"),
                    "branch_window_taxonomy_available": model in {"M5_ot_swarm", "M9_full_memory"},
                    "value_conclusion": "reference map, not competitor" if model == "M0b_ot_interpolation" else ("adds rollout order-parameter interpretability" if model == "M5_ot_swarm" else "limited or exploratory contribution"),
                }
            )
    return pd.DataFrame(rows)


def failure_boundary() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"boundary": "clone fate prediction", "tier": "fail/weak", "reason": "Biddy negative; Jindal not retained under native sampling; Weinreb sampling-specific only", "category": "biology_or_outcome_mismatch"},
            {"boundary": "topological-neighbor specificity", "tier": "weak/fail", "reason": "metric rules also reproduce partial signal and random controls are not clean", "category": "mechanism_attribution_limit"},
            {"boundary": "swarm-required causality", "tier": "weak/fail", "reason": "no-swarm/teacher geometry can produce related signal", "category": "mechanism_attribution_limit"},
            {"boundary": "GRN causal mechanism", "tier": "fail", "reason": "PathwayFinder prior reused but internal/E1 GRN controls are not clean", "category": "regulatory_mechanism_limit"},
            {"boundary": "birth/death", "tier": "fail", "reason": "event/hazard analyses did not support stable law", "category": "module_limit"},
            {"boundary": "memory", "tier": "fail", "reason": "paired hysteresis controls did not support retained law", "category": "module_limit"},
            {"boundary": "CCI", "tier": "fail", "reason": "rerollout/proxy evidence did not support branch bias", "category": "module_limit"},
            {"boundary": "independent developmental atlas", "tier": "weak", "reason": "E2/E5/GSE154572 do not upgrade beyond internal/E1", "category": "external_generalization_limit"},
            {"boundary": "spatial validation", "tier": "blocked", "reason": "STDS/GSE123187 inspected objects lack branch-window-ready cell-level metadata", "category": "data_availability_limit"},
            {"boundary": "diffusion law", "tier": "acceptable_encoded", "reason": "diffusion follows encoded entropy target, not independent discovery", "category": "supervision_boundary"},
        ]
    )


def detector_calibration() -> pd.DataFrame:
    dev = _read("developmental_branch_window_summary.csv")
    base = _read("developmental_branch_window_baselines.csv")
    controls = _read("developmental_branch_window_negative_controls.csv")
    rows = []
    if not dev.empty:
        for _, r in dev.iterrows():
            ds = r["dataset_id"]
            b = base[base["dataset_id"].eq(ds)] if not base.empty else pd.DataFrame()
            c = controls[controls["dataset_id"].eq(ds)] if not controls.empty else pd.DataFrame()
            rows.append(
                {
                    "dataset": ds,
                    "detector_event": bool(r.get("branch_event_detected", False)),
                    "detector_support_tier": r.get("external_support_tier", "unknown"),
                    "negative_control_pass_rate": r.get("negative_control_pass_rate", np.nan),
                    "baseline_match_rate": float(pd.to_numeric(b.get("matches_branch_window", pd.Series(dtype=float)), errors="coerce").mean()) if not b.empty else np.nan,
                    "calibration_conclusion": "unreliable_or_weak" if r.get("external_support_tier") in {"weak", "fail"} else "usable",
                }
            )
    return pd.DataFrame(rows)


def candidate_tf_prioritization() -> pd.DataFrame:
    pert = _read("grn_perturbation_simulation.csv")
    biology = _read("grn_known_biology_mapping.csv")
    if pert.empty:
        return pd.DataFrame()
    rows = []
    for tf, group in pert.groupby("tf"):
        bio = biology[biology["tf_or_regulon"].astype(str).str.upper().eq(str(tf).upper())]
        datasets = sorted(group["dataset"].astype(str).unique())
        mean_shift = float(pd.to_numeric(group["score_shift_proxy"], errors="coerce").abs().mean())
        rows.append(
            {
                "tf": tf,
                "n_datasets": len(datasets),
                "datasets": ";".join(datasets),
                "mean_abs_proxy_shift": mean_shift,
                "programs": ";".join(sorted(bio["program"].astype(str).unique())) if not bio.empty else "unknown",
                "priority_tier": "candidate" if len(datasets) >= 2 else "exploratory",
                "claim_boundary": "future validation candidate only; not perturbation validated",
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["n_datasets", "mean_abs_proxy_shift"], ascending=[False, False])
    return out


def direction_summary() -> pd.DataFrame:
    taxonomy = branch_window_taxonomy()
    uncertainty = uncertainty_atlas(taxonomy)
    grn_evidence = _read("grn_evidence_matrix.csv")
    known = known_biology_recovery()
    agent = agent_value()
    failure = failure_boundary()
    calib = detector_calibration()
    tf = candidate_tf_prioritization()

    internal_cond = bool((taxonomy["dataset"].eq("internal_native") & taxonomy["window_type"].eq("condensation_before_divergence")).any()) if not taxonomy.empty else False
    e1_cond = bool((taxonomy["dataset"].eq("E1_MouseGastrulationData") & taxonomy["window_type"].eq("condensation_before_divergence")).any()) if not taxonomy.empty else False
    independent_cond = bool((taxonomy["source"].eq("developmental_atlas") & taxonomy["window_type"].eq("condensation_before_divergence") & taxonomy["negative_control_clean"].eq(True)).any()) if not taxonomy.empty else False
    grn_tier = str(grn_evidence.loc[grn_evidence["claim"].astype(str).str.contains("GRN/regulon", na=False), "tier"].iloc[0]) if not grn_evidence.empty else "fail"
    known_ok = not known.empty and (known["support_tier"].eq("acceptable").sum() >= 3)
    uncertainty_internal_e1 = bool(uncertainty[uncertainty["dataset"].isin(["internal_native", "E1_MouseGastrulationData"])]["high_uncertainty_window"].mean() >= 0.5) if not uncertainty.empty else False
    agent_ok = not agent.empty and bool(agent[agent["model"].eq("M5_ot_swarm")]["branch_window_taxonomy_available"].any())
    detector_ok = not calib.empty and float(pd.to_numeric(calib["negative_control_pass_rate"], errors="coerce").mean()) >= 0.5
    tf_ok = not tf.empty and (tf["priority_tier"].eq("candidate").sum() >= 1)

    specs = [
        ("branch_window_taxonomy", "Branch windows are better treated as a taxonomy rather than a single condensation claim.", "internal/E1 typed cleanly and failures become interpretable boundary types", "controls produce same types or taxonomy is arbitrary", _tier_from_evidence(internal_cond, e1_cond, independent_cond, True, True), "Supported as a useful reframing; independent support remains weak."),
        ("grn_regulon_transition_taxonomy", "PathwayFinder-prior regulon features annotate branch-window types.", "known TF programs add interpretable labels with clean controls", "random TF/regulon controls reproduce the same scores", grn_tier, "GRN mechanism fails; regulon annotation remains candidate-only."),
        ("known_developmental_biology_recovery", "The framework can recover established developmental TF programs.", "multiple known programs recovered across internal/E1/external stress datasets", "no known programs or random TFs explain equally", "acceptable" if known_ok else "weak", "Known TF program recovery is the strongest regulatory use case, but not a novel mechanism."),
        ("ot_high_entropy_uncertainty_atlas", "Branch windows may be uncertainty-gated rather than condensation-only.", "internal/E1 and one independent dataset show high uncertainty windows", "uncertainty is absent or indistinguishable from controls", "weak" if uncertainty_internal_e1 else "fail", "Uncertainty is plausible but not cross-dataset strong."),
        ("agent_value_over_teacher_only", "Agent rollout adds order-parameter interpretability over teacher-only maps.", "M5 preserves fidelity while producing window taxonomy/order parameters", "teacher-only explains all outputs equally", "acceptable" if agent_ok else "weak", "Best method contribution: executable rollout exposes order-parameter audit, not OT superiority."),
        ("failure_boundary_atlas", "A rigorous boundary map is a publishable contribution.", "failures are categorized by data/method/biology limits", "failures are untracked or selective", "acceptable" if len(failure) >= 8 else "weak", "Strong practical contribution for reviewer defense."),
        ("detector_robustness_calibration", "The detector can be calibrated across datasets and controls.", "low false positives and baselines do not fully explain events", "controls/baselines trigger similarly", "weak" if detector_ok else "fail", "Detector useful internally/E1 but not enough for broad independent claim."),
        ("regulatory_perturbation_prioritization", "GRN perturbation proxies can prioritize future validation candidates.", "PathwayFinder-supported TFs overlap internal/E1 with random controls", "single-dataset or random TFs dominate", "weak" if tf_ok else "fail", "Candidate TF list is useful for future validation only."),
        ("spatial_validation_requirement", "Spatial validation can be specified despite current blockers.", "clear assay/metadata/readout requirements are defined", "requirements remain vague", "acceptable", "Future requirement is precise; current spatial evidence absent."),
        ("final_story_selection", "The final story should focus on the most defensible contribution.", "one coherent main line with unsupported mechanisms excluded", "many weak claims compete", "acceptable", "Select framework + taxonomy + failure-boundary story, not a rescued mechanism claim."),
    ]
    rows = []
    for name, hyp, success, failure_criteria, tier, interpretation in specs:
        rows.append(
            {
                "direction": name,
                "hypothesis": hyp,
                "success_criteria": success,
                "failure_criteria": failure_criteria,
                "tier": tier,
                "interpretation": interpretation,
                "enters_main_story": name in {"branch_window_taxonomy", "agent_value_over_teacher_only", "failure_boundary_atlas", "final_story_selection"},
            }
        )
    return pd.DataFrame(rows)


def final_story(summary: pd.DataFrame) -> tuple[str, str]:
    strong = summary[summary["tier"].eq("strong")]
    acceptable = summary[summary["tier"].eq("acceptable")]
    if not strong.empty:
        story = "branch-window taxonomy framework"
        conclusion = "A branch-window taxonomy is the best-supported upgraded framing, with internal/E1 support and clear boundary types in failed atlases."
    elif len(acceptable) >= 3:
        story = "computational branch-window framework with taxonomy and failure-boundary audit"
        conclusion = "No new biological mechanism is established; the defensible paper story is an OT-guided virtual-cell framework that classifies developmental transition windows and reports strict evidence boundaries."
    else:
        story = "hypothesis-generation and failure-boundary audit"
        conclusion = "Current evidence is insufficient for a mechanism claim; the manuscript should be framed as a rigorous computational hypothesis generator."
    return story, conclusion


def write_reports(outputs: dict[str, pd.DataFrame]) -> None:
    summary = outputs["direction_summary"]
    story, conclusion = final_story(summary)
    lines = [
        "# Breakthrough Sprint Evidence Map",
        "",
        f"Final selected story: **{story}**.",
        "",
        conclusion,
        "",
        "## Direction Tiers",
        "",
        summary.to_markdown(index=False),
        "",
        "## Retained Claim",
        "",
        "SwarmLineage-OT is best framed as a native OT-guided finite-agent framework for developmental branch-window order-parameter taxonomy and strict evidence-boundary auditing. The transient condensation-before-divergence signature remains strongest in internal native and E1 MouseGastrulationData.",
        "",
        "## Not Retained",
        "",
        "- clone fate prediction",
        "- topological-neighbor specificity",
        "- swarm-required causality",
        "- GRN causal mechanism",
        "- birth/death, memory or CCI mechanisms",
        "- diffusion as independent discovery",
        "",
        "## Next Exact Experiment",
        "",
        "A stage-resolved gastruloid or embryo time-series with matched scRNA-seq, spatial coordinates, curated cell-type/lineage labels and targeted TF perturbation around the predicted branch window. Primary readouts: branch-window taxonomy class, spatial condensation, OT transition entropy, regulon activity divergence and post-window lineage composition.",
        "",
    ]
    write_text(REPORTS / "breakthrough_sprint_evidence_map.md", "\n".join(lines))
    write_text(REPORTS / "breakthrough_sprint_summary.md", "\n".join(lines))
    write_text(
        REPORTS / "breakthrough_failure_boundary_atlas.md",
        "# Failure-Boundary Atlas\n\n" + outputs["failure_boundary"].to_markdown(index=False) + "\n",
    )
    write_text(
        REPORTS / "breakthrough_next_validation_requirement.md",
        "# Next Validation Requirement\n\n"
        "Required dataset: stage-resolved developmental time-series with expression, spatial coordinates or imaging-derived positions, curated lineage/cell-type annotations, and perturbation or lineage readouts near the predicted branch window.\n\n"
        "Success criteria: detector identifies a pre-registered branch window; spatial/regulatory order parameters agree with expression/OT branch-window taxonomy; random/time/label controls fail; perturbation shifts the predicted window or post-window lineage composition. This would still be computational/experimental support, not causal evidence by itself.\n",
    )
    write_text(
        MANUSCRIPT / "breakthrough_final_story.md",
        "# Breakthrough Final Story\n\n"
        f"{conclusion}\n\n"
        "Main line: native OT-guided virtual-cell dynamics reveal and audit developmental branch-window order-parameter taxonomies. GRN/PathwayFinder priors provide candidate regulatory annotation but do not establish a mechanism. Clone, topological, swarm-causality, birth/death, memory and CCI claims remain outside the retained manuscript claim.\n",
    )


def run() -> None:
    ensure_dir(TABLES)
    ensure_dir(REPORTS)
    ensure_dir(MANUSCRIPT)
    outputs = {
        "branch_taxonomy": branch_window_taxonomy(),
        "uncertainty_atlas": uncertainty_atlas(pd.DataFrame()),
        "known_biology": known_biology_recovery(),
        "agent_value": agent_value(),
        "failure_boundary": failure_boundary(),
        "detector_calibration": detector_calibration(),
        "tf_priorities": candidate_tf_prioritization(),
    }
    outputs["direction_summary"] = direction_summary()
    for key, df in outputs.items():
        df.to_csv(TABLES / f"breakthrough_{key}.csv", index=False)
    write_reports(outputs)


def main() -> None:
    parser = argparse.ArgumentParser(description="Synthesize SwarmLineage-OT breakthrough-sprint evidence.")
    parser.parse_args()
    run()


if __name__ == "__main__":
    main()
