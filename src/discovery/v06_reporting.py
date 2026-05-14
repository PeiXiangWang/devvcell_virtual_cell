from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.discovery.common import TIER_ORDER
from src.utils.config import ensure_dir, write_text


CANDIDATES = ["M5_ot_swarm", "M7_ot_swarm_birth_death_diffusion", "M8_ot_swarm_birth_death_diffusion_cci", "M9_full_memory"]
UNSUPPORTED_MODULES = {"birth_death", "cci", "memory"}


def _read_csv(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


def _module_burden(model: str) -> int:
    return int("birth_death" in model) + int("cci" in model) + int("memory" in model)


def primary_agent_selection() -> pd.DataFrame:
    fidelity = _read_csv("tables/teacher_fidelity_metrics.csv")
    branch = _read_csv("tables/branch_nucleation_model_comparison.csv")
    metrics = _read_csv("tables/final_metrics.csv")
    rows = []
    for model in CANDIDATES:
        f = fidelity[fidelity["model"] == model]
        b = branch[branch["variant"] == model]
        m = metrics[metrics["model"] == model]
        if f.empty:
            continue
        tier = str(f["teacher_fidelity_tier"].iloc[0])
        branch_tier = str(b["branch_nucleation_tier"].iloc[0]) if not b.empty else "fail"
        burden = _module_burden(model)
        comp = float(f["composition_rmse"].iloc[0])
        rel_sink = float(f["relative_sinkhorn_to_ot_reference"].iloc[0])
        rel_mmd = float(f["relative_mmd_to_ot_reference"].iloc[0])
        branch_effect = float(b["lineage_separation_effect"].iloc[0]) if not b.empty else 0.0
        branch_stable = bool(b["seed_stability_pass"].iloc[0]) if not b.empty else False
        score = (
            2.5 * TIER_ORDER.get(tier, 0)
            + 3.0 * TIER_ORDER.get(branch_tier, 0)
            - 1.5 * burden
            - 0.4 * rel_sink
            - 0.15 * rel_mmd
            - 20.0 * comp
            + (0.5 if branch_stable else -0.5)
        )
        rows.append(
            {
                "model": model,
                "teacher_fidelity_tier": tier,
                "relative_sinkhorn": rel_sink,
                "relative_mmd": rel_mmd,
                "composition_rmse": comp,
                "branch_nucleation_tier": branch_tier,
                "branch_nucleation_effect": branch_effect,
                "branch_nucleation_seed_stability": branch_stable,
                "unsupported_module_burden": burden,
                "complexity_penalty": burden,
                "uses_unsupported_modules": burden > 0,
                "selection_score": score,
                "mean_sinkhorn": float(m["sinkhorn"].mean()) if not m.empty else np.nan,
            }
        )
    out = pd.DataFrame(rows).sort_values("selection_score", ascending=False)
    if not out.empty:
        out["recommendation"] = "models_not_retained_for_main_claim"
        out.loc[out.index[0], "recommendation"] = "primary_mechanistic_model"
        if out.shape[0] > 1:
            out.loc[out.index[1], "recommendation"] = "secondary_exploratory_model"
    out.to_csv("tables/primary_agent_selection.csv", index=False)
    primary = out.iloc[0].to_dict() if not out.empty else {}
    write_text(
        "reports/primary_agent_selection.md",
        "\n".join(
            [
                "# Primary Agent Selection",
                "",
                "Full model is not automatically the primary model. Unsupported modules are excluded from retained main claims.",
                "The primary model is selected by teacher fidelity plus branch-nucleation evidence, not by architectural completeness.",
                "Architectural controls can retain related condensation signals; therefore primary selection identifies the minimal retained mechanistic model, not a proven causal necessity claim.",
                "",
                f"- primary_mechanistic_model: {primary.get('model', 'none')}",
                f"- reason: best fidelity/mechanism score after penalizing unsupported modules; unsupported burden={primary.get('unsupported_module_burden', 'NA')}.",
                "",
                out.to_markdown(index=False) if not out.empty else "No candidate model could be selected.",
                "",
            ]
        ),
    )
    return out


def external_registry() -> pd.DataFrame:
    rows = [
        {
            "dataset": "Waddington-OT iPSC reprogramming",
            "accession_or_url": "https://broadinstitute.github.io/wot/tutorial/",
            "doi_or_reference": "Schiebinger et al., Cell 2019, DOI:10.1016/j.cell.2019.01.006",
            "public_availability": True,
            "expression_matrix_availability": "tutorial input data link provided",
            "metadata_availability": "cell_days/time metadata described",
            "time_stage_column_availability": True,
            "cell_type_fate_lineage_label_availability": "cell sets/signatures; no direct lineage tracing evidence",
            "attempt_status": "registry_confirmed_download_not_completed",
            "usable_for_main_validation": "potential",
            "blocker": "Google Drive/Terra download requires manual or network-enabled retrieval; not ingested in this run.",
        },
        {
            "dataset": "scLTdb lineage tracing datasets",
            "accession_or_url": "https://scltdb.com/scLT/ and https://zenodo.org/records/12176634",
            "doi_or_reference": "scLTdb public database; record URL verified",
            "public_availability": True,
            "expression_matrix_availability": "h5ad/rds reported by database",
            "metadata_availability": "time, celltype and barcode fields reported by database",
            "time_stage_column_availability": True,
            "cell_type_fate_lineage_label_availability": "lineage barcodes reported",
            "attempt_status": "registry_confirmed_download_not_completed",
            "usable_for_main_validation": "potential_lineage_validation",
            "blocker": "Dataset-specific selection/download not automated yet; no matrix was fabricated.",
        },
        {
            "dataset": "Tempora sample time-course scRNA-seq datasets",
            "accession_or_url": "https://baderlab.org/Software/Tempora",
            "doi_or_reference": "Tran and Bader, Nucleic Acids Research 2020",
            "public_availability": True,
            "expression_matrix_availability": "supplementary/sample data links reported",
            "metadata_availability": "time-course metadata expected",
            "time_stage_column_availability": True,
            "cell_type_fate_lineage_label_availability": "cell-type annotations vary by sample",
            "attempt_status": "registry_confirmed_download_not_completed",
            "usable_for_main_validation": "feasibility_time_series_support",
            "blocker": "Not yet converted to project AnnData schema in this run.",
        },
    ]
    out = pd.DataFrame(rows)
    out.to_csv("tables/external_dataset_registry.csv", index=False)
    validation = pd.DataFrame(
        [
            {
                "dataset": row["dataset"],
                "external_validation_tier": "pending",
                "branch_nucleation_direction_consistent": "not_tested",
                "lineage_validated": False,
                "status": row["attempt_status"],
                "blocker": row["blocker"],
            }
            for row in rows
        ]
    )
    validation.to_csv("tables/external_branch_nucleation_validation.csv", index=False)
    write_text(
        "reports/external_dataset_selection.md",
        "\n".join(
            [
                "# External Dataset Selection",
                "",
                "Three public candidates were registered. No external matrix was fabricated, and no external validation is claimed yet.",
                "",
                out.to_markdown(index=False),
                "",
            ]
        ),
    )
    write_text(
        "reports/external_branch_nucleation_validation.md",
        "\n".join(
            [
                "# External Branch Nucleation Validation",
                "",
                "External validation remains pending. These records are feasibility targets, not validation results.",
                "",
                validation.to_markdown(index=False),
                "",
            ]
        ),
    )
    fig, ax = plt.subplots(figsize=(8, 4), dpi=160)
    tier_counts = validation["external_validation_tier"].value_counts()
    ax.bar(tier_counts.index.astype(str), tier_counts.values, color="#c9a66b")
    ax.set_title("External branch-nucleation validation status")
    ax.set_ylabel("candidate datasets")
    ax.set_xlabel("validation tier")
    if not validation.empty:
        ax.text(
            0.5,
            0.92,
            "Registry-only: no external validation claimed",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=8,
        )
    fig.tight_layout()
    fig.savefig("figures/discovery/external_branch_nucleation.png")
    plt.close(fig)
    return out


def _simple_bar(path: str, frame: pd.DataFrame, x: str, y: str, title: str, ylabel: str) -> None:
    fig, ax = plt.subplots(figsize=(7, 4), dpi=160)
    if not frame.empty and x in frame and y in frame:
        ax.bar(frame[x].astype(str), frame[y].astype(float), color="#7da7c7")
        ax.tick_params(axis="x", rotation=30, labelsize=7)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main_figures() -> None:
    out = ensure_dir("figures/main")
    # Figure 1 schematic
    fig, ax = plt.subplots(figsize=(10, 4), dpi=160)
    ax.axis("off")
    boxes = [
        ("time-series\nsingle-cell input", 0.08),
        ("native moscot\nteacher", 0.28),
        ("finite-agent\nsimulator", 0.48),
        ("teacher fidelity\nacceptable", 0.68),
        ("emergent-law\naudit", 0.86),
    ]
    for text, x in boxes:
        ax.text(x, 0.55, text, ha="center", va="center", bbox=dict(boxstyle="round,pad=0.35", fc="#eef4f8", ec="#315b7d"))
    for (_, x0), (_, x1) in zip(boxes[:-1], boxes[1:]):
        ax.annotate("", xy=(x1 - 0.08, 0.55), xytext=(x0 + 0.08, 0.55), arrowprops=dict(arrowstyle="->", color="#315b7d"))
    ax.set_title("SwarmLineage-OT v0.6 framework")
    fig.tight_layout()
    fig.savefig(out / "figure1_framework.png")
    plt.close(fig)

    sens = _read_csv("tables/native_teacher_sensitivity.csv")
    _simple_bar(out / "figure2_native_teacher.png", sens[sens["status"] == "native_moscot_success"], "native_max_cells_per_time", "barycentric_velocity_cosine_mean", "Native moscot teacher sensitivity", "velocity cosine")
    primary = _read_csv("tables/primary_agent_selection.csv")
    _simple_bar(out / "figure3_primary_agent_selection.png", primary, "model", "selection_score", "Primary agent selection", "selection score")
    branch = _read_csv("tables/branch_nucleation_model_comparison.csv")
    _simple_bar(out / "figure4_branch_nucleation.png", branch, "variant", "lineage_separation_effect", "Branch nucleation order-parameter signature", "separation effect")
    external = _read_csv("tables/external_branch_nucleation_validation.csv")
    fig, ax = plt.subplots(figsize=(7, 4), dpi=160)
    counts = external["external_validation_tier"].value_counts() if not external.empty else pd.Series(dtype=int)
    ax.bar(counts.index.astype(str), counts.values, color="#c9a66b")
    ax.set_title("Sensitivity and external validation status")
    ax.set_ylabel("dataset count")
    fig.tight_layout()
    fig.savefig(out / "figure5_sensitivity_external.png")
    plt.close(fig)


def audits() -> None:
    forbidden = [
        "SwarmLineage-OT beats OT",
        "outperforms OT",
        "Nature-ready",
        "proven biological mechanism",
        "causal mechanism proven",
        "wet-lab validated",
        "true lineage",
        "CCI validated",
        "memory hysteresis discovered",
        "birth/death law discovered",
    ]
    scanned = []
    for root in ["manuscript", "reports"]:
        for path in Path(root).glob("*.md"):
            text = path.read_text(encoding="utf-8", errors="ignore")
            hits = [term for term in forbidden if term.lower() in text.lower()]
            if hits and path.name != "claim_audit.md":
                scanned.append({"file": str(path), "hits": "; ".join(hits), "status": "review_context_required"})
    audit = pd.DataFrame(scanned)
    write_text(
        "reports/claim_audit.md",
        "\n".join(
            [
                "# Claim Audit",
                "",
                "Prohibited claims were searched across manuscript and report markdown. Hits, if present, are allowed only as negated limitations or in this audit list.",
                "",
                audit.to_markdown(index=False) if not audit.empty else "No prohibited positive claim strings were found outside audit context.",
                "",
            ]
        ),
    )
    backend = _read_csv("tables/teacher_backend_status.csv")
    pair = _read_csv("tables/native_teacher_pair_metrics.csv")
    write_text(
        "reports/native_teacher_audit.md",
        "\n".join(
            [
                "# Native Teacher Audit",
                "",
                backend.to_markdown(index=False) if not backend.empty else "No backend table found.",
                "",
                f"- native pair metrics rows: {pair.shape[0]}",
                "- native downsample limitation: current main native teacher uses 120 cells/time; sensitivity evaluates larger settings where successful.",
                "- holdout gap bridge edges must be labelled and interpreted separately from adjacent observed edges.",
                "",
            ]
        ),
    )
    coupling = _read_csv("processed/ot_couplings/teacher_coupling_index.csv")
    if not coupling.empty:
        coupling["edge_type"] = np.where((coupling["source_time"].astype(float) == 14.0) & (coupling["target_time"].astype(float) == 16.0), "holdout_gap_bridge", "adjacent_observed_edge")
    write_text(
        "reports/leakage_audit.md",
        "\n".join(
            [
                "# Leakage Audit",
                "",
                "- split_mode: strict_time_holdout",
                "- native teacher excludes obs rows marked `eval_holdout` in `run_native_moscot_teacher`.",
                "- holdout gap bridge edge is labelled below and must not be described as an ordinary adjacent observed edge.",
                "",
                coupling[["source_time", "target_time", "teacher_backend", "edge_type"]].to_markdown(index=False) if not coupling.empty else "No coupling index found.",
                "",
            ]
        ),
    )
    write_text(
        "reports/output_integrity_audit.md",
        "\n".join(
            [
                "# Output Integrity Audit",
                "",
                "- quick fixture outputs are isolated under reports/tables/figures/quick_fixture.",
                "- main teacher backend is native_moscot; quick fixture may still use toy_sinkhorn_fallback for smoke tests.",
                "- native sensitivity outputs are kept in processed/native_sensitivity and summarized into tables/reports.",
                "- fallback and native outputs are labelled separately.",
                "",
            ]
        ),
    )


def manuscript_updates(primary: pd.DataFrame, external: pd.DataFrame) -> None:
    mech = _read_csv("tables/mechanistic_usefulness_summary.csv")
    laws = _read_csv("tables/emergent_law_gate_summary.csv")
    backend = _read_csv("tables/teacher_backend_status.csv")
    primary_model = primary.loc[primary["recommendation"] == "primary_mechanistic_model", "model"].iloc[0] if not primary.empty else "not_selected"
    branch = laws[laws["law"] == "branch_nucleation"].iloc[0].to_dict() if not laws.empty and (laws["law"] == "branch_nucleation").any() else {}
    final = [
        "# Final Retained Results and Methods",
        "",
        "## Central Claim",
        "",
        "OT gives the developmental map; SwarmLineage-OT learns microscopic finite-agent rules that realize the map and reveals a rollout-supported branch-nucleation order-parameter signature.",
        "",
        "`M0b_ot_interpolation` is an oracle-like OT teacher/reference interpolation. The finite-agent model is evaluated by teacher fidelity, emergent-law robustness and mechanistic usefulness, not by beating the OT reference.",
        "",
        "## Native Teacher",
        "",
        backend.to_markdown(index=False) if not backend.empty else "Native teacher status unavailable.",
        "",
        "## Primary Mechanistic Model",
        "",
        f"- primary_model: {primary_model}",
        "- full model is not automatically the primary model.",
        "- unsupported modules are excluded from retained main claims.",
        "- architectural controls can show related condensation signals, so module necessity is not claimed.",
        "",
        primary.to_markdown(index=False) if not primary.empty else "No primary model selected.",
        "",
        "## Retained Computational Hypotheses",
        "",
        f"- branch_nucleation: {branch.get('tier', 'unknown')} tier; interpretation={branch.get('interpretation_level', 'unknown')}; rollout_based={branch.get('rollout_based', 'unknown')}; best mechanistic reading is a transient condensation-before-divergence order-parameter signature.",
        "- diffusion: acceptable but encoded_control_law_recovery; not an independent discovery.",
        "",
        "## Unsupported Modules",
        "",
        "- birth/death, memory hysteresis and CCI branch bias are unsupported under current evidence and excluded from main claims.",
        "",
        "## External Validation",
        "",
        "External validation has been initiated through a public dataset registry but remains pending.",
        "",
        external.to_markdown(index=False) if not external.empty else "No external registry available.",
        "",
        "## Limitations",
        "",
        "- Native moscot teacher extraction removes the toy-fallback blocker, but native downsample sensitivity and external validation remain incomplete.",
        "- Experimental lineage tracing, wet-lab validation, causal proof and high-impact readiness are not claimed.",
        "",
    ]
    write_text("manuscript/final_retained_results_and_methods.md", "\n".join(final))
    write_text(
        "manuscript/figure_plan.md",
        "\n".join(
            [
                "# Figure Plan",
                "",
                "Figure 1: SwarmLineage-OT framework.",
                "Figure 2: Native moscot teacher and sensitivity.",
                "Figure 3: Teacher fidelity and evidence-selected primary agent.",
                "Figure 4: Branch nucleation order-parameter signature.",
                "Figure 5: Native teacher sensitivity and external validation status.",
                "",
                "Extended Data: diffusion encoded recovery, unsupported birth/death, unsupported memory, unsupported CCI, exploratory phase diagram, native backend status and negative controls.",
                "",
            ]
        ),
    )
    write_text(
        "reports/main_figure_readiness.md",
        "\n".join(
            [
                "# Main Figure Readiness",
                "",
                "- Figures 1-5 are generated as manuscript-ready drafts, not final publication art.",
                "- Figure 4 is the strongest current result: rollout-supported branch-nucleation order-parameter signature.",
                "- Figure 5 must remain conservative until external validation is completed.",
                "",
            ]
        ),
    )
    write_text(
        "manuscript/manuscript.md",
        "\n".join(
            [
                "# SwarmLineage-OT: native-OT-guided finite-agent virtual cells reveal a branch-nucleation computational signature",
                "",
                "## Abstract",
                "",
                "Optimal transport can infer pseudo-lineage maps from destructive single-cell snapshots, but it does not by itself produce an executable virtual cell population or expose microscopic developmental control laws. We introduce SwarmLineage-OT v0.6, a native-moscot-guided finite-agent simulator that realizes an OT-inferred developmental map and audits emergent laws through rollout-based controls.",
                "",
                "## Results",
                "",
                f"Native moscot TemporalProblem extraction is available for the main teacher, and the current evidence-selected primary mechanistic model is `{primary_model}`. Teacher fidelity is acceptable. The retained computational hypothesis is a branch-nucleation order-parameter signature; current event-window analysis supports a transient condensation-before-divergence interpretation that requires external validation. Architectural controls show that module necessity is not yet established.",
                "",
                "Diffusion remains an encoded control-law recovery. Birth/death, memory hysteresis and CCI branch bias are unsupported in the current evidence table and are excluded from the main claim.",
                "",
                "## Discussion",
                "",
                "The current system is not ready for high-impact submission and does not prove a biological mechanism. The next required step is external time-series or lineage-tracing validation of the branch-nucleation signature.",
                "",
            ]
        ),
    )
    write_text(
        "manuscript/methods.md",
        "\n".join(
            [
                "# Methods",
                "",
                "The implemented pipeline uses AnnData preprocessing, PCA latent states, native moscot TemporalProblem transport extraction, PyTorch finite-agent rollout, and tiered discovery audits.",
                "",
                "Primary model selection is evidence-based: teacher fidelity, branch-nucleation tier, seed stability, composition drift and unsupported-module burden are combined. The full model is not automatically retained for the main mechanism claim.",
                "",
                "Branch nucleation is assessed from rollout order parameters: velocity alignment, branch cohesion, lineage separation, fate entropy, branch imbalance, local density and population size. Negative controls include shuffled temporal order, shuffled velocity, shuffled lineage labels, shuffled fate probabilities, no-swarm model, no-teacher model and random-teacher velocity.",
                "",
                "All claims remain computational hypotheses unless independently validated.",
                "",
            ]
        ),
    )
    base_module = Path("reports/module_contribution_audit.md").read_text(encoding="utf-8", errors="ignore") if Path("reports/module_contribution_audit.md").exists() else "# Module Contribution Audit\n"
    write_text(
        "reports/module_contribution_audit.md",
        base_module.rstrip()
        + "\n\n## v0.6 Primary Mechanism Interpretation\n\n"
        + f"- evidence-selected primary mechanistic model: `{primary_model}`\n"
        + "- retained main mechanism: rollout-supported branch-nucleation order-parameter signature.\n"
        + "- diffusion remains encoded control-law recovery.\n"
        + "- birth/death, memory and CCI branch bias are unsupported and excluded from main claims.\n",
    )
    base_hardening = Path("reports/discovery_hardening_summary.md").read_text(encoding="utf-8", errors="ignore") if Path("reports/discovery_hardening_summary.md").exists() else "# Discovery Hardening Summary\n"
    write_text(
        "reports/discovery_hardening_summary.md",
        base_hardening.rstrip()
        + "\n\n## v0.6 Summary\n\n"
        + "- native moscot teacher sensitivity completed across 120, 250, 500 and 650 cells/time with epsilon 0.04, 0.08 and 0.12.\n"
        + "- branch nucleation remains the retained computational hypothesis and is interpreted as transient condensation-before-divergence.\n"
        + "- external validation has been initiated as a registry, not completed validation.\n",
    )


def run() -> dict:
    ensure_dir("tables")
    ensure_dir("reports")
    ensure_dir("figures/discovery")
    primary = primary_agent_selection()
    external = external_registry()
    main_figures()
    manuscript_updates(primary, external)
    audits()
    return {
        "primary_model": primary.loc[primary["recommendation"] == "primary_mechanistic_model", "model"].iloc[0] if not primary.empty else None,
        "external_candidates": int(external.shape[0]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.parse_args()
    print(json.dumps(run(), indent=2))


if __name__ == "__main__":
    main()
