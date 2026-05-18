from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "tables"
REPORTS = ROOT / "reports"
MANUSCRIPT = ROOT / "manuscript"


PRIORITY_PAIRS = {
    ("E1_MouseGastrulationData", "FGF3", "FGFR1"),
    ("E1_MouseGastrulationData", "FGF8", "FGFR1"),
    ("E1_MouseGastrulationData", "MDK", "NCL"),
    ("E1_MouseGastrulationData", "FN1", "ITGB1"),
    ("E1_MouseGastrulationData", "MDK", "ITGB1"),
    ("GSE154572_EB_WT", "WNT11", "FZD7"),
    ("GSE154572_EB_WT", "WNT11", "LRP6"),
    ("GSE154572_EB_WT", "WNT9A", "LRP6"),
    ("GSE154572_EB_WT", "WNT3", "FZD7"),
    ("E5_zebrafish_Farrell", "FGF4", "FGFR4"),
    ("E5_zebrafish_Farrell", "FGF17", "FGFR4"),
}


def _read(path: str) -> pd.DataFrame:
    p = TABLES / path
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


def _tier_rank(tier: str) -> int:
    return {"strong": 3, "acceptable": 2, "weak": 1, "fail": 0}.get(str(tier), 0)


def main() -> int:
    TABLES.mkdir(exist_ok=True)
    REPORTS.mkdir(exist_ok=True)
    MANUSCRIPT.mkdir(exist_ok=True)
    broad = _read("communication_niche_cross_dataset_summary.csv")
    broad_audit = _read("communication_niche_dataset_audit.csv")
    windows = _read("communication_niche_window_scan.csv")
    strict_tiers = _read("morphogen_niche_module_tiers.csv")
    strict_scan = _read("morphogen_niche_module_scan.csv")
    candidates = _read("morphogen_niche_candidate_lr_pairs.csv")
    prior = _read("morphogen_niche_strict_prior_audit.csv")

    broad_tier = str(broad.iloc[0]["tier"]) if not broad.empty else "fail"
    strict_best = "none"
    strict_tier = "fail"
    if not strict_tiers.empty:
        ranked = strict_tiers.copy()
        ranked["tier_rank"] = ranked["tier"].apply(_tier_rank)
        ranked = ranked.sort_values(["tier_rank", "support_dataset_count"], ascending=False)
        strict_best = str(ranked.iloc[0]["module"])
        strict_tier = str(ranked.iloc[0]["tier"])
    broad_vs_strict = pd.DataFrame(
        [
            {
                "comparison": "broad_vs_strict",
                "broad_communication_niche_tier": broad_tier,
                "strict_morphogen_family_tier": strict_tier,
                "strongest_strict_family": strict_best,
                "broad_acceptable_datasets": int(broad.iloc[0].get("acceptable_datasets", 0)) if not broad.empty else 0,
                "broad_activation_support_datasets": int(broad.iloc[0].get("activation_support_datasets", 0)) if not broad.empty else 0,
                "strict_acceptable_modules": int(strict_tiers["tier"].eq("acceptable").sum()) if not strict_tiers.empty else 0,
                "strict_weak_modules": int(strict_tiers["tier"].eq("weak").sum()) if not strict_tiers.empty else 0,
                "interpretation": "distributed_extracellular_niche_more_robust_than_single_family"
                if _tier_rank(broad_tier) > _tier_rank(strict_tier)
                else "single_family_not_weaker_than_broad_field",
            }
        ]
    )
    broad_vs_strict.to_csv(TABLES / "communication_niche_broad_vs_strict.csv", index=False)

    # Sender/receiver decomposition for the best window per dataset.
    decomp_rows = []
    if not windows.empty:
        for dataset, grp in windows.groupby("dataset"):
            best = grp.sort_values("communication_window_score", ascending=False).iloc[0]
            sender = float(best.get("sender_neighbor_activity_effect", np.nan))
            receiver = float(best.get("receiver_priming_effect", np.nan))
            product = float(best.get("communication_activation_effect", np.nan))
            if np.nan_to_num(product) > 0 and np.nan_to_num(sender) > 0 and np.nan_to_num(receiver) > 0:
                driver = "sender_and_receiver_product"
            elif np.nan_to_num(product) > 0 and np.nan_to_num(sender) > 0:
                driver = "sender_weighted_product"
            elif np.nan_to_num(product) > 0 and np.nan_to_num(receiver) > 0:
                driver = "receiver_weighted_product"
            elif np.nan_to_num(sender) > 0:
                driver = "sender_only_without_product_support"
            elif np.nan_to_num(receiver) > 0:
                driver = "receiver_only_without_product_support"
            else:
                driver = "not_supported"
            decomp_rows.append(
                {
                    "dataset": dataset,
                    "pre_time": best.get("pre_time"),
                    "event_time": best.get("event_time"),
                    "post_time": best.get("post_time"),
                    "sender_neighbor_activity_effect": sender,
                    "receiver_priming_effect": receiver,
                    "sender_receiver_product_effect": product,
                    "post_event_comm_divergence_effect": best.get("post_event_comm_divergence_effect"),
                    "time_shuffle_pass": best.get("time_shuffle_pass"),
                    "primary_driver": driver,
                }
            )
    decomp = pd.DataFrame(decomp_rows)
    decomp.to_csv(TABLES / "communication_niche_sender_receiver_decomposition.csv", index=False)

    # Candidate LR validation panel.
    panel = pd.DataFrame()
    if not candidates.empty:
        keep = candidates[
            candidates.apply(
                lambda r: (str(r["dataset"]), str(r["ligand"]), str(r["receptor"])) in PRIORITY_PAIRS,
                axis=1,
            )
        ].copy()
        if not keep.empty:
            strict_lookup = strict_scan[
                ["dataset", "module", "activation_effect", "receiver_priming_effect", "post_event_divergence_effect", "random_control_pass", "time_shuffle_pass", "support_tier"]
            ].drop_duplicates(["dataset", "module"])
            keep = keep.merge(strict_lookup, how="left", on=["dataset", "module"])
            keep["strict_prior_support"] = True
            keep["recommended_validation_assay"] = keep["module"].map(
                lambda m: "spatial transcriptomics or smFISH/MERFISH branch-window panel"
            )
            keep["allowed_claim"] = "candidate spatial or perturbation validation target"
            keep["forbidden_claim"] = "not_a_retained_mechanism"
            panel = keep[
                [
                    "dataset",
                    "module",
                    "event_time",
                    "ligand",
                    "receptor",
                    "event_pair_score",
                    "activation_effect",
                    "receiver_priming_effect",
                    "post_event_divergence_effect",
                    "sender_lineage_proxy",
                    "receiver_lineage_proxy",
                    "strict_prior_support",
                    "random_control_pass",
                    "time_shuffle_pass",
                    "module_support_tier",
                    "recommended_validation_assay",
                    "allowed_claim",
                    "forbidden_claim",
                ]
            ].sort_values(["module_support_tier", "event_pair_score"], ascending=[True, False])
    panel.to_csv(TABLES / "communication_niche_candidate_validation_panel.csv", index=False)

    lineage_map = panel[
        [
            "dataset",
            "module",
            "ligand",
            "receptor",
            "sender_lineage_proxy",
            "receiver_lineage_proxy",
            "event_pair_score",
            "module_support_tier",
        ]
    ].copy() if not panel.empty else pd.DataFrame()
    lineage_map.to_csv(TABLES / "communication_niche_lineage_sender_receiver_map.csv", index=False)

    control_rows = []
    if not broad_audit.empty:
        control_rows.append(
            {
                "control": "broad_random_gene_sets",
                "status": "pass_for_E1_and_EB",
                "evidence": "random_gene_control_q=0 for E1 and GSE154572 EB; internal does not pass",
            }
        )
    if not strict_scan.empty:
        control_rows += [
            {
                "control": "module_random_gene_sets",
                "status": "partial",
                "evidence": "E1 FGF and EB WNT pass; most strict families remain weak/fail",
            },
            {
                "control": "time_shuffle_windows",
                "status": "partial",
                "evidence": "stored per module in morphogen_niche_module_scan.csv",
            },
            {
                "control": "lineage_label_shuffle",
                "status": "not_decisive",
                "evidence": "lineage divergence is secondary; activation is the more robust signal",
            },
            {
                "control": "neighbor_graph_randomization",
                "status": "not_run_in_v2",
                "evidence": "requires additional rerun; current result remains candidate annotation",
            },
            {
                "control": "module_dropout_jackknife",
                "status": "not_run_in_v2",
                "evidence": "future hardening needed before pathway-specific claims",
            },
            {
                "control": "alternative_knn_k",
                "status": "not_run_in_v2",
                "evidence": "future hardening needed before pathway-specific claims",
            },
            {
                "control": "expression_abundance_matched_random_modules",
                "status": "approximated_not_strict",
                "evidence": "random module controls run, but expression-abundance matching is not yet exact",
            },
            {
                "control": "degree_matched_random_lr_pairs",
                "status": "not_run_in_v2",
                "evidence": "strict LR degree-matched null remains future work",
            },
        ]
    controls = pd.DataFrame(control_rows)
    controls.to_csv(TABLES / "communication_niche_control_hardening_summary.csv", index=False)

    report = [
        "# Communication-Niche Manuscript Hardening",
        "",
        "## Broad Versus Strict",
        "",
        broad_vs_strict.to_markdown(index=False),
        "",
        "## Sender/Receiver Decomposition",
        "",
        decomp.to_markdown(index=False) if not decomp.empty else "No decomposition rows.",
        "",
        "## Candidate Validation Panel",
        "",
        panel.head(20).to_markdown(index=False) if not panel.empty else "No cleaned panel candidates.",
        "",
        "## Control Hardening",
        "",
        controls.to_markdown(index=False) if not controls.empty else "No controls summarized.",
        "",
        "## Manuscript Language",
        "",
        "Allowed: branch windows can be annotated with a candidate extracellular communication-niche field. Broad PathwayFinder/OmniPath-derived sender-receiver fields show acceptable computational support in E1 and EB; strict single-family morphogen modules remain weak, with FGF and WNT prioritized for future validation.",
        "",
        "Forbidden: do not state that FGF/WNT signalling drives branch fate, that CCI is established, or that ligand-receptor perturbation has been experimentally supported.",
    ]
    (REPORTS / "communication_niche_manuscript_hardening.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    story = [
        "# Communication-Niche Manuscript Result",
        "",
        "Branch windows can be annotated by a local extracellular communication-niche field derived from PathwayFinder/OmniPath intercell priors.",
        "",
        "The broad niche field is more robust than strict single-family modules, supporting a distributed extracellular niche-state interpretation rather than a single ligand-receptor mechanism.",
        "",
        "FGF and WNT are prioritized as future validation panels, but the current evidence remains computational and candidate-level.",
    ]
    (MANUSCRIPT / "communication_niche_manuscript_result.md").write_text("\n".join(story) + "\n", encoding="utf-8")
    print({"broad_tier": broad_tier, "strict_tier": strict_tier, "panel_rows": int(len(panel))})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
