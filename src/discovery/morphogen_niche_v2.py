from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from src.discovery.communication_niche_search import (
    DATASETS,
    OMNIPATH_DIR,
    REPORTS,
    ROOT,
    TABLES,
    _choose_windows,
    _gene_index,
    _indices_for_prefixes,
    _latent,
    _lineage_key,
    _local_neighbor_mean,
    _module_scan,
    _normalize_score,
    _read_symbol_columns,
    _time_key,
    _vector_for_indices,
    _window_metrics,
)


MANUSCRIPT = ROOT / "manuscript"


@dataclass(frozen=True)
class StrictModule:
    name: str
    ligand_prefixes: tuple[str, ...]
    receptor_prefixes: tuple[str, ...]
    role: str


STRICT_MODULES = [
    StrictModule("FGF_niche", ("FGF",), ("FGFR",), "primary_morphogen"),
    StrictModule("WNT_niche", ("WNT",), ("FZD", "LRP5", "LRP6", "ROR", "RYK"), "primary_morphogen"),
    StrictModule("BMP_niche", ("BMP", "GDF"), ("BMPR", "ACVR"), "primary_morphogen"),
    StrictModule("TGF_NODAL_ACTIVIN_niche", ("TGFB", "NODAL", "INHBA", "INHBB"), ("TGFBR", "ACVR"), "primary_morphogen"),
    StrictModule("SHH_niche", ("SHH", "IHH", "DHH"), ("PTCH", "SMO", "BOC", "CDON"), "primary_morphogen"),
    StrictModule("Notch_Delta_niche", ("DLL", "JAG"), ("NOTCH",), "primary_contact"),
    StrictModule(
        "ECM_adhesion_guidance_niche",
        ("COL", "LAMA", "LAMB", "LAMC", "FN", "THBS", "NTN", "SEMA", "EFN", "MDK", "PTN"),
        ("ITG", "CD47", "DCC", "UNC5", "NRP", "PLXN", "EPH", "NCL", "SDC"),
        "secondary_ecm_guidance",
    ),
    StrictModule(
        "chemokine_growth_factor_niche",
        ("CXCL", "CCL", "KITL", "VEGF", "PDGF", "CSF", "IL"),
        ("CXCR", "CCR", "KIT", "KDR", "FLT", "PDGFR", "CSF", "IL"),
        "secondary_growth",
    ),
]


SUSPICIOUS_PREFIXES = ("HSP", "PCNA", "TP53", "RPL", "RPS", "GAPDH", "ACTB", "TUB", "CALM", "HNRN")
SUSPICIOUS_EXACT = {"HSP90AA1", "HSP90AB1", "PCNA", "TP53", "GAPDH", "ACTB"}


def _load_strict_prior() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    inter = pd.read_csv(OMNIPATH_DIR / "omnipath_intercell.tsv", sep="\t", low_memory=False)
    interactions = pd.read_csv(OMNIPATH_DIR / "omnipath_interactions.tsv", sep="\t", low_memory=False)
    inter["genesymbol"] = inter["genesymbol"].astype(str).str.upper()
    suspicious = inter[
        inter["genesymbol"].isin(SUSPICIOUS_EXACT)
        | inter["genesymbol"].apply(lambda x: any(str(x).startswith(prefix) for prefix in SUSPICIOUS_PREFIXES))
    ][["genesymbol", "category", "transmitter", "receiver", "secreted", "plasma_membrane_transmembrane"]].drop_duplicates()

    annot = inter[
        [
            "uniprot",
            "genesymbol",
            "category",
            "transmitter",
            "receiver",
            "secreted",
            "plasma_membrane_transmembrane",
            "plasma_membrane_peripheral",
        ]
    ].dropna(subset=["uniprot", "genesymbol"]).drop_duplicates()
    src = annot.rename(
        columns={
            "uniprot": "source",
            "genesymbol": "ligand",
            "category": "ligand_category",
            "transmitter": "ligand_transmitter",
            "receiver": "ligand_receiver",
            "secreted": "ligand_secreted",
            "plasma_membrane_transmembrane": "ligand_membrane_tm",
            "plasma_membrane_peripheral": "ligand_membrane_peripheral",
        }
    )
    tgt = annot.rename(
        columns={
            "uniprot": "target",
            "genesymbol": "receptor",
            "category": "receptor_category",
            "transmitter": "receptor_transmitter",
            "receiver": "receptor_receiver",
            "secreted": "receptor_secreted",
            "plasma_membrane_transmembrane": "receptor_membrane_tm",
            "plasma_membrane_peripheral": "receptor_membrane_peripheral",
        }
    )
    joined = interactions.merge(src, on="source", how="inner").merge(tgt, on="target", how="inner")
    curated_ligand_prefixes = tuple(sorted({p for m in STRICT_MODULES for p in m.ligand_prefixes}))
    curated_receptor_prefixes = tuple(sorted({p for m in STRICT_MODULES for p in m.receptor_prefixes}))
    ligand_curated = joined["ligand"].astype(str).apply(lambda x: any(x.startswith(p) for p in curated_ligand_prefixes))
    receptor_curated = joined["receptor"].astype(str).apply(lambda x: any(x.startswith(p) for p in curated_receptor_prefixes))
    ligand_ok = (joined["ligand_secreted"] == True) | (joined["ligand_transmitter"] == True) | ligand_curated
    receptor_ok = (
        (joined["receptor_receiver"] == True)
        | (joined["receptor_membrane_tm"] == True)
        | (joined["receptor_membrane_peripheral"] == True)
        | receptor_curated
    )
    suspicious_pair = (
        joined["ligand"].isin(SUSPICIOUS_EXACT)
        | joined["receptor"].isin(SUSPICIOUS_EXACT)
        | joined["ligand"].apply(lambda x: any(str(x).startswith(prefix) for prefix in SUSPICIOUS_PREFIXES))
        | joined["receptor"].apply(lambda x: any(str(x).startswith(prefix) for prefix in SUSPICIOUS_PREFIXES))
    )
    strict = joined[ligand_ok & receptor_ok & ~suspicious_pair].copy()
    strict = strict[["ligand", "receptor", "is_stimulation", "is_inhibition", "ligand_secreted", "ligand_transmitter", "receptor_receiver", "receptor_membrane_tm", "receptor_membrane_peripheral"]].drop_duplicates()
    return joined, strict, suspicious


def _module_genes_from_prior(strict: pd.DataFrame, module: StrictModule, adata: ad.AnnData) -> tuple[list[int], list[int], pd.DataFrame]:
    pairs = strict[
        strict["ligand"].astype(str).apply(lambda x: any(x.startswith(prefix) for prefix in module.ligand_prefixes))
        & strict["receptor"].astype(str).apply(lambda x: any(x.startswith(prefix) for prefix in module.receptor_prefixes))
    ].drop_duplicates(["ligand", "receptor"])
    index = _gene_index(adata)
    ligand_idx = []
    receptor_idx = []
    rows = []
    for row in pairs.itertuples(index=False):
        lig = str(row.ligand).lower()
        rec = str(row.receptor).lower()
        if lig in index and rec in index:
            ligand_idx.append(index[lig])
            receptor_idx.append(index[rec])
            rows.append({"ligand": str(row.ligand), "receptor": str(row.receptor)})
    return sorted(set(ligand_idx)), sorted(set(receptor_idx)), pd.DataFrame(rows)


def _mean_expr(adata: ad.AnnData, idx: int) -> np.ndarray:
    x = adata.X[:, idx]
    if sparse.issparse(x):
        return np.asarray(x.todense()).ravel()
    return np.asarray(x).ravel()


def _top_sender_receiver_lineage(values: np.ndarray, lineages: np.ndarray) -> str:
    rows = []
    for lab in sorted(set(lineages.astype(str))):
        mask = lineages.astype(str) == lab
        if mask.sum() >= 5:
            rows.append((lab, float(values[mask].mean())))
    if not rows:
        return "unresolved"
    return max(rows, key=lambda x: x[1])[0]


def _candidate_pairs(
    adata: ad.AnnData,
    module: StrictModule,
    pairs: pd.DataFrame,
    times: np.ndarray,
    lineages: np.ndarray,
    event_time: float,
) -> pd.DataFrame:
    if pairs.empty:
        return pd.DataFrame()
    index = _gene_index(adata)
    event_mask = times == event_time
    if event_mask.sum() < 10:
        event_mask = np.ones_like(times, dtype=bool)
    rows = []
    for row in pairs.drop_duplicates(["ligand", "receptor"]).itertuples(index=False):
        lig_key = str(row.ligand).lower()
        rec_key = str(row.receptor).lower()
        if lig_key not in index or rec_key not in index:
            continue
        lig = _normalize_score(_mean_expr(adata, index[lig_key]))
        rec = _normalize_score(_mean_expr(adata, index[rec_key]))
        pair_score = float(np.sqrt(np.maximum(lig[event_mask], 0).mean() * np.maximum(rec[event_mask], 0).mean()))
        if pair_score <= 0:
            continue
        rows.append(
            {
                "module": module.name,
                "ligand": str(row.ligand),
                "receptor": str(row.receptor),
                "event_pair_score": pair_score,
                "sender_lineage_proxy": _top_sender_receiver_lineage(lig[event_mask], lineages[event_mask]),
                "receiver_lineage_proxy": _top_sender_receiver_lineage(rec[event_mask], lineages[event_mask]),
                "known_developmental_relevance": _known_relevance(str(row.ligand), str(row.receptor), module.name),
                "confidence_tier": "candidate",
                "allowed_claim": "candidate spatial or perturbation target",
                "forbidden_claim": "not_a_retained_mechanism",
            }
        )
    return pd.DataFrame(rows).sort_values("event_pair_score", ascending=False).head(10)


def _known_relevance(ligand: str, receptor: str, module: str) -> str:
    if module.startswith("FGF"):
        return "FGF signalling is a canonical developmental patterning and mesoderm/primitive-streak-associated pathway."
    if module.startswith("WNT"):
        return "WNT signalling is a canonical developmental patterning and lineage commitment pathway."
    if module.startswith("BMP"):
        return "BMP/GDF signalling is a canonical germ-layer and patterning pathway."
    if module.startswith("TGF") or ligand.startswith("NODAL"):
        return "TGF/NODAL/ACTIVIN signalling is a canonical early developmental patterning pathway."
    if module.startswith("Notch"):
        return "Notch/Delta signalling is a canonical contact-dependent fate patterning pathway."
    if "ECM" in module:
        return "ECM/adhesion/guidance cues are plausible niche and migration-associated branch-window signals."
    return "developmental relevance should be checked before prioritizing wet-lab validation."


def _module_random_controls(
    adata: ad.AnnData,
    z: np.ndarray,
    times: np.ndarray,
    lineages: np.ndarray,
    pre: float,
    event: float,
    post: float,
    n_ligands: int,
    n_receptors: int,
    observed_activation: float,
    observed_score: float,
    n_permutations: int = 40,
) -> tuple[float, float]:
    rng = np.random.default_rng(91357)
    all_idx = np.arange(adata.n_vars)
    if n_ligands <= 0 or n_receptors <= 0:
        return np.nan, np.nan
    activation_ge = 0
    score_ge = 0
    valid = 0
    for _ in range(n_permutations):
        lig_idx = rng.choice(all_idx, size=min(n_ligands, len(all_idx)), replace=False).tolist()
        rec_idx = rng.choice(all_idx, size=min(n_receptors, len(all_idx)), replace=False).tolist()
        lig = _vector_for_indices(adata, lig_idx)
        rec = _vector_for_indices(adata, rec_idx)
        field = _local_neighbor_mean(z, lig, times) * rec
        df = pd.DataFrame({"time": times, "lineage": lineages, "comm_field": field, "receiver_score": rec}).dropna(
            subset=["time"]
        )
        res = _window_metrics(df, pre, event, post, shuffled=False)
        if res.get("valid"):
            valid += 1
            activation_ge += int(res["communication_activation_effect"] >= observed_activation)
            score_ge += int(res["communication_window_score"] >= observed_score)
    if valid == 0:
        return np.nan, np.nan
    return activation_ge / valid, score_ge / valid


def analyze() -> None:
    TABLES.mkdir(exist_ok=True)
    REPORTS.mkdir(exist_ok=True)
    MANUSCRIPT.mkdir(exist_ok=True)
    joined, strict, suspicious = _load_strict_prior()
    prior_rows = [
        {"metric": "total_omnipath_joined_intercell_edges", "count": int(len(joined))},
        {"metric": "strict_extracellular_edges", "count": int(len(strict))},
        {"metric": "removed_intracellular_or_generic_edges", "count": int(max(len(joined) - len(strict), 0))},
        {"metric": "suspicious_intercell_genes_flagged", "count": int(suspicious["genesymbol"].nunique())},
    ]
    for module in STRICT_MODULES:
        module_edges = strict[
            strict["ligand"].astype(str).apply(lambda x: any(x.startswith(prefix) for prefix in module.ligand_prefixes))
            & strict["receptor"].astype(str).apply(lambda x: any(x.startswith(prefix) for prefix in module.receptor_prefixes))
        ]
        prior_rows.append({"metric": f"{module.name}_strict_edges", "count": int(len(module_edges.drop_duplicates(['ligand', 'receptor'])))})
    pd.DataFrame(prior_rows).to_csv(TABLES / "morphogen_niche_strict_prior_audit.csv", index=False)
    suspicious.head(50).to_csv(TABLES / "morphogen_niche_removed_suspicious_examples.csv", index=False)

    scans = []
    candidates = []
    for spec in DATASETS:
        if not spec.path.exists():
            continue
        adata = ad.read_h5ad(spec.path)
        tkey = _time_key(adata)
        lkey = _lineage_key(adata)
        if tkey is None or lkey is None:
            continue
        times = pd.to_numeric(adata.obs[tkey], errors="coerce").to_numpy()
        lineages = adata.obs[lkey].astype(str).to_numpy()
        z = _latent(adata)
        windows = _choose_windows(times, spec.event_time)
        if not windows:
            continue
        # Prefer the pre-registered window, otherwise scan for the highest module score per dataset.
        pre, event, post, source = windows[0]
        for module in STRICT_MODULES:
            lig_idx, rec_idx, pair_df = _module_genes_from_prior(strict, module, adata)
            if not lig_idx or not rec_idx:
                scans.append(
                    {
                        "dataset": spec.dataset_id,
                        "module": module.name,
                        "module_role": module.role,
                        "valid": False,
                        "matched_ligand_genes": len(lig_idx),
                        "matched_receptor_genes": len(rec_idx),
                        "support_tier": "fail",
                        "failure_reason": "no_strict_ligand_or_receptor_match",
                    }
                )
                continue
            lig = _vector_for_indices(adata, lig_idx)
            rec = _vector_for_indices(adata, rec_idx)
            field = _local_neighbor_mean(z, lig, times) * rec
            df = pd.DataFrame({"time": times, "lineage": lineages, "comm_field": field, "receiver_score": rec}).dropna(
                subset=["time"]
            )
            result = _window_metrics(df, pre, event, post, shuffled=False)
            time_control = _window_metrics(df, pre, event, post, shuffled=True)
            if not result.get("valid"):
                scans.append(
                    {
                        "dataset": spec.dataset_id,
                        "module": module.name,
                        "module_role": module.role,
                        "valid": False,
                        "matched_ligand_genes": len(lig_idx),
                        "matched_receptor_genes": len(rec_idx),
                        "support_tier": "fail",
                        "failure_reason": "window_not_valid",
                    }
                )
                continue
            activation_q, score_q = _module_random_controls(
                adata,
                z,
                times,
                lineages,
                pre,
                event,
                post,
                len(lig_idx),
                len(rec_idx),
                float(result["communication_activation_effect"]),
                float(result["communication_window_score"]),
            )
            time_shuffle_score = float(time_control.get("communication_window_score", np.nan)) if time_control.get("valid") else np.nan
            activation_positive = float(result["communication_activation_effect"]) > 0
            receiver_positive = float(result["receiver_priming_effect"]) > 0
            divergence_positive = float(result["post_event_comm_divergence_effect"]) >= 0
            random_clean = pd.isna(activation_q) or activation_q <= 0.10
            time_clean = pd.isna(time_shuffle_score) or float(result["communication_window_score"]) > time_shuffle_score
            if activation_positive and random_clean and time_clean:
                tier = "acceptable" if spec.dataset_id in {"E1_MouseGastrulationData", "internal_native", "GSE154572_EB_WT"} else "weak"
            elif activation_positive:
                tier = "weak"
            else:
                tier = "fail"
            row = {
                "dataset": spec.dataset_id,
                "role": spec.role,
                "independence_tier": spec.independence_tier,
                "expected_window_type": spec.expected_window_type,
                "module": module.name,
                "module_role": module.role,
                "valid": True,
                "pre_time": pre,
                "event_time": event,
                "post_time": post,
                "window_source": source,
                "matched_ligand_genes": len(lig_idx),
                "matched_receptor_genes": len(rec_idx),
                "strict_pair_count": int(len(pair_df)),
                "activation_effect": float(result["communication_activation_effect"]),
                "receiver_priming_effect": float(result["receiver_priming_effect"]),
                "post_event_divergence_effect": float(result["post_event_comm_divergence_effect"]),
                "communication_window_score": float(result["communication_window_score"]),
                "random_activation_q": activation_q,
                "random_score_q": score_q,
                "time_shuffle_score": time_shuffle_score,
                "activation_positive": bool(activation_positive),
                "receiver_priming_positive": bool(receiver_positive),
                "post_event_divergence_positive": bool(divergence_positive),
                "random_control_pass": bool(random_clean),
                "time_shuffle_pass": bool(time_clean),
                "support_tier": tier,
                "failure_reason": "" if tier != "fail" else "no_activation_or_controls_failed",
            }
            scans.append(row)
            if tier in {"acceptable", "weak"}:
                cands = _candidate_pairs(adata, module, pair_df, times, lineages, event)
                if not cands.empty:
                    cands.insert(0, "dataset", spec.dataset_id)
                    cands["event_time"] = event
                    cands["module_support_tier"] = tier
                    candidates.append(cands)
    scan_df = pd.DataFrame(scans)
    scan_df.to_csv(TABLES / "morphogen_niche_module_scan.csv", index=False)
    cand_df = pd.concat(candidates, ignore_index=True) if candidates else pd.DataFrame()
    if not cand_df.empty:
        cand_df.to_csv(TABLES / "morphogen_niche_candidate_lr_pairs.csv", index=False)

    tier_rows = []
    for module, grp in scan_df[scan_df["valid"].astype(str).str.lower().isin(["true", "1"])].groupby("module"):
        acceptable = grp[grp["support_tier"].eq("acceptable")]
        weak_or_better = grp[grp["support_tier"].isin(["acceptable", "weak"])]
        datasets = weak_or_better["dataset"].astype(str).tolist()
        has_internal = "internal_native" in datasets
        has_e1 = "E1_MouseGastrulationData" in datasets
        has_independent = any(d in datasets for d in ["GSE154572_EB_WT", "E5_zebrafish_Farrell", "E2_GSE212050_gastruloid"])
        if has_internal and has_e1 and has_independent and len(acceptable) >= 2:
            tier = "acceptable"
            interpretation = "strict module survives as a candidate branch-window niche signal"
        elif has_e1 and has_independent:
            tier = "weak"
            interpretation = "module has cross-dataset direction but limited internal or control support"
        elif len(weak_or_better) >= 1:
            tier = "weak"
            interpretation = "dataset-specific candidate only"
        else:
            tier = "fail"
            interpretation = "no reliable module-specific branch-window niche signal"
        tier_rows.append(
            {
                "module": module,
                "tier": tier,
                "support_dataset_count": int(len(weak_or_better)),
                "support_datasets": ";".join(datasets),
                "acceptable_datasets": ";".join(acceptable["dataset"].astype(str).tolist()),
                "has_internal": has_internal,
                "has_e1": has_e1,
                "has_independent": has_independent,
                "mean_activation_effect": float(pd.to_numeric(grp["activation_effect"], errors="coerce").mean()),
                "mean_post_event_divergence_effect": float(pd.to_numeric(grp["post_event_divergence_effect"], errors="coerce").mean()),
                "interpretation": interpretation,
            }
        )
    tier_df = pd.DataFrame(tier_rows).sort_values(["tier", "module"])
    tier_df.to_csv(TABLES / "morphogen_niche_module_tiers.csv", index=False)
    ranked = tier_df.copy()
    if not ranked.empty:
        ranked["tier_rank"] = ranked["tier"].map({"strong": 3, "acceptable": 2, "weak": 1, "fail": 0}).fillna(0)
        ranked["e1_independent_bonus"] = ranked["has_e1"].astype(int) + ranked["has_independent"].astype(int)
        ranked = ranked.sort_values(["tier_rank", "support_dataset_count", "e1_independent_bonus"], ascending=False)
    best = ranked[ranked["tier"].isin(["acceptable", "strong"])].head(1)
    final_tier = "acceptable" if not best.empty else ("weak" if (tier_df["tier"] == "weak").any() else "fail")
    strongest = str(best.iloc[0]["module"]) if not best.empty else (
        str(ranked[ranked["tier"].eq("weak")].iloc[0]["module"]) if (ranked["tier"] == "weak").any() else "none"
    )

    report = [
        "# Morphogen Communication-Niche v2",
        "",
        "## Strict Prior Audit",
        "",
        pd.DataFrame(prior_rows).to_markdown(index=False),
        "",
        "## Cross-Dataset Module Tiers",
        "",
        tier_df.to_markdown(index=False) if not tier_df.empty else "No valid modules.",
        "",
        "## Final Interpretation",
        "",
        f"- final_communication_niche_tier: `{final_tier}`",
        f"- strongest_module: `{strongest}`",
        "- allowed_claim: strict extracellular morphogen/communication-niche priming is a candidate branch-window annotation.",
        "- forbidden_claim: do not describe this as confirmed signalling, communication-driven fate control, or experimental validation.",
        "",
        "## Future Spatial/Perturbation Design",
        "",
        "Primary hypothesis: morphogen communication-niche priming occurs at the branch window and can be spatially observed as local ligand-producing sender neighborhoods around receptor-competent receiver cells.",
        "Required assay: spatial transcriptomics, MERFISH, seqFISH or smFISH in a gastruloid/embryoid-body time course with branch-window stages, ligand/receptor readouts, cell type labels and optional FGF/WNT/BMP/TGF/NODAL perturbations.",
        "Primary readouts: sender ligand density, receiver receptor competence, local niche field, branch-window order parameter and post-event lineage divergence.",
    ]
    (REPORTS / "morphogen_niche_v2_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    story = [
        "# Morphogen Communication-Niche Hypothesis",
        "",
        f"Final tier: `{final_tier}`.",
        "",
        f"The strongest strict extracellular module is `{strongest}`. The result should be framed as a candidate branch-window niche annotation, not as established ligand-receptor biology.",
        "",
        "This direction preserves the SwarmLineage-OT methodology: native OT defines developmental flow, finite-agent rollout defines branch-window order parameters, and PathwayFinder/OmniPath intercell priors annotate local sender--receiver niche fields.",
    ]
    (MANUSCRIPT / "morphogen_niche_hypothesis.md").write_text("\n".join(story) + "\n", encoding="utf-8")
    print(json.dumps({"final_tier": final_tier, "strongest_module": strongest, "strict_edges": int(len(strict))}, indent=2))


def main() -> int:
    analyze()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
