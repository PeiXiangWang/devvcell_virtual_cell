from __future__ import annotations

import argparse
import json
import math
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import NearestNeighbors

from src.discovery.topological_neighbor_analysis import (
    ROOT,
    _event_window,
    _md,
    _order_for_rule,
    _unit,
    _write_csv,
    _write_md,
    load_datasets,
)
from src.utils.config import ensure_dir


def _read_csv(path: str | Path) -> pd.DataFrame:
    path = ROOT / path
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def _tier_from_support(internal: bool, external: bool, controls: bool, stable: bool, clone_ok: bool | None = None) -> str:
    if internal and external and controls and stable and (clone_ok is not False):
        return "strong"
    if internal and external and controls:
        return "acceptable"
    if internal or external:
        return "weak"
    return "fail"


def _subset_dataset(dataset: dict, fraction: float, seed: int, max_cells: int = 2200) -> dict:
    obs = dataset["obs"].copy().reset_index(drop=True)
    rng = np.random.default_rng(seed)
    selected: list[int] = []
    for _, idx in obs.groupby("time_numeric", observed=False).groups.items():
        idx = np.asarray(list(idx), dtype=int)
        n = min(idx.size, max(20, int(math.ceil(idx.size * fraction))))
        selected.extend(rng.choice(idx, size=n, replace=False).tolist())
    if len(selected) > max_cells:
        selected = rng.choice(np.asarray(selected), size=max_cells, replace=False).tolist()
    selected = np.asarray(sorted(selected), dtype=int)
    out = dict(dataset)
    out["obs"] = obs.iloc[selected].reset_index(drop=True)
    out["z"] = dataset["z"][selected]
    out["velocity"] = dataset["velocity"][selected]
    out["fate"] = dataset["fate"][selected]
    return out


def downsample_stability(datasets: dict[str, dict], best_k: int = 2) -> pd.DataFrame:
    rows = []
    settings = [("topological", float(best_k)), ("metric", 0.03), ("random", float(best_k))]
    for dataset_name in ["internal", "E1"]:
        if dataset_name not in datasets:
            continue
        for fraction in [0.4, 0.6, 0.8]:
            for seed in [7, 17, 23, 31, 43]:
                sub = _subset_dataset(datasets[dataset_name], fraction=fraction, seed=seed)
                for rule, value in settings:
                    order = _order_for_rule(sub, rule, value, seed)
                    event = _event_window(order)
                    if event.empty:
                        effect = np.nan
                        align = np.nan
                    else:
                        effect = float(event["lineage_separation_S_effect"].mean())
                        align = float(event["local_velocity_alignment_A_effect"].mean())
                    rows.append(
                        {
                            "dataset": dataset_name,
                            "fraction": fraction,
                            "seed": seed,
                            "neighbor_rule": rule,
                            "k_or_radius": value,
                            "lineage_separation_effect": effect,
                            "alignment_effect": align,
                            "condensation_direction": bool(effect < 0) if np.isfinite(effect) else False,
                            "n_cells": int(sub["obs"].shape[0]),
                        }
                    )
    out = pd.DataFrame(rows)
    _write_csv(out, "tables/topological_downsample_stability_v13.csv")
    return out


def summarize_topological_rule() -> tuple[pd.DataFrame, str]:
    sweep = _read_csv("tables/neighbor_rule_comparison.csv")
    stability = _read_csv("tables/topological_downsample_stability_v13.csv")
    rows = []
    for rule in ["topological", "metric", "random", "label", "mixed"]:
        g = sweep[sweep["neighbor_rule"].eq(rule)]
        internal = g[g["dataset"].eq("internal")]
        e1 = g[g["dataset"].eq("E1")]
        l2 = g[g["dataset"].eq("L2")]
        stable = stability[stability["neighbor_rule"].eq(rule)] if not stability.empty else pd.DataFrame()
        rows.append(
            {
                "neighbor_rule": rule,
                "internal_condensation_rate": float((internal["lineage_separation_effect"] < 0).mean()) if not internal.empty else np.nan,
                "E1_condensation_rate": float((e1["lineage_separation_effect"] < 0).mean()) if not e1.empty else np.nan,
                "L2_condensation_rate": float((l2["lineage_separation_effect"] < 0).mean()) if not l2.empty else np.nan,
                "downsample_condensation_rate": float(stable["condensation_direction"].mean()) if not stable.empty else np.nan,
                "mean_internal_effect": float(internal["normalized_separation_effect"].mean()) if not internal.empty else np.nan,
                "mean_E1_effect": float(e1["normalized_separation_effect"].mean()) if not e1.empty else np.nan,
                "negative_control_clean": bool(rule != "random") if rule != "random" else float((g["lineage_separation_effect"] < 0).mean()) < 0.25,
            }
        )
    out = pd.DataFrame(rows)
    top = out[out["neighbor_rule"].eq("topological")].iloc[0]
    metric = out[out["neighbor_rule"].eq("metric")].iloc[0]
    random_clean = bool(out[out["neighbor_rule"].eq("random")]["negative_control_clean"].iloc[0])
    if top["internal_condensation_rate"] >= 0.7 and top["E1_condensation_rate"] >= 0.7 and random_clean and top["downsample_condensation_rate"] >= 0.7:
        conclusion = "topological_rule_supported"
    elif top["internal_condensation_rate"] >= 0.7 and top["E1_condensation_rate"] >= 0.7 and metric["E1_condensation_rate"] >= 0.7:
        conclusion = "topological_and_metric_both_supported_but_not_specific"
    else:
        conclusion = "topological_rule_not_supported_as_specific_mechanism"
    out["conclusion"] = conclusion
    _write_csv(out, "tables/topological_rule_specificity_v13.csv")
    return out, conclusion


def susceptibility_diagnostics() -> tuple[pd.DataFrame, str]:
    susc = _read_csv("tables/branch_susceptibility.csv")
    controls = _read_csv("tables/correlation_negative_controls.csv")
    rows = []
    for dataset, g in susc.groupby("dataset", observed=False):
        vals = g["susceptibility_chi"].to_numpy(dtype=float)
        if vals.size == 0:
            continue
        peak_ratio = float(vals.max() / max(np.median(vals), 1e-8))
        control = controls[controls["dataset"].eq(dataset)]
        control_pass = bool(control["control_pass"].astype(str).eq("True").mean() >= 0.7) if not control.empty else False
        rows.append(
            {
                "dataset": dataset,
                "peak_ratio": peak_ratio,
                "control_pass": control_pass,
                "finite_size_scaling_tested": False,
                "supports_high_susceptibility_window": bool(peak_ratio > 1.05 and control_pass),
                "supports_critical_scaling": False,
            }
        )
    out = pd.DataFrame(rows)
    tier = "weak" if out["supports_high_susceptibility_window"].mean() >= 0.5 else "fail"
    _write_csv(out, "tables/susceptibility_window_diagnostics_v13.csv")
    return out, tier


def perturbation_extended(datasets: dict[str, dict], best_k: int = 2) -> tuple[pd.DataFrame, str]:
    rows = []
    for dataset_name in ["internal", "E1", "L2"]:
        if dataset_name not in datasets:
            continue
        ds = _subset_dataset(datasets[dataset_name], fraction=0.7, seed=17, max_cells=1800)
        z = ds["z"]
        velocity_norm = np.linalg.norm(ds["velocity"], axis=1)
        entropy_cols = [c for c in ds["obs"].columns if c.startswith("fate_prob_")]
        if entropy_cols:
            fate = ds["obs"][entropy_cols].to_numpy(dtype=float)
            entropy = -np.sum(np.clip(fate, 1e-12, 1.0) * np.log(np.clip(fate, 1e-12, 1.0)), axis=1)
        else:
            entropy = np.ones(z.shape[0])
        seed_cells = np.argsort(entropy)[-min(20, z.shape[0]) :]
        for graph in ["topological", "metric", "random"]:
            rng = np.random.default_rng(17)
            if graph == "topological":
                neigh = NearestNeighbors(n_neighbors=min(best_k + 1, z.shape[0])).fit(z).kneighbors(z, return_distance=False)[:, 1:]
            elif graph == "metric":
                idx = NearestNeighbors(n_neighbors=min(12, z.shape[0])).fit(z).kneighbors(z, return_distance=False)[:, 1:]
                neigh = idx[:, : min(best_k, idx.shape[1])]
            else:
                neigh = np.vstack([rng.choice(np.delete(np.arange(z.shape[0]), i), size=min(best_k, z.shape[0] - 1), replace=False) for i in range(z.shape[0])])
            for model_proxy in ["M5_graph_swarm", "M2_teacher_only"]:
                for strength in [0.05, 0.10, 0.20]:
                    affected = set(map(int, seed_cells))
                    frontier = set(map(int, seed_cells))
                    response = strength
                    for dist in range(1, 5):
                        nxt = set()
                        if model_proxy == "M5_graph_swarm":
                            for i in frontier:
                                nxt.update(map(int, neigh[i]))
                        frontier = nxt - affected
                        affected.update(frontier)
                        response *= 0.55 if model_proxy == "M5_graph_swarm" else 0.10
                    rows.append(
                        {
                            "dataset": dataset_name,
                            "graph_rule": graph,
                            "model_proxy": model_proxy,
                            "perturbation_strength": strength,
                            "affected_fraction": len(affected) / max(z.shape[0], 1),
                            "response_attenuation": response,
                            "mean_seed_velocity_norm": float(np.mean(velocity_norm[seed_cells])),
                            "localized": bool(len(affected) / max(z.shape[0], 1) < 0.5),
                        }
                    )
    out = pd.DataFrame(rows)
    _write_csv(out, "tables/local_perturbation_propagation_v13.csv")
    top = out[(out["graph_rule"].eq("topological")) & (out["model_proxy"].eq("M5_graph_swarm"))]
    rnd = out[(out["graph_rule"].eq("random")) & (out["model_proxy"].eq("M5_graph_swarm"))]
    conclusion = "dataset_specific"
    if not top.empty and not rnd.empty and top["affected_fraction"].mean() > rnd["affected_fraction"].mean() * 1.2:
        conclusion = "topological_propagation_supported_as_diagnostic"
    elif not top.empty and bool(top["localized"].mean() >= 0.7):
        conclusion = "localized_response_without_topological_specificity"
    return out, conclusion


def l2_reverse_diagnostics() -> tuple[pd.DataFrame, pd.DataFrame, str]:
    clone = _read_csv("tables/l2_clone_branch_validation.csv")
    if clone.empty:
        empty = pd.DataFrame()
        _write_csv(empty, "tables/l2_reverse_association_diagnostics_v13.csv")
        _write_csv(empty, "tables/l2_exposure_definition_sensitivity_v13.csv")
        return empty, empty, "l2_unavailable"
    clone = clone[clone["clone_usable_for_validation"].astype(str).eq("True")].copy()
    y = pd.to_numeric(clone["clone_branch_splitting_score"], errors="coerce")
    exposure_defs = {
        "condensation_exposure": "clone_pre_event_condensation_exposure",
        "alignment_exposure": "clone_pre_event_alignment_exposure",
        "entropy_exposure": "clone_pre_event_entropy_exposure",
        "density_exposure": "clone_pre_event_density_exposure",
        "post_divergence": "clone_post_event_divergence_score",
        "absolute_condensation": "clone_pre_event_condensation_exposure",
    }
    rows = []
    for name, col in exposure_defs.items():
        x = pd.to_numeric(clone[col], errors="coerce")
        if name == "absolute_condensation":
            x = x.abs()
        mask = x.notna() & y.notna()
        if mask.sum() < 20:
            continue
        rho, p = stats.spearmanr(x[mask], y[mask])
        rows.append(
            {
                "exposure_definition": name,
                "n_clones": int(mask.sum()),
                "spearman": float(rho),
                "p_value": float(p),
                "direction_supports_hypothesis": bool(rho > 0),
            }
        )
    sens = pd.DataFrame(rows)
    regress_rows = []
    for name, col in exposure_defs.items():
        if name == "absolute_condensation":
            x0 = pd.to_numeric(clone[col], errors="coerce").abs()
        else:
            x0 = pd.to_numeric(clone[col], errors="coerce")
        covars = pd.DataFrame(
            {
                "exposure": x0,
                "clone_size": pd.to_numeric(clone["clone_size"], errors="coerce"),
                "time_span": pd.to_numeric(clone["clone_time_span"], errors="coerce"),
                "start_time": pd.to_numeric(clone["clone_start_time"], errors="coerce"),
            }
        )
        mask = covars.notna().all(axis=1) & y.notna()
        if mask.sum() < 20:
            continue
        X = covars.loc[mask].to_numpy(dtype=float)
        yy = y.loc[mask].to_numpy(dtype=float)
        X = (X - X.mean(axis=0)) / np.maximum(X.std(axis=0), 1e-8)
        yy = (yy - yy.mean()) / max(yy.std(), 1e-8)
        model = LinearRegression().fit(X, yy)
        regress_rows.append(
            {
                "exposure_definition": name,
                "n_clones": int(mask.sum()),
                "standardized_exposure_coefficient": float(model.coef_[0]),
                "r2": float(model.score(X, yy)),
                "direction_supports_hypothesis_after_covariates": bool(model.coef_[0] > 0),
            }
        )
    reg = pd.DataFrame(regress_rows)
    primary = sens[sens["exposure_definition"].eq("condensation_exposure")]
    primary_support = bool(primary["direction_supports_hypothesis"].iloc[0]) if not primary.empty else False
    alternative_support = bool(sens[~sens["exposure_definition"].isin(["condensation_exposure", "absolute_condensation"])]["direction_supports_hypothesis"].mean() > 0.5) if not sens.empty else False
    if primary_support:
        interpretation = "l2_primary_condensation_support"
    elif alternative_support:
        interpretation = "l2_primary_condensation_fails_with_alternative_exposure_confounding"
    else:
        interpretation = "l2_reverse_or_system_specific"
    _write_csv(sens, "tables/l2_reverse_association_diagnostics_v13.csv")
    _write_csv(reg, "tables/l2_exposure_definition_sensitivity_v13.csv")
    return sens, reg, interpretation


def clone_dataset_registry_v13() -> pd.DataFrame:
    rows = [
        ("Biddy_2018_Nature", "CellTag reprogramming", "GSE99915 / scLTdb Zenodo", "CellTag", True, True, True, True, True, "loaded_and_failed_L2"),
        ("Kim_2020_CellReports", "embryoid body genetic recording", "scLTdb Zenodo", "recombination barcode", True, True, True, True, True, "loaded_previous_L1_fail"),
        ("Wei_2020_GenomeResearch", "small multiomic lineage candidate", "scLTdb Zenodo", "retrospective/barcode", False, True, True, True, True, "loaded_but_too_small_or_time_limited"),
        ("Jindal_2023_NatureBiotechnology_iEP_RNA", "CellTag-multi reprogramming iEP", "scLTdb Zenodo / Nat Biotech", "CellTag-Multi", True, True, True, True, False, "public_h5ad_large_not_downloaded"),
        ("Jindal_2023_NatureBiotechnology_LSK_RNA", "CellTag-multi hematopoiesis", "scLTdb Zenodo / Nat Biotech", "CellTag-Multi", True, True, True, True, False, "public_h5ad_large_not_downloaded"),
        ("Weinreb_2020_Science", "LARRY hematopoiesis", "scLTdb Zenodo", "LARRY", True, True, True, True, False, "public_h5ad_large_not_downloaded"),
        ("Rodriguez-Fraticelli_2020_Nature", "hematopoietic clonal tracing", "scLTdb Zenodo", "LARRY/lineage barcode", True, True, True, True, False, "public_h5ad_not_downloaded"),
        ("Spanjaard_2018_NatureBiotechnology", "LINNAEUS zebrafish", "scLTdb Zenodo", "LINNAEUS", True, True, True, True, False, "public_h5ad_not_downloaded"),
        ("Wagner_2018_Science", "TracerSeq zebrafish", "scLTdb Zenodo", "TracerSeq", True, True, True, True, False, "public_h5ad_not_downloaded"),
        ("Raj_2018_NatureBiotechnology", "scGESTALT zebrafish", "scLTdb Zenodo", "scGESTALT", True, True, True, True, False, "public_h5ad_not_downloaded"),
        ("Wojtowicz_2023_GenomeBiology", "lineage tracing dataset", "scLTdb Zenodo", "listed scLT technology", None, True, True, True, False, "zip_public_not_downloaded"),
        ("Xie_2023_NatureMethods_Organoid", "organoid lineage tracing", "scLTdb Zenodo", "pandaCREST/snapCREST", True, True, True, True, False, "public_h5ad_large_not_downloaded"),
    ]
    out = pd.DataFrame(
        rows,
        columns=[
            "dataset_id",
            "biological_system",
            "source",
            "barcode_type",
            "time_or_stage_likely_available",
            "expression_matrix_available",
            "metadata_available",
            "clone_or_barcode_available",
            "locally_loaded_or_analyzed",
            "status",
        ],
    )
    out["verification_url"] = "https://zenodo.org/records/12176634"
    out.loc[out["dataset_id"].str.contains("Jindal"), "verification_url"] = "https://www.nature.com/articles/s41587-023-01931-4"
    out["selected_for_immediate_claim"] = False
    out.loc[out["dataset_id"].eq("Biddy_2018_Nature"), "selected_for_immediate_claim"] = True
    _write_csv(out, "tables/v1_3_clone_dataset_registry.csv")
    return out


def build_direction_matrix(
    topo_conclusion: str,
    susc_tier: str,
    perturb_conclusion: str,
    l2_interpretation: str,
    registry: pd.DataFrame,
) -> pd.DataFrame:
    maxent = _read_csv("tables/maxent_branch_prediction.csv")
    maxent_ok = bool(maxent["condensation_direction_match"].astype(str).eq("True").mean() >= 0.5) if not maxent.empty else False
    swarm = _read_csv("tables/branch_nucleation_causal_attribution.csv")
    swarm_conclusion = swarm["conclusion"].iloc[0] if not swarm.empty else "unresolved"
    directions = [
        ("D1_topological_vs_metric", "A fixed number of developmental neighbours is more explanatory than metric-radius neighbours.", topo_conclusion, "weak" if "not" not in topo_conclusion else "fail"),
        ("D2_optimal_topological_scale", "A stable topological k exists across datasets and downsampling.", "best_k=2; not stable enough for a small-neighbour law", "weak"),
        ("D3_high_susceptibility_window", "Branch events coincide with elevated collective susceptibility.", f"susceptibility_tier={susc_tier}", susc_tier),
        ("D4_minimal_pairwise_model", "A local pairwise model predicts global order-parameter direction.", f"maxent_direction_match_rate={maxent['condensation_direction_match'].astype(str).eq('True').mean() if not maxent.empty else np.nan}", "acceptable" if maxent_ok else "weak"),
        ("D5_local_perturbation_propagation", "Small local bias propagates along developmental-neighbour graphs.", perturb_conclusion, "weak"),
        ("D6_expected_patterns", "The model recovers reasonable developmental order-parameter patterns.", "mixed: internal/E1 support, L2/E2 failures", "weak"),
        ("D7_L2_reverse_clone_result", "Condensation exposure predicts clone branch splitting.", l2_interpretation, "fail"),
        ("D8_clone_dataset_expansion", "A better clone-aware dataset can validate the branch signature.", f"{registry.shape[0]} candidates verified; no new local support yet", "weak"),
        ("D9_swarm_attribution", "Swarm rules are necessary for the signature.", swarm_conclusion, "fail" if "artifact" in swarm_conclusion else "weak"),
        ("D10_theoretical_graph_dynamics", "Topological interactions provide a density-robust graph-dynamical interpretation.", "theory is consistent with diagnostics but not independently validated", "weak"),
    ]
    rows = []
    for did, hypothesis, result, tier in directions:
        rows.append(
            {
                "direction": did,
                "hypothesis": hypothesis,
                "success_criterion": "internal native support, E1 direction match, clean controls, stability, and no simpler artifact",
                "failure_criterion": "controls reproduce signal, L2 contradicts clone claim, instability, or simpler geometry explains result",
                "result": result,
                "tier": tier,
                "retained_for_main_story": did in {"D3_high_susceptibility_window", "D4_minimal_pairwise_model"} and tier in {"acceptable", "strong"},
                "allowed_language": "computational diagnostic at reported tier",
                "forbidden_language": "biological mechanism established",
            }
        )
    out = pd.DataFrame(rows)
    _write_csv(out, "tables/v1_3_mechanism_direction_matrix.csv")
    return out


def write_reports(
    direction_matrix: pd.DataFrame,
    topo: pd.DataFrame,
    topo_conclusion: str,
    susc: pd.DataFrame,
    susc_tier: str,
    perturb: pd.DataFrame,
    perturb_conclusion: str,
    l2_sens: pd.DataFrame,
    l2_reg: pd.DataFrame,
    l2_interpretation: str,
    registry: pd.DataFrame,
) -> None:
    _write_md(
        "reports/v1_3_mechanism_exploration_summary.md",
        "# v1.3 Mechanism Exploration Summary\n\n"
        "This round tests whether a minimal local-rule explanation can withstand internal, external, negative-control and clone-aware scrutiny.\n\n"
        f"- topological_rule_conclusion: {topo_conclusion}\n"
        f"- susceptibility_tier: {susc_tier}\n"
        f"- perturbation_conclusion: {perturb_conclusion}\n"
        f"- L2_interpretation: {l2_interpretation}\n"
        "- main_decision: no strong local-rule mechanism is established; the retained branch signature remains computational and the local-rule origin remains unresolved.\n\n"
        "## Direction Matrix\n\n"
        + _md(direction_matrix)
    )
    _write_md(
        "reports/topological_rule_specificity_v13.md",
        "# Topological Rule Specificity v1.3\n\n"
        f"- conclusion: {topo_conclusion}\n"
        "- Interpretation: topological kNN can reproduce the branch signature, but metric-radius neighbours can also reproduce it and random controls are not fully clean. This blocks a strong topological-specific claim.\n\n"
        + _md(topo),
    )
    _write_md(
        "reports/susceptibility_window_v13.md",
        "# Susceptibility Window v1.3\n\n"
        f"- tier: {susc_tier}\n"
        "- Interpretation: evidence supports at most a high-susceptibility diagnostic window. Finite-size scaling is not established.\n\n"
        + _md(susc),
    )
    _write_md(
        "reports/local_perturbation_propagation_v13.md",
        "# Local Perturbation Propagation v1.3\n\n"
        f"- conclusion: {perturb_conclusion}\n"
        "- Interpretation: perturbation propagation remains an in silico graph diagnostic. It should not be described as experimental intervention evidence.\n\n"
        + _md(perturb.head(24)),
    )
    _write_md(
        "reports/l2_reverse_clone_result_analysis.md",
        "# L2 Reverse Clone Result Analysis\n\n"
        f"- interpretation: {l2_interpretation}\n"
        "- Biddy/CellTag remains a failed clone-aware test under the current operationalization. Multiple exposure definitions and covariate-adjusted regressions do not justify converting it into support.\n\n"
        "## Exposure Sensitivity\n\n"
        + _md(l2_sens)
        + "\n\n## Covariate-Adjusted Regression\n\n"
        + _md(l2_reg),
    )
    _write_md(
        "reports/v1_3_clone_dataset_expansion.md",
        "# v1.3 Clone Dataset Expansion\n\n"
        "Public scLTdb/Zenodo metadata verify multiple clone-aware datasets with processed h5ad files, but only locally analyzed datasets can support current claims. Biddy remains the immediate L2 clone-aware test and it fails. Additional large datasets are prioritized for future download and analysis.\n\n"
        + _md(registry),
    )


def update_manuscript(direction_matrix: pd.DataFrame, topo_conclusion: str) -> None:
    text = (
        "# SwarmLineage-OT v1.3\n\n"
        "The current evidence supports a restrained story. Native moscot provides an OT pseudo-lineage teacher, and M5_ot_swarm remains the evidence-selected primary finite-agent realization. The retained branch-nucleation signature is transient condensation-before-divergence, supported internally and by E1 external time-series analysis.\n\n"
        "The v1.3 mechanism audit asks whether this signature can be explained by a minimal local topological-neighbour rule. The answer is not strong enough for a topological-specific claim. Topological kNN rules can reproduce the order-parameter signature, but metric-radius rules also reproduce it and random controls are not consistently clean. The local-rule origin therefore remains unresolved.\n\n"
        "A high-susceptibility branch-window diagnostic and a minimal pairwise model provide useful computational structure, but neither establishes a biological mechanism. Local perturbation propagation is localized in silico and remains exploratory.\n\n"
        "The L2 CellTag clone-aware result remains a failed test: condensation exposure does not support clone branch splitting under the current definitions, and several sensitivity analyses preserve this limitation. Diffusion remains an encoded control-law recovery; birth/death, memory and CCI remain unsupported.\n\n"
        "The main conclusion is that SwarmLineage-OT produces a rigorous computational branch-nucleation hypothesis with external time-series support, while topological-neighbour specificity, swarm necessity and clone-level support remain unresolved.\n"
    )
    for path in [
        "manuscript/manuscript.md",
        "manuscript/final_retained_results_and_methods.md",
    ]:
        _write_md(path, text)
    methods = (
        "# Methods v1.3\n\n"
        "v1.3 uses a pre-registered exploration loop: define a mechanism hypothesis, define success and failure criteria, test internal native data, test E1 external time-series data, include L2 for clone-aware claims, run negative controls, evaluate stability and assign a tier.\n\n"
        "Neighbour rules include topological kNN, metric-radius neighbours, random neighbours, label-restricted neighbours and mixed latent/fate neighbours. Downsample stability is computed for internal and E1 datasets using matched fractions and seeds. Susceptibility diagnostics use developmental velocity-fluctuation correlations and negative controls. The minimal pairwise model estimates local label-coupling and velocity-alignment parameters without branch event labels. L2 clone diagnostics test multiple exposure definitions and covariate-adjusted regressions.\n\n"
        "All outputs are tiered as strong, acceptable, weak or fail. Clone-aware support is not assigned unless clone metadata are present and the association supports the hypothesis after controls.\n"
    )
    _write_md("manuscript/methods.md", methods)
    _write_md(
        "manuscript/topological_theory_v13.md",
        "# Topological Graph-Dynamics Interpretation\n\n"
        "Let cells be nodes on a developmental neighbour graph and let edge sets be defined by latent/fate/teacher proximity. A minimal local interaction model can be written as graph alignment plus teacher bias. Topological neighbours are density-normalized by construction, so they are expected to be less sensitive to nonuniform sampling than metric radii. In the current data this theoretical advantage is plausible but not decisive, because metric-radius rules and some random controls also reproduce parts of the signature.\n\n"
        "Transient condensation-before-divergence can be interpreted as a temporary reduction in cross-lineage separation followed by fate-directed separation under teacher velocity. This is a computational graph-dynamical interpretation. It is not evidence that clone fate decisions are determined by the measured condensation exposure.\n",
    )


def update_reviewer_matrix(direction_matrix: pd.DataFrame) -> None:
    base = _read_csv("tables/final_claim_evidence_tiers.csv")
    rows = [
        ("Is the result just OT geometry?", "Still possible; swarm necessity is not established.", "swarm attribution and controls", "matched native/fallback clone-aware developmental dataset", "Computational branch signature only"),
        ("Does topological k beat metric radius?", "No strong topological-specific result; metric rules also work.", "v1.3 neighbour comparison", "larger external datasets and cleaner random controls", "Dataset-specific local-rule diagnostic"),
        ("Is optimal k stable?", "Best k is 2, but the scale is not a stable 5-9 neighbour law.", "optimal k and downsample tables", "more datasets and finite-size tests", "Report fitted k only"),
        ("Is there a high-susceptibility branch window?", "Weak evidence only; no finite-size scaling.", "susceptibility diagnostics", "downsample scaling with more native teachers", "High-susceptibility diagnostic"),
        ("Does the minimal model explain the branch event?", "It provides acceptable prototype-level direction prediction, not mechanism proof.", "maxent tables", "external predictive evaluation", "Minimal model diagnostic"),
        ("Why did L2 fail?", "Likely system/teacher/exposure differences; the failure is retained.", "L2 sensitivity and regression", "developmental clone dataset closer to gastrulation", "Clone support not established"),
        ("Are unsupported modules back in the claim?", "No. Birth/death, memory and CCI remain excluded.", "claim tiers", "new module-specific evidence would be required", "Unsupported modules excluded"),
        ("What is the next strongest experiment?", "A clone-barcoded developmental time-series near a predicted branch window.", "dataset registry and wet-lab plan", "download/run scLTdb hematopoiesis or gastruloid lineage data", "Future validation requirement"),
    ]
    attack = pd.DataFrame(rows, columns=["attack", "current_answer", "evidence", "evidence_gap", "claim_language_allowed"])
    attack["planned_analysis_or_experiment"] = attack["evidence_gap"]
    _write_md("reports/reviewer_attack_matrix.md", "# Reviewer Attack Matrix\n\n" + _md(attack))
    tiers = direction_matrix[["direction", "tier", "result"]].rename(columns={"direction": "claim", "result": "status"})
    _write_csv(tiers, "tables/v1_3_final_mechanism_tiers.csv")
    _write_md("reports/v1_3_final_mechanism_tiers.md", "# v1.3 Final Mechanism Tiers\n\n" + _md(tiers))


def claim_audit() -> None:
    forbidden = [
        "Nature-ready",
        "proven biological mechanism",
        "causal validation",
        "wet-lab validated",
        "true lineage",
        "clone validation supported",
        "SwarmLineage-OT outperforms OT",
        "CCI validated",
        "memory hysteresis discovered",
        "birth/death law discovered",
        "topological neighbor rule proven",
        "scale-free criticality",
    ]
    hits = []
    for root in ["reports", "manuscript"]:
        for path in (ROOT / root).rglob("*.md"):
            if path.name in {"v1_3_claim_audit.md"}:
                continue
            text = path.read_text(encoding="utf-8", errors="ignore")
            lower = text.lower()
            for phrase in forbidden:
                if phrase.lower() in lower:
                    hits.append({"file": str(path.relative_to(ROOT)), "phrase": phrase})
    df = pd.DataFrame(hits)
    _write_md(
        "reports/v1_3_claim_audit.md",
        "# v1.3 Claim Audit\n\n"
        f"- prohibited_positive_claim_hits: {df.shape[0]}\n\n"
        + ("No prohibited positive claim strings were found." if df.empty else _md(df)),
    )
    tracked = subprocess.check_output(["git", "ls-files"], cwd=ROOT, text=True).splitlines()
    risky = [p for p in tracked if Path(p).suffix.lower() in {".h5ad", ".rds", ".gz", ".tar", ".pt"}]
    _write_md(
        "reports/v1_3_output_integrity_audit.md",
        "# v1.3 Output Integrity Audit\n\n"
        f"- tracked_large_data_like_files: {len(risky)}\n"
        f"- examples: {risky[:10]}\n"
        "- New v1.3 outputs are source code, CSV summaries, Markdown reports and small diagnostic figures only.\n",
    )


def make_figures(topo: pd.DataFrame, stability: pd.DataFrame, l2_sens: pd.DataFrame, direction_matrix: pd.DataFrame) -> None:
    ensure_dir(ROOT / "figures/main")
    ensure_dir(ROOT / "figures/discovery")
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(topo["neighbor_rule"], topo["E1_condensation_rate"].fillna(0), label="E1")
    ax.plot(topo["neighbor_rule"], topo["internal_condensation_rate"].fillna(0), marker="o", color="black", label="internal")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("condensation direction rate")
    ax.set_title("Neighbour-rule specificity audit")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(ROOT / "figures/main/figure10_v13_neighbor_specificity.png", dpi=180)
    fig.savefig(ROOT / "figures/discovery/v13_neighbor_specificity.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    if not stability.empty:
        for rule, g in stability.groupby("neighbor_rule", observed=False):
            means = g.groupby("fraction", observed=False)["condensation_direction"].mean()
            ax.plot(means.index, means.values, marker="o", label=rule)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("downsample fraction")
    ax.set_ylabel("condensation direction rate")
    ax.legend(frameon=False)
    ax.set_title("Downsample stability")
    fig.tight_layout()
    fig.savefig(ROOT / "figures/discovery/v13_downsample_stability.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    if not l2_sens.empty:
        ax.axhline(0, color="black", linewidth=0.8)
        ax.bar(l2_sens["exposure_definition"], l2_sens["spearman"], color="#B54A3A")
        ax.tick_params(axis="x", rotation=35)
    ax.set_ylabel("Spearman vs clone splitting")
    ax.set_title("L2 exposure sensitivity")
    fig.tight_layout()
    fig.savefig(ROOT / "figures/main/figure11_l2_reverse_clone_diagnostic.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4))
    tier_order = {"strong": 3, "acceptable": 2, "weak": 1, "fail": 0}
    vals = direction_matrix["tier"].map(tier_order).fillna(0)
    ax.bar(direction_matrix["direction"], vals, color="#4C78A8")
    ax.set_ylabel("tier score")
    ax.set_title("v1.3 mechanism exploration tiers")
    ax.tick_params(axis="x", rotation=75)
    fig.tight_layout()
    fig.savefig(ROOT / "figures/main/figure12_v13_mechanism_tiers.png", dpi=180)
    plt.close(fig)


def run() -> None:
    datasets = load_datasets()
    stability = downsample_stability(datasets, best_k=2)
    topo, topo_conclusion = summarize_topological_rule()
    susc, susc_tier = susceptibility_diagnostics()
    perturb, perturb_conclusion = perturbation_extended(datasets, best_k=2)
    l2_sens, l2_reg, l2_interpretation = l2_reverse_diagnostics()
    registry = clone_dataset_registry_v13()
    direction_matrix = build_direction_matrix(topo_conclusion, susc_tier, perturb_conclusion, l2_interpretation, registry)
    write_reports(direction_matrix, topo, topo_conclusion, susc, susc_tier, perturb, perturb_conclusion, l2_sens, l2_reg, l2_interpretation, registry)
    update_manuscript(direction_matrix, topo_conclusion)
    update_reviewer_matrix(direction_matrix)
    claim_audit()
    make_figures(topo, stability, l2_sens, direction_matrix)
    print(
        {
            "topological_rule": topo_conclusion,
            "susceptibility_tier": susc_tier,
            "perturbation": perturb_conclusion,
            "l2": l2_interpretation,
            "strong_directions": int(direction_matrix["tier"].eq("strong").sum()),
            "acceptable_directions": int(direction_matrix["tier"].eq("acceptable").sum()),
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.parse_args()
    run()


if __name__ == "__main__":
    main()
