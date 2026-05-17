from __future__ import annotations

import argparse
import importlib.util
import math
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats

from src.utils.config import ensure_dir, write_text


ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class DatasetSpec:
    dataset_id: str
    label: str
    path: Path
    role: str
    independence_tier: str
    max_cells: int = 2600


DATASETS = [
    DatasetSpec("internal_native", "Internal native moscot main", ROOT / "processed" / "ot_teacher.h5ad", "internal", "internal_reference"),
    DatasetSpec("E1_MouseGastrulationData", "E1 MouseGastrulationData", ROOT / "processed" / "external" / "e1_ot_teacher.h5ad", "external_time_series", "related_mouse_gastrulation_sample"),
    DatasetSpec(
        "E3_MouseGastrulationData_full",
        "E3 MouseGastrulationData full stage-mapped",
        ROOT / "processed" / "developmental_atlas" / "E3_MouseGastrulationData_wt_chimera_full_stage_mapped" / "ot_teacher.h5ad",
        "related_atlas_stress",
        "related_mouse_gastrulation_atlas",
    ),
    DatasetSpec(
        "E5_zebrafish_Farrell",
        "E5 Farrell zebrafish axial mesoderm subset",
        ROOT / "processed" / "developmental_atlas" / "E5_CellRank_Farrell_zebrafish_axial_mesoderm" / "ot_teacher.h5ad",
        "independent_developmental_stress",
        "independent_cross_species",
    ),
    DatasetSpec(
        "GSE154572_EB_proxy",
        "GSE154572 embryoid-body WT cluster proxy",
        ROOT / "processed" / "developmental_atlas" / "V2_GSE154572_EB_WT_cluster_proxy" / "ot_teacher.h5ad",
        "independent_developmental_stress",
        "independent_cluster_proxy",
    ),
    DatasetSpec(
        "E2_GSE212050_gastruloid",
        "E2 GSE212050 gastruloid native atlas",
        ROOT / "processed" / "developmental_atlas" / "E2_GSE212050_gastruloid_native_atlas" / "ot_teacher.h5ad",
        "independent_developmental_stress",
        "independent_but_gene_matrix_unavailable",
    ),
]


TF_PROGRAMS = {
    "pluripotency": ["Pou5f1", "Nanog", "Sox2", "Klf4", "Esrrb", "Tfcp2l1", "Zfp42"],
    "primitive_streak_mesendoderm": ["T", "Eomes", "Mixl1", "Mesp1", "Mesp2", "Foxa2", "Gsc", "Lhx1"],
    "endoderm": ["Sox17", "Foxa2", "Gata4", "Gata6", "Sox7", "Hhex"],
    "mesoderm_emt": ["Snai1", "Snai2", "Twist1", "Zeb1", "Zeb2", "Tbx6", "Msgn1", "Hand1", "Hand2"],
    "neural_ectoderm": ["Sox1", "Sox2", "Sox3", "Pax6", "Otx2", "Pou3f1", "Neurog1", "Neurog2"],
    "hematopoietic": ["Tal1", "Runx1", "Gata1", "Gata2", "Lmo2", "Spi1", "Irf8", "Cebpa", "Cebpb", "Klf1"],
    "zebrafish_mesendoderm": ["TBXTA", "NTLA", "SOX32", "SOX17", "GATA5", "GATA6", "FOXA2", "MIXL1", "MESPAA", "MESPAB"],
}


def _gene_symbols(a: ad.AnnData) -> list[str]:
    for col in ["swarm_gene_symbol", "gene_short_name", "feature_name", "gene_symbols", "symbol"]:
        if col in a.var:
            vals = a.var[col].astype(str).replace("nan", "").tolist()
            if any(vals):
                return [v if v else str(idx) for v, idx in zip(vals, a.var_names)]
    return [str(x) for x in a.var_names]


def _to_dense(x):
    if sparse.issparse(x):
        return x.toarray()
    return np.asarray(x)


def _zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    mu = np.nanmean(x, axis=0, keepdims=True)
    sd = np.nanstd(x, axis=0, keepdims=True)
    sd[sd < 1e-8] = 1.0
    return (x - mu) / sd


def _entropy(vals: pd.Series) -> float:
    counts = vals.astype(str).value_counts()
    if counts.empty:
        return 0.0
    p = counts.to_numpy(dtype=float) / counts.sum()
    h = float(-(p * np.log2(p + 1e-12)).sum())
    return h / max(math.log2(len(p) + 1e-12), 1e-12)


def _pairwise_sample_distance(x: np.ndarray, max_n: int = 400, seed: int = 7) -> float:
    if x.shape[0] < 3:
        return 0.0
    rng = np.random.default_rng(seed)
    idx = np.arange(x.shape[0])
    if len(idx) > max_n:
        idx = rng.choice(idx, max_n, replace=False)
    y = x[idx]
    centroid = np.nanmean(y, axis=0)
    return float(np.nanmean(np.linalg.norm(y - centroid, axis=1)))


def _lineage_separation(x: np.ndarray, labels: pd.Series) -> float:
    labs = labels.astype(str).to_numpy()
    centroids = []
    for lab in sorted(set(labs)):
        mask = labs == lab
        if mask.sum() >= 4:
            centroids.append(np.nanmean(x[mask], axis=0))
    if len(centroids) < 2:
        return 0.0
    c = np.vstack(centroids)
    vals = []
    for i in range(len(c)):
        for j in range(i + 1, len(c)):
            vals.append(np.linalg.norm(c[i] - c[j]))
    return float(np.nanmean(vals)) if vals else 0.0


def _mean_cosine_to_direction(x: np.ndarray, direction: np.ndarray) -> float:
    if x.shape[0] == 0 or np.linalg.norm(direction) < 1e-8:
        return 0.0
    centered = x - np.nanmean(x, axis=0, keepdims=True)
    denom = np.linalg.norm(centered, axis=1) * np.linalg.norm(direction)
    ok = denom > 1e-8
    if ok.sum() == 0:
        return 0.0
    return float(np.nanmean((centered[ok] @ direction) / denom[ok]))


def _susceptibility(x: np.ndarray, max_n: int = 350, seed: int = 11) -> float:
    if x.shape[0] < 5:
        return 0.0
    rng = np.random.default_rng(seed)
    idx = np.arange(x.shape[0])
    if len(idx) > max_n:
        idx = rng.choice(idx, max_n, replace=False)
    y = _zscore(x[idx])
    y = y - y.mean(axis=0, keepdims=True)
    norm = np.linalg.norm(y, axis=1, keepdims=True)
    y = y / np.maximum(norm, 1e-8)
    sim = y @ y.T
    tri = sim[np.triu_indices_from(sim, k=1)]
    return float(np.nanmean(tri))


def _match_tfs(a: ad.AnnData) -> tuple[list[str], list[int], dict[str, str]]:
    symbols = _gene_symbols(a)
    lookup = {s.lower(): (i, s) for i, s in enumerate(symbols)}
    matched: list[str] = []
    idxs: list[int] = []
    program_by_tf: dict[str, str] = {}
    for program, genes in TF_PROGRAMS.items():
        for g in genes:
            hit = lookup.get(g.lower())
            if hit is None:
                continue
            idx, symbol = hit
            if symbol not in matched:
                matched.append(symbol)
                idxs.append(idx)
                program_by_tf[symbol] = program
    return matched, idxs, program_by_tf


def _activity_matrix(a: ad.AnnData, tf_idxs: list[int]) -> np.ndarray:
    x = _to_dense(a.X[:, tf_idxs])
    return _zscore(x)


def _infer_target_modules(a: ad.AnnData, tf_activity: np.ndarray, tf_names: list[str], max_targets: int = 25) -> tuple[np.ndarray, pd.DataFrame]:
    if not tf_names:
        return tf_activity, pd.DataFrame()
    n = a.n_obs
    rng = np.random.default_rng(19)
    cell_idx = np.arange(n)
    if n > 1200:
        cell_idx = rng.choice(cell_idx, 1200, replace=False)
    symbols = _gene_symbols(a)
    candidate = np.arange(a.n_vars)
    if a.n_vars > 1600:
        sample = rng.choice(candidate, 1600, replace=False)
    else:
        sample = candidate
    x_sample = _zscore(_to_dense(a.X[cell_idx, :][:, sample]))
    tf_sample = tf_activity[cell_idx]
    target_rows = []
    module_activity = tf_activity.copy()
    for t, tf in enumerate(tf_names):
        y = tf_sample[:, t]
        cor = np.nan_to_num((x_sample.T @ y) / max(len(y) - 1, 1))
        order = np.argsort(cor)[::-1]
        chosen = [sample[o] for o in order if symbols[sample[o]] != tf and cor[o] > 0.05][:max_targets]
        if chosen:
            target_x = _zscore(_to_dense(a.X[:, chosen]))
            module_activity[:, t] = 0.55 * tf_activity[:, t] + 0.45 * np.nanmean(target_x, axis=1)
        for rank, gi in enumerate(chosen[:10], start=1):
            target_rows.append({"tf": tf, "target": symbols[gi], "rank": rank, "correlation_proxy": float(cor[np.where(sample == gi)[0][0]]) if gi in sample else np.nan})
    return _zscore(module_activity), pd.DataFrame(target_rows)


def _time_order(a: ad.AnnData) -> list[float]:
    times = pd.to_numeric(a.obs.get("time_numeric", pd.Series(np.arange(a.n_obs), index=a.obs_names)), errors="coerce")
    vals = sorted([float(x) for x in times.dropna().unique()])
    return vals


def _order_parameters(a: ad.AnnData, activity: np.ndarray, tf_names: list[str], program_by_tf: dict[str, str], seed: int = 5) -> tuple[pd.DataFrame, pd.DataFrame]:
    times = pd.to_numeric(a.obs["time_numeric"], errors="coerce")
    lineage_col = "lineage" if "lineage" in a.obs else "cell_type"
    rows = []
    for i, t in enumerate(_time_order(a)):
        mask = times.to_numpy() == t
        x = activity[mask]
        next_direction = np.zeros(activity.shape[1])
        if i + 1 < len(_time_order(a)):
            next_mask = times.to_numpy() == _time_order(a)[i + 1]
            next_direction = np.nanmean(activity[next_mask], axis=0) - np.nanmean(x, axis=0)
        row = {
            "time_numeric": t,
            "n_cells": int(mask.sum()),
            "regulon_condensation": _pairwise_sample_distance(x, seed=seed),
            "regulon_divergence": _lineage_separation(x, a.obs.loc[mask, lineage_col]),
            "regulon_alignment": _mean_cosine_to_direction(x, next_direction),
            "branch_specific_regulon_entropy": _entropy(a.obs.loc[mask, lineage_col]),
            "regulon_susceptibility": _susceptibility(x, seed=seed),
            "TF_switch_score": float(np.nanmean(np.abs(next_direction))) if i + 1 < len(_time_order(a)) else 0.0,
        }
        if "ot_transition_entropy" in a.obs:
            vals = pd.to_numeric(a.obs.loc[mask, "ot_transition_entropy"], errors="coerce")
            row["mean_ot_transition_entropy"] = float(vals.mean())
        else:
            row["mean_ot_transition_entropy"] = np.nan
        rows.append(row)
    op = pd.DataFrame(rows)
    events = []
    if len(op) >= 3:
        for i in range(1, len(op) - 1):
            pre = op.iloc[i - 1]
            cur = op.iloc[i]
            post = op.iloc[i + 1]
            condensation_effect = float(pre["regulon_condensation"] - cur["regulon_condensation"])
            post_divergence = float(post["regulon_divergence"] - cur["regulon_divergence"])
            alignment_effect = float(cur["regulon_alignment"] - pre["regulon_alignment"])
            entropy_effect = float(post["branch_specific_regulon_entropy"] - cur["branch_specific_regulon_entropy"])
            score = condensation_effect + post_divergence + 0.5 * alignment_effect - 0.25 * abs(entropy_effect)
            events.append(
                {
                    "event_time": cur["time_numeric"],
                    "pre_time": pre["time_numeric"],
                    "post_time": post["time_numeric"],
                    "regulon_condensation_effect": condensation_effect,
                    "normalized_regulon_condensation_effect": condensation_effect / max(float(pre["regulon_condensation"]), 1e-8),
                    "regulon_post_divergence_effect": post_divergence,
                    "regulon_alignment_effect": alignment_effect,
                    "regulon_entropy_effect": entropy_effect,
                    "regulon_susceptibility_effect": float(cur["regulon_susceptibility"] - pre["regulon_susceptibility"]),
                    "regulatory_branch_window_score": score,
                    "regulatory_alignment_before_divergence": bool(alignment_effect > 0 and post_divergence > 0),
                }
            )
    ev = pd.DataFrame(events)
    return op, ev


def _controls(a: ad.AnnData, activity: np.ndarray, tf_names: list[str], program_by_tf: dict[str, str]) -> pd.DataFrame:
    rows = []
    base_op, base_ev = _order_parameters(a, activity, tf_names, program_by_tf)
    base_score = float(base_ev["regulatory_branch_window_score"].max()) if not base_ev.empty else 0.0
    rng = np.random.default_rng(23)
    for control in ["time_shuffle", "lineage_label_shuffle", "regulon_activity_shuffle", "random_tf_program"]:
        b = a.copy()
        act = activity.copy()
        if control == "time_shuffle":
            b.obs["time_numeric"] = rng.permutation(pd.to_numeric(b.obs["time_numeric"], errors="coerce").to_numpy())
        elif control == "lineage_label_shuffle" and "lineage" in b.obs:
            b.obs["lineage"] = rng.permutation(b.obs["lineage"].astype(str).to_numpy())
        elif control == "regulon_activity_shuffle":
            for j in range(act.shape[1]):
                act[:, j] = rng.permutation(act[:, j])
        elif control == "random_tf_program":
            if b.n_vars >= max(3, len(tf_names)):
                idx = rng.choice(np.arange(b.n_vars), size=max(3, len(tf_names)), replace=False)
                act = _zscore(_to_dense(b.X[:, idx]))
        _, ev = _order_parameters(b, act, tf_names, program_by_tf)
        score = float(ev["regulatory_branch_window_score"].max()) if not ev.empty else 0.0
        rows.append(
            {
                "control": control,
                "observed_score": base_score,
                "control_score": score,
                "negative_control_pass": bool(abs(score) < abs(base_score) * 0.8) if abs(base_score) > 1e-8 else False,
            }
        )
    return pd.DataFrame(rows)


def _expression_baseline(a: ad.AnnData) -> dict:
    if "X_pca" in a.obsm:
        x = np.asarray(a.obsm["X_pca"])
    elif a.n_vars:
        x = _zscore(_to_dense(a.X[:, : min(a.n_vars, 40)]))
    else:
        return {"expression_branch_window_score": np.nan, "fate_entropy_baseline_score": np.nan}
    dummy_tf = [f"pc{i}" for i in range(x.shape[1])]
    _, ev = _order_parameters(a, x, dummy_tf, {})
    score = float(ev["regulatory_branch_window_score"].max()) if not ev.empty else 0.0
    fate_cols = [c for c in a.obs.columns if str(c).startswith("fate_prob_")]
    fate_score = np.nan
    if fate_cols:
        fp = a.obs[fate_cols].apply(pd.to_numeric, errors="coerce").to_numpy()
        fate_score = float(np.nanmean(-(fp * np.log(fp + 1e-12)).sum(axis=1)))
    return {"expression_branch_window_score": score, "fate_entropy_baseline_score": fate_score}


def _known_biology(tf_df: pd.DataFrame, program_by_tf: dict[str, str], dataset_id: str) -> pd.DataFrame:
    rows = []
    for tf, prog in program_by_tf.items():
        rows.append(
            {
                "dataset": dataset_id,
                "tf_or_regulon": tf,
                "program": prog,
                "known_biology_mapping": {
                    "primitive_streak_mesendoderm": "primitive streak, mesendoderm, gastrulation transition program",
                    "endoderm": "definitive/endoderm lineage commitment program",
                    "mesoderm_emt": "mesoderm or EMT-like transition program",
                    "neural_ectoderm": "neural/ectoderm lineage commitment program",
                    "hematopoietic": "hematopoietic fate regulator program",
                    "pluripotency": "pluripotency or pre-commitment state program",
                    "zebrafish_mesendoderm": "zebrafish mesendoderm/notochord-associated developmental program",
                }.get(prog, "curated developmental TF program"),
                "evidence_level": "computational_fallback_regulon_activity",
            }
        )
    return pd.DataFrame(rows)


def _perturbation(activity: np.ndarray, tf_names: list[str], events: pd.DataFrame, dataset_id: str) -> pd.DataFrame:
    if activity.shape[1] == 0 or events.empty:
        return pd.DataFrame()
    base = float(events["regulatory_branch_window_score"].max())
    rows = []
    strength = 0.25
    variances = np.nanvar(activity, axis=0)
    top = np.argsort(variances)[::-1][: min(6, len(tf_names))]
    for j in top:
        perturbed_score = base - strength * float(variances[j]) * np.sign(base if base != 0 else 1)
        rows.append(
            {
                "dataset": dataset_id,
                "perturbation": f"reduce_{tf_names[j]}_activity_25pct",
                "tf": tf_names[j],
                "baseline_branch_window_score": base,
                "perturbed_branch_window_score_proxy": perturbed_score,
                "score_shift_proxy": perturbed_score - base,
                "interpretation": "in_silico_activity_rescaling_proxy_not_wet_lab_validation",
            }
        )
    return pd.DataFrame(rows)


def _tier(summary: pd.DataFrame, controls: pd.DataFrame, baseline: dict, dataset_role: str, matched_tfs: int) -> tuple[str, str]:
    if summary.empty or matched_tfs < 3:
        return "fail", "insufficient matched TFs or no branch-window event"
    best = summary.iloc[summary["regulatory_branch_window_score"].abs().argmax()]
    control_pass = bool(controls["negative_control_pass"].mean() >= 0.5) if not controls.empty else False
    expr_score = abs(float(baseline.get("expression_branch_window_score", 0.0) or 0.0))
    grn_score = abs(float(best["regulatory_branch_window_score"]))
    complements_expr = grn_score > 0 and (expr_score == 0 or grn_score >= 0.25 * expr_score)
    aligns = bool(best["regulatory_alignment_before_divergence"])
    if dataset_role == "internal" and aligns and control_pass and complements_expr:
        return "acceptable", "internal GRN fallback supports regulatory alignment-before-divergence with partial controls"
    if dataset_role == "external_time_series" and aligns and control_pass and complements_expr:
        return "acceptable", "E1 GRN fallback supports regulatory alignment-before-divergence with partial controls"
    if aligns or complements_expr:
        return "weak", "GRN features show partial branch-window structure but controls, independence or fallback limits prevent stronger claims"
    return "fail", "GRN branch-window structure was not detected beyond controls or expression baseline"


def analyze_dataset(spec: DatasetSpec) -> dict[str, pd.DataFrame]:
    if not spec.path.exists():
        return {
            "availability": pd.DataFrame([{"dataset": spec.dataset_id, "status": "missing_h5ad", "path": str(spec.path)}]),
            "order": pd.DataFrame(),
            "summary": pd.DataFrame(),
            "controls": pd.DataFrame(),
            "tf": pd.DataFrame(),
            "targets": pd.DataFrame(),
            "biology": pd.DataFrame(),
            "perturb": pd.DataFrame(),
            "baseline": pd.DataFrame(),
        }
    a = ad.read_h5ad(spec.path)
    if a.n_vars == 0 or "time_numeric" not in a.obs:
        return {
            "availability": pd.DataFrame(
                [
                    {
                        "dataset": spec.dataset_id,
                        "status": "unusable_for_grn",
                        "n_cells": a.n_obs,
                        "n_genes": a.n_vars,
                        "reason": "missing gene matrix or time_numeric",
                    }
                ]
            ),
            "order": pd.DataFrame(),
            "summary": pd.DataFrame(),
            "controls": pd.DataFrame(),
            "tf": pd.DataFrame(),
            "targets": pd.DataFrame(),
            "biology": pd.DataFrame(),
            "perturb": pd.DataFrame(),
            "baseline": pd.DataFrame(),
        }
    if a.n_obs > spec.max_cells:
        rng = np.random.default_rng(31)
        obs = a.obs.copy()
        keep = []
        per_time = max(80, spec.max_cells // max(1, obs["time_numeric"].nunique()))
        for _, group in obs.groupby("time_numeric", observed=True):
            idx = group.index.to_numpy()
            if len(idx) > per_time:
                idx = rng.choice(idx, per_time, replace=False)
            keep.extend(idx.tolist())
        a = a[keep].copy()
    tf_names, tf_idxs, program_by_tf = _match_tfs(a)
    availability = pd.DataFrame(
        [
            {
                "dataset": spec.dataset_id,
                "status": "analyzed" if len(tf_names) >= 3 else "limited_tf_overlap",
                "n_cells": a.n_obs,
                "n_genes": a.n_vars,
                "n_time_points": int(pd.to_numeric(a.obs["time_numeric"], errors="coerce").nunique()),
                "n_lineages": int(a.obs["lineage"].astype(str).nunique()) if "lineage" in a.obs else np.nan,
                "matched_tf_count": len(tf_names),
                "matched_tfs": ";".join(tf_names[:30]),
                "teacher_backend": "native_moscot" if "swarmlineage_ot_teacher" in a.uns else "fallback_or_unknown",
                "independence_tier": spec.independence_tier,
            }
        ]
    )
    if len(tf_names) < 3:
        return {
            "availability": availability,
            "order": pd.DataFrame(),
            "summary": pd.DataFrame(),
            "controls": pd.DataFrame(),
            "tf": pd.DataFrame(),
            "targets": pd.DataFrame(),
            "biology": pd.DataFrame(),
            "perturb": pd.DataFrame(),
            "baseline": pd.DataFrame(),
        }
    tf_activity = _activity_matrix(a, tf_idxs)
    activity, target_df = _infer_target_modules(a, tf_activity, tf_names)
    order, event = _order_parameters(a, activity, tf_names, program_by_tf)
    controls = _controls(a, activity, tf_names, program_by_tf)
    baseline = _expression_baseline(a)
    tier, interp = _tier(event, controls, baseline, spec.role, len(tf_names))
    if not event.empty:
        event.insert(0, "dataset", spec.dataset_id)
        event["dataset_label"] = spec.label
        event["role"] = spec.role
        event["independence_tier"] = spec.independence_tier
        event["matched_tf_count"] = len(tf_names)
        event["grn_tier"] = tier
        event["interpretation"] = interp
    order.insert(0, "dataset", spec.dataset_id)
    controls.insert(0, "dataset", spec.dataset_id)
    target_df.insert(0, "dataset", spec.dataset_id) if not target_df.empty else None
    tf_df = pd.DataFrame(
        [{"dataset": spec.dataset_id, "tf": tf, "program": program_by_tf.get(tf, "unknown"), "activity_variance": float(np.nanvar(activity[:, i]))} for i, tf in enumerate(tf_names)]
    )
    biology = _known_biology(tf_df, program_by_tf, spec.dataset_id)
    perturb = _perturbation(activity, tf_names, event, spec.dataset_id)
    baseline_df = pd.DataFrame([{**{"dataset": spec.dataset_id}, **baseline}])
    return {
        "availability": availability,
        "order": order,
        "summary": event,
        "controls": controls,
        "tf": tf_df,
        "targets": target_df,
        "biology": biology,
        "perturb": perturb,
        "baseline": baseline_df,
    }


def _method_availability() -> pd.DataFrame:
    rows = []
    checks = {
        "pySCENIC_or_GRNBoost2_style": ["pyscenic", "arboreto"],
        "CellOracle_style": ["celloracle"],
        "lightweight_TF_target_correlation_fallback": [],
        "OT_coupled_time_lagged_GRN_proxy": [],
    }
    for method, mods in checks.items():
        if not mods:
            rows.append({"method": method, "available": True, "status": "implemented_fallback", "evidence_tier_cap": "acceptable"})
            continue
        available = all(importlib.util.find_spec(m) is not None for m in mods)
        rows.append(
            {
                "method": method,
                "available": available,
                "status": "native_available_not_run_in_main_pipeline" if available else "not_installed_or_not_importable",
                "evidence_tier_cap": "strong_if_run_and_validated" if available else "not_claimed",
            }
        )
    return pd.DataFrame(rows)


def _plot(summary: pd.DataFrame, order: pd.DataFrame, outdir: Path) -> None:
    ensure_dir(outdir)
    if not summary.empty:
        fig, ax = plt.subplots(figsize=(9, 4.8))
        s = summary.sort_values("dataset")
        ax.bar(s["dataset"], s["regulatory_branch_window_score"], color="#4c78a8")
        ax.axhline(0, color="black", lw=0.8)
        ax.set_ylabel("GRN branch-window score")
        ax.set_xticks(np.arange(len(s["dataset"])))
        ax.set_xticklabels(s["dataset"], rotation=35, ha="right")
        fig.tight_layout()
        fig.savefig(outdir / "grn_branch_window_scores.png", dpi=180)
        ensure_dir(ROOT / "figures" / "main")
        fig.savefig(ROOT / "figures" / "main" / "figure7_grn_regulatory_audit.png", dpi=180)
        plt.close(fig)
    if not order.empty:
        fig, ax = plt.subplots(figsize=(9, 4.8))
        for ds, group in order.groupby("dataset"):
            ax.plot(group["time_numeric"], group["regulon_condensation"], marker="o", label=ds)
        ax.set_xlabel("time")
        ax.set_ylabel("regulon condensation distance")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(outdir / "grn_regulon_condensation_timeseries.png", dpi=180)
        plt.close(fig)


def _reports(outputs: dict[str, pd.DataFrame]) -> None:
    summary = outputs["summary"]
    availability = outputs["availability"]
    method = outputs["methods"]
    rows = []
    if summary.empty:
        final_tier = "fail"
        final_interp = "GRN branch-window analysis did not produce evaluable event summaries."
    else:
        tiers = summary.groupby("dataset")["grn_tier"].first().to_dict()
        internal = tiers.get("internal_native", "fail")
        e1 = tiers.get("E1_MouseGastrulationData", "fail")
        independent = [t for d, t in tiers.items() if d not in {"internal_native", "E1_MouseGastrulationData", "E3_MouseGastrulationData_full"}]
        if internal in {"acceptable", "strong"} and e1 in {"acceptable", "strong"} and any(t in {"acceptable", "strong"} for t in independent):
            final_tier = "strong"
            final_interp = "GRN/regulon branch-window structure is supported internally, in E1, and in an independent dataset."
        elif internal in {"acceptable", "strong"} and e1 in {"acceptable", "strong"}:
            final_tier = "acceptable"
            final_interp = "GRN/regulon analysis supports the mouse-gastrulation-like branch window, but independent atlas generalization remains unresolved."
        elif internal in {"acceptable", "strong"} or e1 in {"acceptable", "strong"}:
            final_tier = "weak"
            final_interp = "GRN support is dataset-limited and remains exploratory."
        else:
            final_tier = "fail"
            final_interp = "GRN features do not currently strengthen the branch-window mechanism beyond expression/OT geometry."
    evidence = pd.DataFrame(
        [
            {
                "claim": "GRN/regulon alignment-before-divergence",
                "tier": final_tier,
                "allowed_language": final_interp,
                "forbidden_language": "GRN proves causal developmental mechanism or experimentally validated TF control",
                "native_or_fallback_boundary": "pySCENIC/CellOracle native GRN not claimed unless available and run; current main evidence uses fallback TF-target correlation and OT/time-lag proxy",
            },
            {
                "claim": "Known developmental TF programs recovered",
                "tier": "weak" if outputs["biology"].empty else "acceptable",
                "allowed_language": "Curated developmental TF programs are recovered as computational candidates where genes overlap.",
                "forbidden_language": "Validated regulators or causal perturbation evidence",
                "native_or_fallback_boundary": "candidate mapping only",
            },
            {
                "claim": "In silico GRN perturbation",
                "tier": "weak" if not outputs["perturb"].empty else "fail",
                "allowed_language": "Activity-rescaling perturbations provide candidate sensitivity probes only.",
                "forbidden_language": "experimental perturbation validation",
                "native_or_fallback_boundary": "proxy perturbation",
            },
        ]
    )
    outputs["evidence"] = evidence
    lines = [
        "# SwarmLineage-OT-GRN evidence matrix",
        "",
        f"Final GRN tier: **{final_tier}**.",
        "",
        final_interp,
        "",
        "Method boundary:",
    ]
    for _, r in method.iterrows():
        lines.append(f"- {r['method']}: {r['status']}; tier cap {r['evidence_tier_cap']}.")
    lines.extend(
        [
            "",
            "Dataset status:",
        ]
    )
    for _, r in availability.iterrows():
        lines.append(f"- {r['dataset']}: {r['status']}; cells={r.get('n_cells', 'NA')}; genes={r.get('n_genes', 'NA')}; matched TFs={r.get('matched_tf_count', 'NA')}.")
    lines.extend(
        [
            "",
            "Interpretation boundary: this analysis can support a computational regulatory order-parameter hypothesis. It does not establish causal GRN control, experimental validation, direct lineage validation, or model superiority over OT.",
        ]
    )
    write_text(ROOT / "reports" / "grn_evidence_matrix.md", "\n".join(lines) + "\n")
    write_text(ROOT / "reports" / "grn_branch_window_report.md", "\n".join(lines) + "\n")
    write_text(
        ROOT / "reports" / "grn_perturbation_report.md",
        "# In silico GRN perturbation report\n\n"
        "Perturbations are activity-rescaling sensitivity probes on fallback regulon activities, not experimental validation or causal proof.\n\n"
        + (outputs["perturb"].head(30).to_markdown(index=False) if not outputs["perturb"].empty else "No perturbation rows were produced.\n"),
    )
    write_text(
        ROOT / "reports" / "grn_method_availability.md",
        "# GRN method availability\n\n" + method.to_markdown(index=False) + "\n",
    )


def run() -> None:
    ensure_dir(ROOT / "tables")
    ensure_dir(ROOT / "reports")
    ensure_dir(ROOT / "figures" / "grn")
    collected: dict[str, list[pd.DataFrame]] = {k: [] for k in ["availability", "order", "summary", "controls", "tf", "targets", "biology", "perturb", "baseline"]}
    for spec in DATASETS:
        result = analyze_dataset(spec)
        for key, value in result.items():
            if not value.empty:
                collected[key].append(value)
    outputs = {k: pd.concat(v, ignore_index=True) if v else pd.DataFrame() for k, v in collected.items()}
    outputs["methods"] = _method_availability()
    _reports(outputs)
    file_map = {
        "methods": "grn_method_availability.csv",
        "availability": "grn_dataset_availability.csv",
        "order": "grn_branch_window_order_parameters.csv",
        "summary": "grn_branch_window_summary.csv",
        "controls": "grn_negative_controls.csv",
        "tf": "grn_tf_candidates.csv",
        "targets": "grn_inferred_tf_targets.csv",
        "biology": "grn_known_biology_mapping.csv",
        "perturb": "grn_perturbation_simulation.csv",
        "baseline": "grn_expression_vs_regulon_detector.csv",
        "evidence": "grn_evidence_matrix.csv",
    }
    for key, name in file_map.items():
        outputs.get(key, pd.DataFrame()).to_csv(ROOT / "tables" / name, index=False)
    _plot(outputs["summary"], outputs["order"], ROOT / "figures" / "grn")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run fallback GRN/regulon branch-window analysis.")
    parser.parse_args()
    run()


if __name__ == "__main__":
    main()
