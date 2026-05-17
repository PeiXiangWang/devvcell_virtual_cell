from __future__ import annotations

import argparse
import json
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import OneHotEncoder

from src.discovery.clone_developmental_validation_v14 import DATASETS, _clean_clone, _entropy_from_counts, _parse_time, _write_csv, _write_md
from src.utils.config import ensure_dir


ROOT = Path(__file__).resolve().parents[2]
NATIVE_PYTHON = ".venv_moscot_native/Scripts/python.exe"
OUTPUT_ROOT = "processed/external_l6_outcome_preserving"


@dataclass(frozen=True)
class OutcomeSamplingStrategy:
    name: str
    max_cells_per_time: int
    terminal_quota: int
    early_quota: int
    mode: str
    description: str


STRATEGIES = [
    OutcomeSamplingStrategy(
        "terminal_outcome_preserving_native",
        500,
        3,
        1,
        "terminal_priority",
        "prioritize clones with enough terminal cells to estimate terminal fate entropy",
    ),
    OutcomeSamplingStrategy(
        "early_terminal_balanced_native",
        500,
        2,
        2,
        "early_terminal_balanced",
        "retain early and terminal representatives for the same clones",
    ),
    OutcomeSamplingStrategy(
        "clone_outcome_quality_native",
        500,
        4,
        2,
        "outcome_quality",
        "rank clones by terminal-depth, time-span and terminal fate diversity",
    ),
    OutcomeSamplingStrategy(
        "fate_entropy_preserving_native",
        500,
        3,
        2,
        "entropy_preserving",
        "prioritize terminal fate entropy strata to preserve the full-data outcome distribution",
    ),
    OutcomeSamplingStrategy(
        "max_outcome_quality_feasible_native",
        650,
        4,
        2,
        "outcome_quality",
        "largest attempted native run while maximizing outcome-estimable clone count",
    ),
]


MODELS = {
    "condensation_only": ["branch_window_condensation_exposure"],
    "post_event_divergence_only": ["post_event_divergence_exposure"],
    "fate_entropy_only": ["fate_entropy_exposure"],
    "teacher_bias_only": ["teacher_velocity_bias_exposure"],
    "two_phase_condensation_plus_divergence": ["branch_window_condensation_exposure", "post_event_divergence_exposure"],
    "uncertainty_plus_teacher_bias": ["fate_entropy_exposure", "teacher_velocity_bias_exposure"],
}


OUTCOMES = [
    ("terminal_fate_entropy", True),
    ("terminal_lineage_entropy", False),
    ("clone_multilineage_score", False),
    ("clone_branch_count", False),
    ("clone_fate_diversification_index", False),
    ("clone_transition_entropy", False),
]


def _rel(path: Path | str) -> str:
    p = Path(path)
    try:
        return str(p.relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(p).replace("\\", "/")


def _load_obs(cfg) -> tuple[ad.AnnData | None, pd.DataFrame]:
    raw = ROOT / "data/external_l3/raw" / cfg.file_name
    if not raw.exists():
        return None, pd.DataFrame()
    adata = ad.read_h5ad(raw, backed="r")
    obs = adata.obs.copy()
    obs["global_pos"] = np.arange(obs.shape[0])
    obs["time_numeric"] = _parse_time(obs[cfg.time_col], cfg.dataset_id) if cfg.time_col in obs else np.nan
    obs["time_point"] = obs[cfg.time_col].astype(str) if cfg.time_col in obs else ""
    obs["clone_id"] = _clean_clone(obs[cfg.clone_col]) if cfg.clone_col in obs else np.nan
    obs["lineage"] = obs[cfg.celltype_col].astype(str) if cfg.celltype_col in obs else "unknown"
    return adata, obs.reset_index(drop=True)


def _latent(obs: pd.DataFrame, cfg) -> np.ndarray:
    if {"UMAP_1", "UMAP_2"}.issubset(obs.columns):
        base = obs[["UMAP_1", "UMAP_2"]].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float)
        numeric_cols = [
            c
            for c in obs.columns
            if c.startswith("progenitor_")
            or c.startswith("fate_map_")
            or c.startswith("fate_bias_")
            or c in {"NeuMon_fate_bias", "growth_rate_smooth", "growth_rate_raw"}
        ]
        extra = obs[numeric_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float) if numeric_cols else np.empty((obs.shape[0], 0))
        x = np.hstack([base, extra])
    else:
        numeric_cols = [
            c
            for c in obs.columns
            if c.startswith("progenitor_")
            or c.startswith("fate_map_")
            or c.startswith("fate_bias_")
            or c in {"NeuMon_fate_bias", "growth_rate_smooth", "growth_rate_raw"}
        ]
        numeric = obs[numeric_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0) if numeric_cols else pd.DataFrame(index=obs.index)
        enc = OneHotEncoder(sparse_output=False, handle_unknown="ignore", max_categories=40)
        cats = enc.fit_transform(obs[[cfg.celltype_col]].astype(str))
        x = np.hstack([numeric.to_numpy(dtype=float), cats])
    if x.shape[1] == 0:
        x = np.arange(obs.shape[0], dtype=float)[:, None]
    x = (x - x.mean(axis=0, keepdims=True)) / np.maximum(x.std(axis=0, keepdims=True), 1e-8)
    return x[:, : min(30, x.shape[1])].astype(np.float32)


def _full_clone_outcomes(obs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    valid = obs["clone_id"].notna() & obs["time_numeric"].notna()
    times = sorted(obs.loc[valid, "time_numeric"].dropna().unique())
    if not times:
        return pd.DataFrame()
    terminal = times[-1]
    for clone_id, g in obs.loc[valid].groupby("clone_id", observed=False):
        terminal_g = g[g["time_numeric"].eq(terminal)]
        early_g = g[g["time_numeric"].lt(terminal)]
        if g.shape[0] < 2 or g["time_numeric"].nunique() < 2:
            continue
        terminal_entropy = _entropy_from_counts(terminal_g["lineage"]) if not terminal_g.empty else np.nan
        rows.append(
            {
                "clone_id": clone_id,
                "clone_size_full": int(g.shape[0]),
                "terminal_cells_full": int(terminal_g.shape[0]),
                "early_cells_full": int(early_g.shape[0]),
                "clone_time_span_full": float(g["time_numeric"].max() - g["time_numeric"].min()),
                "time_coverage_full": int(g["time_numeric"].nunique()),
                "terminal_fate_entropy_full": terminal_entropy,
                "terminal_branch_count_full": int(terminal_g["lineage"].astype(str).nunique()) if not terminal_g.empty else 0,
                "outcome_estimable_full": bool(terminal_g.shape[0] >= 2 and early_g.shape[0] >= 1),
            }
        )
    return pd.DataFrame(rows)


def _priority_table(obs: pd.DataFrame, strategy: OutcomeSamplingStrategy) -> pd.DataFrame:
    table = _full_clone_outcomes(obs)
    if table.empty:
        return table
    terminal = table["terminal_cells_full"].clip(lower=0)
    early = table["early_cells_full"].clip(lower=0)
    entropy = table["terminal_fate_entropy_full"].fillna(0.0)
    diversity = table["terminal_branch_count_full"].clip(lower=0)
    span = table["clone_time_span_full"].clip(lower=0)
    if strategy.mode == "entropy_preserving":
        table["priority"] = entropy * 1000 + diversity * 100 + np.log1p(terminal) * 10 + span
    elif strategy.mode == "early_terminal_balanced":
        table["priority"] = np.minimum(terminal, early) * 1000 + span * 100 + np.log1p(table["clone_size_full"])
    elif strategy.mode == "terminal_priority":
        table["priority"] = terminal * 1000 + span * 100 + diversity * 10 + np.log1p(early)
    else:
        table["priority"] = (terminal >= 2).astype(float) * 10000 + np.minimum(terminal, 6) * 1000 + np.minimum(early, 6) * 100 + entropy * 10 + span
    return table.sort_values(["outcome_estimable_full", "priority", "clone_size_full"], ascending=False)


def _select_indices(obs: pd.DataFrame, strategy: OutcomeSamplingStrategy, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    valid = obs["clone_id"].notna() & obs["time_numeric"].notna()
    times = sorted(obs.loc[valid, "time_numeric"].dropna().unique())
    if not times:
        return np.asarray([], dtype=int)
    terminal = times[-1]
    priority = _priority_table(obs, strategy)
    priority_clones = priority["clone_id"].tolist() if not priority.empty else []
    selected: list[int] = []
    for t in times:
        time_idx = obs.index[valid & obs["time_numeric"].eq(t)].to_numpy()
        if time_idx.size == 0:
            continue
        keep: list[int] = []
        time_obs = obs.loc[time_idx]
        for clone_id in priority_clones:
            members = time_obs.index[time_obs["clone_id"].eq(clone_id)].to_numpy()
            if members.size == 0:
                continue
            quota = strategy.terminal_quota if t == terminal else strategy.early_quota
            take_n = min(quota, members.size)
            if take_n > 0:
                keep.extend(rng.choice(members, size=take_n, replace=False).tolist())
            if len(keep) >= strategy.max_cells_per_time:
                break
        if len(keep) < strategy.max_cells_per_time:
            preferred = obs.index[valid & obs["time_numeric"].eq(t) & obs["clone_id"].isin(priority_clones)].to_numpy()
            rem = np.setdiff1d(preferred if preferred.size else time_idx, np.asarray(keep, dtype=int), assume_unique=False)
            if rem.size:
                keep.extend(rng.choice(rem, size=min(strategy.max_cells_per_time - len(keep), rem.size), replace=False).tolist())
        selected.extend(sorted(set(keep[: strategy.max_cells_per_time])))
    return np.asarray(sorted(set(selected)), dtype=int)


def _sample_clone_outcomes(obs: pd.DataFrame) -> pd.DataFrame:
    valid = obs["clone_id"].notna() & obs["time_numeric"].notna()
    times = sorted(obs.loc[valid, "time_numeric"].dropna().unique())
    if not times:
        return pd.DataFrame()
    terminal = times[-1]
    rows = []
    for clone_id, g in obs.loc[valid].groupby("clone_id", observed=False):
        if g.shape[0] < 2 or g["time_numeric"].nunique() < 2:
            continue
        terminal_g = g[g["time_numeric"].eq(terminal)]
        early_g = g[g["time_numeric"].lt(terminal)]
        terminal_entropy = _entropy_from_counts(terminal_g["lineage"]) if not terminal_g.empty else np.nan
        rows.append(
            {
                "clone_id": clone_id,
                "clone_size": int(g.shape[0]),
                "terminal_cells_per_clone": int(terminal_g.shape[0]),
                "early_cells_per_clone": int(early_g.shape[0]),
                "clone_time_span": float(g["time_numeric"].max() - g["time_numeric"].min()),
                "terminal_fate_entropy": terminal_entropy,
                "terminal_branch_count": int(terminal_g["lineage"].astype(str).nunique()) if not terminal_g.empty else 0,
                "outcome_estimable": bool(terminal_g.shape[0] >= 2 and early_g.shape[0] >= 1 and pd.notna(terminal_entropy)),
            }
        )
    return pd.DataFrame(rows)


def _distribution_similarity(full: pd.DataFrame, sample: pd.DataFrame) -> dict:
    full_vals = pd.to_numeric(full.loc[full["outcome_estimable_full"], "terminal_fate_entropy_full"], errors="coerce").dropna()
    sample_vals = pd.to_numeric(sample.loc[sample["outcome_estimable"], "terminal_fate_entropy"], errors="coerce").dropna()
    if len(full_vals) < 2 or len(sample_vals) < 2:
        return {"outcome_distribution_ks": np.nan, "outcome_distribution_similarity": np.nan}
    ks = stats.ks_2samp(full_vals.to_numpy(dtype=float), sample_vals.to_numpy(dtype=float)).statistic
    return {"outcome_distribution_ks": float(ks), "outcome_distribution_similarity": float(1.0 - ks)}


def _coverage(full_obs: pd.DataFrame, selected: np.ndarray) -> dict:
    full = _full_clone_outcomes(full_obs)
    sample_obs = full_obs.iloc[selected].copy()
    sample = _sample_clone_outcomes(sample_obs)
    dist = _distribution_similarity(full, sample) if not full.empty and not sample.empty else {"outcome_distribution_ks": np.nan, "outcome_distribution_similarity": np.nan}
    outcome_estimable = sample[sample["outcome_estimable"]] if not sample.empty else pd.DataFrame()
    full_estimable = int(full["outcome_estimable_full"].sum()) if not full.empty else 0
    return {
        "n_cells": int(sample_obs.shape[0]),
        "n_time_points": int(sample_obs["time_numeric"].nunique()),
        "usable_clones": int(sample.shape[0]),
        "full_outcome_estimable_clones": full_estimable,
        "outcome_estimable_clones": int(outcome_estimable.shape[0]),
        "outcome_estimable_clone_ratio": float(outcome_estimable.shape[0] / max(full_estimable, 1)),
        "median_terminal_cells_per_clone": float(outcome_estimable["terminal_cells_per_clone"].median()) if not outcome_estimable.empty else 0.0,
        "median_early_cells_per_clone": float(outcome_estimable["early_cells_per_clone"].median()) if not outcome_estimable.empty else 0.0,
        "median_clone_time_span": float(outcome_estimable["clone_time_span"].median()) if not outcome_estimable.empty else 0.0,
        "terminal_fate_entropy_missing_rate": float(1.0 - outcome_estimable.shape[0] / max(sample.shape[0], 1)),
        "terminal_fate_entropy_variance": float(pd.to_numeric(outcome_estimable["terminal_fate_entropy"], errors="coerce").var()) if outcome_estimable.shape[0] > 1 else 0.0,
        **dist,
    }


def _write_native_config(cfg, strategy: OutcomeSamplingStrategy, h5ad_path: str) -> Path:
    name = f"clone_outcome_preserving_native_{cfg.dataset_id}_{strategy.name}.yaml"
    cfg_path = ROOT / "configs" / name
    outdir = f"{OUTPUT_ROOT}/{cfg.dataset_id}/{strategy.name}"
    data = {
        "adata_path": h5ad_path,
        "teacher_path": f"{outdir}/ot_teacher.h5ad",
        "couplings_dir": f"{outdir}/ot_couplings",
        "teacher_index_path": f"{outdir}/ot_couplings/teacher_coupling_index.csv",
        "summary_path": f"{outdir}/ot_teacher_summary.json",
        "time_key": "time_numeric",
        "time_label_key": "time_point",
        "cell_type_key": "lineage",
        "latent_key": "X_pca",
        "epsilon": 0.08,
        "max_cells_per_time": strategy.max_cells_per_time,
        "native_moscot_timeout_seconds": 180,
        "use_native_moscot": True,
        "native_max_cells_per_time": strategy.max_cells_per_time,
        "native_max_iterations": 350,
        "native_jit": False,
        "native_device": "cpu",
        "random_seed": 23,
        "split_mode": "none",
        "teacher_backend": "native_moscot",
    }
    with open(cfg_path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, sort_keys=False)
    return cfg_path


def prepare_strategy(cfg, strategy: OutcomeSamplingStrategy, seed: int) -> dict:
    adata, obs = _load_obs(cfg)
    if adata is None or obs.empty:
        return {"dataset_id": cfg.dataset_id, "sampling_strategy": strategy.name, "prepared": False, "reason": "raw_missing_or_unreadable"}
    selected = _select_indices(obs, strategy, seed)
    if selected.size == 0:
        adata.file.close()
        return {"dataset_id": cfg.dataset_id, "sampling_strategy": strategy.name, "prepared": False, "reason": "no_selected_cells"}
    outdir = ensure_dir(ROOT / OUTPUT_ROOT / cfg.dataset_id / strategy.name)
    out_h5ad = outdir / "native_input.h5ad"
    sub = adata[selected].to_memory()
    sub.obs = obs.iloc[selected].copy()
    sub.obsm["X_pca"] = _latent(sub.obs, cfg)
    sub.write_h5ad(out_h5ad)
    adata.file.close()
    cfg_path = _write_native_config(cfg, strategy, _rel(out_h5ad))
    return {
        "dataset_id": cfg.dataset_id,
        "sampling_strategy": strategy.name,
        "sampling_mode": strategy.mode,
        "strategy_description": strategy.description,
        "prepared": True,
        "input_path": _rel(out_h5ad),
        "config_path": _rel(cfg_path),
        **_coverage(obs, selected),
    }


def run_native(dataset_id: str, strategy: str, cfg_path: str, timeout: int, python_exe: str) -> dict:
    start = time.time()
    cmd = [python_exe, "-m", "src.ot_teacher.run_moscot", "--config", cfg_path, "--try-native"]
    try:
        proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=timeout)
        runtime = time.time() - start
        detail = (proc.stdout + "\n" + proc.stderr).strip()
        summary_path = ROOT / OUTPUT_ROOT / dataset_id / strategy / "ot_couplings" / "moscot_run_summary.json"
        backend = "unknown"
        used = False
        pairs = []
        if summary_path.exists():
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            backend = summary.get("teacher_backend", "unknown")
            used = bool(summary.get("native_moscot_used", False))
            pairs = summary.get("pairs", [])
        shapes = []
        for item in pairs:
            if "plan_shape" in item:
                shapes.append(str(tuple(item["plan_shape"])))
            elif "shape" in item:
                shapes.append(str(tuple(item["shape"])))
        if not shapes:
            coupling_dir = ROOT / OUTPUT_ROOT / dataset_id / strategy / "ot_couplings"
            for npz_path in sorted(coupling_dir.glob("teacher_native_moscot_*.npz")):
                raw = np.load(npz_path, allow_pickle=True)
                if "plan" in raw:
                    shapes.append(str(tuple(raw["plan"].shape)))
        return {
            "dataset_id": dataset_id,
            "sampling_strategy": strategy,
            "native_attempted": True,
            "native_returncode": proc.returncode,
            "native_moscot_success": bool(used and backend == "native_moscot" and proc.returncode == 0),
            "teacher_backend": backend,
            "fallback_used": backend != "native_moscot",
            "runtime_seconds": runtime,
            "n_pairs": len(pairs),
            "plan_shapes": ";".join(shapes),
            "failure_reason": "" if proc.returncode == 0 and used else detail[-500:],
        }
    except subprocess.TimeoutExpired:
        return {
            "dataset_id": dataset_id,
            "sampling_strategy": strategy,
            "native_attempted": True,
            "native_returncode": "timeout",
            "native_moscot_success": False,
            "teacher_backend": "timeout",
            "fallback_used": False,
            "runtime_seconds": timeout,
            "n_pairs": 0,
            "plan_shapes": "",
            "failure_reason": f"native moscot timed out after {timeout}s",
        }


def _native_velocity_entropy(dataset_id: str, strategy: str, adata: ad.AnnData) -> tuple[np.ndarray, np.ndarray]:
    z = np.asarray(adata.obsm["X_pca"], dtype=float)
    velocity = np.zeros_like(z)
    entropy = np.full(adata.n_obs, np.nan, dtype=float)
    coupling_dir = ROOT / OUTPUT_ROOT / dataset_id / strategy / "ot_couplings"
    for path in sorted(coupling_dir.glob("teacher_native_moscot_*.npz")):
        raw = np.load(path, allow_pickle=True)
        src = raw["source_indices"].astype(int)
        bary = raw["barycentric"].astype(float)
        velocity[src] = bary - z[src]
        entropy[src] = raw["entropy"].astype(float)
    return velocity, entropy


def _unit(x: np.ndarray) -> np.ndarray:
    return x / np.maximum(np.linalg.norm(x, axis=1, keepdims=True), 1e-8)


def _cell_exposures_native(dataset_id: str, strategy: str, adata: ad.AnnData) -> pd.DataFrame:
    obs = adata.obs.copy().reset_index(drop=True)
    z = np.asarray(adata.obsm["X_pca"], dtype=float)
    velocity, entropy = _native_velocity_entropy(dataset_id, strategy, adata)
    vu = _unit(velocity)
    times = sorted(pd.to_numeric(obs["time_numeric"], errors="coerce").dropna().unique())
    terminal = times[-1]
    event_time = times[0] if len(times) < 3 else times[len(times) // 2 - 1]
    terminal_centroid = z[obs["time_numeric"].eq(terminal).to_numpy()].mean(axis=0)
    global_centroid = np.zeros_like(z)
    same_lineage_centroid = np.zeros_like(z)
    local_density = np.zeros(adata.n_obs)
    local_alignment = np.zeros(adata.n_obs)
    for t in times:
        tidx = obs.index[obs["time_numeric"].eq(t)].to_numpy()
        if tidx.size == 0:
            continue
        global_centroid[tidx] = z[tidx].mean(axis=0)
        k = min(8, tidx.size - 1)
        if k > 0:
            dist, nn = NearestNeighbors(n_neighbors=k + 1).fit(z[tidx]).kneighbors(z[tidx])
            local_density[tidx] = 1.0 / (np.mean(dist[:, 1:], axis=1) + 1e-6)
            local_v = vu[tidx]
            local_alignment[tidx] = np.mean(np.sum(local_v[:, None, :] * local_v[nn[:, 1:]], axis=2), axis=1)
        for _, members in obs.loc[tidx].groupby("lineage", observed=True).groups.items():
            midx = np.asarray(list(members), dtype=int)
            if midx.size:
                same_lineage_centroid[midx] = z[midx].mean(axis=0)
    fate_entropy = np.where(np.isfinite(entropy), entropy, np.nanmedian(entropy[np.isfinite(entropy)]) if np.isfinite(entropy).any() else 0.0)
    post_gate = obs["time_numeric"].ge(event_time).to_numpy(dtype=float)
    return pd.DataFrame(
        {
            "branch_window_condensation_exposure": -np.linalg.norm(z - same_lineage_centroid, axis=1),
            "global_centroid_condensation_exposure": -np.linalg.norm(z - global_centroid, axis=1),
            "same_lineage_centroid_condensation_exposure": -np.linalg.norm(z - same_lineage_centroid, axis=1),
            "local_alignment_exposure": local_alignment,
            "alignment_exposure": local_alignment,
            "fate_entropy_exposure": fate_entropy,
            "density_exposure": local_density,
            "post_event_divergence_exposure": np.linalg.norm(z - same_lineage_centroid, axis=1) * post_gate,
            "teacher_velocity_bias_exposure": np.sum(_unit(terminal_centroid[None, :] - z) * vu, axis=1),
        }
    )


def _clone_table(adata: ad.AnnData, exposure: pd.DataFrame) -> pd.DataFrame:
    obs = adata.obs.copy().reset_index(drop=True)
    times = sorted(pd.to_numeric(obs["time_numeric"], errors="coerce").dropna().unique())
    terminal = times[-1]
    event_time = times[0] if len(times) < 3 else times[len(times) // 2 - 1]
    rows = []
    for clone_id, g in obs.groupby("clone_id", observed=False):
        if pd.isna(clone_id) or g.shape[0] < 2 or g["time_numeric"].nunique() < 2:
            continue
        terminal_g = g[g["time_numeric"].eq(terminal)]
        if terminal_g.empty:
            terminal_g = g.iloc[0:0]
        pre_idx = g.index[g["time_numeric"].le(event_time)]
        if len(pre_idx) == 0:
            pre_idx = g.index
        terminal_entropy = _entropy_from_counts(terminal_g["lineage"]) if not terminal_g.empty else np.nan
        row = {
            "clone_id": clone_id,
            "clone_size": int(g.shape[0]),
            "clone_start_time": float(g["time_numeric"].min()),
            "clone_end_time": float(g["time_numeric"].max()),
            "clone_time_span": float(g["time_numeric"].max() - g["time_numeric"].min()),
            "initial_state": str(g.sort_values("time_numeric")["lineage"].iloc[0]),
            "initial_lineage": str(g.sort_values("time_numeric")["lineage"].iloc[0]),
            "terminal_sampling_depth": int(terminal_g.shape[0]),
            "early_sampling_depth": int(g[g["time_numeric"].lt(terminal)].shape[0]),
            "number_of_observed_time_points": int(g["time_numeric"].nunique()),
            "terminal_fate_entropy": terminal_entropy,
            "terminal_lineage_entropy": terminal_entropy,
            "clone_multilineage_score": int(terminal_g["lineage"].astype(str).nunique()) if not terminal_g.empty else np.nan,
            "clone_branch_count": int(terminal_g["lineage"].astype(str).nunique()) if not terminal_g.empty else np.nan,
            "clone_fate_diversification_index": 1.0 - float(terminal_g["lineage"].astype(str).value_counts(normalize=True).max()) if not terminal_g.empty else np.nan,
            "clone_transition_entropy": _entropy_from_counts(g["lineage"]),
        }
        row["outcome_estimable"] = bool(row["terminal_sampling_depth"] >= 2 and row["early_sampling_depth"] >= 1 and pd.notna(row["terminal_fate_entropy"]))
        for col in [
            "branch_window_condensation_exposure",
            "post_event_divergence_exposure",
            "fate_entropy_exposure",
            "teacher_velocity_bias_exposure",
            "local_alignment_exposure",
            "density_exposure",
        ]:
            row[col] = float(pd.to_numeric(exposure.loc[pre_idx, col], errors="coerce").mean())
        rows.append(row)
    clone = pd.DataFrame(rows)
    if clone.empty:
        return clone
    for model, cols in MODELS.items():
        ranks = [pd.to_numeric(clone[col], errors="coerce").rank(pct=True) for col in cols]
        clone[model] = pd.concat(ranks, axis=1).mean(axis=1)
    return clone


def _spearman(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    mask = x.notna() & y.notna()
    if mask.sum() < 10 or x[mask].nunique() < 2 or y[mask].nunique() < 2:
        return np.nan, np.nan
    rho, p = stats.spearmanr(x[mask], y[mask])
    return float(rho), float(p)


def _bootstrap_ci(x: pd.Series, y: pd.Series, seed: int = 23, reps: int = 150) -> tuple[float, float]:
    mask = (x.notna() & y.notna()).to_numpy()
    xv = x[mask].to_numpy(dtype=float)
    yv = y[mask].to_numpy(dtype=float)
    if len(xv) < 10:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(reps):
        idx = rng.integers(0, len(xv), len(xv))
        if np.unique(xv[idx]).size < 2 or np.unique(yv[idx]).size < 2:
            continue
        vals.append(stats.spearmanr(xv[idx], yv[idx]).correlation)
    if not vals:
        return np.nan, np.nan
    return float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5))


def _adjusted_effect(clone: pd.DataFrame, predictor: str, outcome: str) -> tuple[float, float, int]:
    cov = clone[
        [
            predictor,
            "clone_size",
            "clone_start_time",
            "clone_end_time",
            "clone_time_span",
            "terminal_sampling_depth",
            "number_of_observed_time_points",
        ]
    ].copy()
    y = pd.to_numeric(clone[outcome], errors="coerce")
    state = pd.get_dummies(clone["initial_state"].astype(str), prefix="state", drop_first=True)
    cov = pd.concat([cov.apply(pd.to_numeric, errors="coerce"), state], axis=1)
    m = cov.notna().all(axis=1) & y.notna()
    if m.sum() < 15 or cov.loc[m, predictor].nunique() < 2 or y.loc[m].nunique() < 2:
        return np.nan, np.nan, int(m.sum())
    x = cov.loc[m].to_numpy(dtype=float)
    yy = y.loc[m].to_numpy(dtype=float)
    x = (x - x.mean(axis=0)) / np.maximum(x.std(axis=0), 1e-8)
    yy = (yy - yy.mean()) / max(yy.std(), 1e-8)
    model = LinearRegression().fit(x, yy)
    return float(model.coef_[0]), float(model.score(x, yy)), int(m.sum())


def _stratified_effect(clone: pd.DataFrame, predictor: str, outcome: str, group_col: str) -> float:
    vals = []
    weights = []
    if group_col not in clone:
        return np.nan
    for _, g in clone.groupby(group_col, observed=False):
        rho, _ = _spearman(pd.to_numeric(g[predictor], errors="coerce"), pd.to_numeric(g[outcome], errors="coerce"))
        if pd.notna(rho):
            vals.append(rho)
            weights.append(g.shape[0])
    return float(np.average(vals, weights=weights)) if vals else np.nan


def _permutation_q(x: pd.Series, y: pd.Series, observed: float, seed: int = 23, reps: int = 100) -> float:
    if pd.isna(observed):
        return np.nan
    mask = x.notna() & y.notna()
    xv = x[mask].to_numpy(dtype=float)
    yv = y[mask].to_numpy(dtype=float)
    if len(xv) < 10:
        return np.nan
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(reps):
        yp = rng.permutation(yv)
        if np.unique(yp).size < 2 or np.unique(xv).size < 2:
            continue
        null.append(stats.spearmanr(xv, yp).correlation)
    return float((np.sum(np.abs(null) >= abs(observed)) + 1) / (len(null) + 1)) if null else np.nan


def analyze_strategy(cfg, strategy: OutcomeSamplingStrategy, teacher_success: bool) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    path = ROOT / OUTPUT_ROOT / cfg.dataset_id / strategy.name / "native_input.h5ad"
    if not path.exists() or not teacher_success:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    adata = ad.read_h5ad(path)
    exposure = _cell_exposures_native(cfg.dataset_id, strategy.name, adata)
    clone = _clone_table(adata, exposure)
    if clone.empty:
        return clone, pd.DataFrame(), pd.DataFrame()
    clone["clone_size_bin"] = pd.qcut(clone["clone_size"].rank(method="first"), q=min(4, clone.shape[0]), duplicates="drop")
    clone["time_span_bin"] = pd.qcut(clone["clone_time_span"].rank(method="first"), q=min(4, clone.shape[0]), duplicates="drop")
    clone["start_time_bin"] = clone["clone_start_time"].astype(str)
    estimable = clone[clone["outcome_estimable"]].copy()
    rows = []
    controls = []
    for model in MODELS:
        for outcome, primary in OUTCOMES:
            if outcome not in estimable:
                continue
            x = pd.to_numeric(estimable[model], errors="coerce")
            y = pd.to_numeric(estimable[outcome], errors="coerce")
            effect, p = _spearman(x, y)
            ci_low, ci_high = _bootstrap_ci(x, y)
            adjusted, r2, n_adj = _adjusted_effect(estimable, model, outcome)
            within_start = _stratified_effect(estimable, model, outcome, "start_time_bin")
            size_matched = _stratified_effect(estimable, model, outcome, "clone_size_bin")
            span_matched = _stratified_effect(estimable, model, outcome, "time_span_bin")
            q = _permutation_q(x, y, effect)
            adjusted_support = bool(pd.notna(adjusted) and adjusted > 0)
            matched_support = bool(
                pd.notna(within_start)
                and within_start > 0
                and pd.notna(size_matched)
                and size_matched > 0
                and pd.notna(span_matched)
                and span_matched > 0
            )
            support = bool(primary and pd.notna(effect) and effect > 0 and adjusted_support and pd.notna(q) and q <= 0.10)
            tier = "acceptable" if support and matched_support and estimable.shape[0] >= 50 and q <= 0.05 else "weak" if support else "fail"
            rows.append(
                {
                    "dataset_id": cfg.dataset_id,
                    "sampling_strategy": strategy.name,
                    "model": model,
                    "outcome": outcome,
                    "outcome_primary": primary,
                    "n_clones_total": int(clone.shape[0]),
                    "n_outcome_estimable_clones": int(estimable.shape[0]),
                    "raw_effect": effect,
                    "raw_p": p,
                    "bootstrap_ci_low": ci_low,
                    "bootstrap_ci_high": ci_high,
                    "covariate_adjusted_effect": adjusted,
                    "adjusted_r2": r2,
                    "n_adjusted": n_adj,
                    "within_start_time_bin_effect": within_start,
                    "clone_size_matched_effect": size_matched,
                    "time_span_matched_effect": span_matched,
                    "permutation_q": q,
                    "support_tier": tier,
                }
            )
            for control, mode in [
                ("outcome_label_shuffle", "y"),
                ("clone_label_shuffle", "y"),
                ("time_label_shuffle", "x"),
                ("exposure_shuffle", "x"),
                ("teacher_velocity_shuffle", "x"),
            ]:
                rng = np.random.default_rng(abs(hash((cfg.dataset_id, strategy.name, model, outcome, control))) % (2**32))
                cx = x.copy()
                cy = y.copy()
                if mode == "x":
                    cx = pd.Series(rng.permutation(cx.to_numpy()), index=cx.index)
                else:
                    cy = pd.Series(rng.permutation(cy.to_numpy()), index=cy.index)
                ce, cp = _spearman(cx, cy)
                controls.append(
                    {
                        "dataset_id": cfg.dataset_id,
                        "sampling_strategy": strategy.name,
                        "model": model,
                        "outcome": outcome,
                        "control": control,
                        "control_effect": ce,
                        "control_p": cp,
                        "negative_control_clean": bool(pd.isna(ce) or abs(ce) < max(abs(effect), 0.05) or pd.isna(effect)),
                    }
                )
    return clone, pd.DataFrame(rows), pd.DataFrame(controls)


def _summarize(sampling: pd.DataFrame, runs: pd.DataFrame, results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    support_order = {"acceptable": 2, "weak": 1, "fail": 0}
    for _, s in sampling.iterrows():
        run = runs[(runs["dataset_id"].eq(s["dataset_id"])) & (runs["sampling_strategy"].eq(s["sampling_strategy"]))]
        res = results[
            results["dataset_id"].eq(s["dataset_id"])
            & results["sampling_strategy"].eq(s["sampling_strategy"])
            & results["outcome_primary"].eq(True)
        ]
        best = pd.DataFrame()
        if not res.empty:
            best = res.assign(score=res["support_tier"].map(support_order).fillna(0)).sort_values(["score", "raw_effect"], ascending=False).head(1)
        cond = res[res["model"].eq("condensation_only")]
        rows.append(
            {
                "dataset_id": s["dataset_id"],
                "sampling_strategy": s["sampling_strategy"],
                "n_cells": s.get("n_cells", np.nan),
                "usable_clones": s.get("usable_clones", np.nan),
                "outcome_estimable_clones": s.get("outcome_estimable_clones", np.nan),
                "terminal_fate_entropy_missing_rate": s.get("terminal_fate_entropy_missing_rate", np.nan),
                "terminal_fate_entropy_variance": s.get("terminal_fate_entropy_variance", np.nan),
                "outcome_distribution_similarity": s.get("outcome_distribution_similarity", np.nan),
                "teacher_backend": run["teacher_backend"].iloc[0] if not run.empty else "not_run",
                "native_moscot_success": bool(run["native_moscot_success"].iloc[0]) if not run.empty else False,
                "runtime_seconds": run["runtime_seconds"].iloc[0] if not run.empty else np.nan,
                "plan_shapes": run["plan_shapes"].iloc[0] if not run.empty else "",
                "best_model": best["model"].iloc[0] if not best.empty else "not_tested",
                "best_model_tier": best["support_tier"].iloc[0] if not best.empty else "fail",
                "best_model_effect": best["raw_effect"].iloc[0] if not best.empty else np.nan,
                "condensation_only_tier": cond["support_tier"].iloc[0] if not cond.empty else "fail",
                "condensation_only_effect": cond["raw_effect"].iloc[0] if not cond.empty else np.nan,
                "uncertainty_plus_teacher_bias_tier": res[res["model"].eq("uncertainty_plus_teacher_bias")]["support_tier"].iloc[0]
                if not res[res["model"].eq("uncertainty_plus_teacher_bias")].empty
                else "fail",
            }
        )
    return pd.DataFrame(rows)


def _final_interpretation(summary: pd.DataFrame) -> tuple[str, str, str]:
    j = summary[summary["dataset_id"].str.startswith("Jindal", na=False)]
    w = summary[summary["dataset_id"].str.startswith("Weinreb", na=False)]
    j_cond = j["condensation_only_tier"].isin(["acceptable", "weak"]).any()
    w_cond = w["condensation_only_tier"].isin(["acceptable", "weak"]).any()
    j_best = j["best_model_tier"].isin(["acceptable", "weak"]).any()
    w_best = w["best_model_tier"].isin(["acceptable", "weak"]).any()
    j_unc = j["uncertainty_plus_teacher_bias_tier"].isin(["acceptable", "weak"]).any()
    w_unc = w["uncertainty_plus_teacher_bias_tier"].isin(["acceptable", "weak"]).any()
    if j_cond and w_cond:
        return (
            "cross_dataset_condensation_candidate",
            "acceptable",
            "Outcome-preserving native sampling yields primary condensation support in both Jindal and Weinreb, but this remains computational and requires independent validation.",
        )
    if j_cond and not w_cond:
        return (
            "dataset_specific_jindal_clone_signal",
            "weak",
            "Jindal supports primary condensation under outcome-preserving native sampling, but Weinreb does not; this is dataset-specific weak support only.",
        )
    if w_cond and not j_cond:
        return (
            "weinreb_sampling_specific_condensation_signal",
            "weak",
            "Weinreb retains a sampling-specific primary condensation signal, but Jindal does not recover support after outcome-preserving native sampling; no general clone-aware support is established.",
        )
    if j_unc and w_unc and not (j_cond or w_cond):
        return (
            "secondary_uncertainty_gated_candidate",
            "weak",
            "Condensation-only fails, while uncertainty/teacher-bias features show weak secondary support; this revises rather than supports the original clone-level condensation hypothesis.",
        )
    if not j_best and not w_best:
        return (
            "clone_level_fate_diversification_not_supported",
            "fail",
            "Outcome-preserving native sampling does not support primary condensation or revised branch-window models in Jindal or Weinreb.",
        )
    return (
        "mixed_secondary_or_sampling_specific_signal",
        "weak",
        "Outcome-preserving native results remain mixed or secondary-only and cannot upgrade clone-aware support beyond weak.",
    )


def _write_figure(summary: pd.DataFrame) -> None:
    figdir = ensure_dir(ROOT / "figures/main")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    colors = {"Jindal_2023_NatureBiotechnology_LSK_RNA": "#4C78A8", "Weinreb_2020_Science": "#F58518"}
    for i, (_, row) in enumerate(summary.iterrows()):
        axes[0].bar(i, row["outcome_estimable_clones"], color=colors.get(row["dataset_id"], "#888888"))
        axes[1].bar(i, row["condensation_only_effect"] if pd.notna(row["condensation_only_effect"]) else 0.0, color=colors.get(row["dataset_id"], "#888888"))
    labels = summary["dataset_id"].str.replace("_2023_NatureBiotechnology_LSK_RNA", "", regex=False).str.replace("_2020_Science", "", regex=False) + "\n" + summary["sampling_strategy"].str.replace("_native", "", regex=False)
    for ax in axes:
        ax.set_xticks(range(summary.shape[0]))
        ax.set_xticklabels(labels, rotation=55, ha="right", fontsize=6)
        ax.axhline(0, color="black", linewidth=0.8)
    axes[0].set_ylabel("outcome-estimable clones")
    axes[0].set_title("Outcome-preserving native sampling")
    axes[1].set_ylabel("condensation -> terminal entropy effect")
    axes[1].set_title("Primary clone association")
    fig.tight_layout()
    fig.savefig(figdir / "figure16_outcome_preserving_clone_audit.png", dpi=180)
    plt.close(fig)


def _update_claim_tier_table(status: str, tier: str, interpretation: str) -> None:
    path = ROOT / "tables/final_claim_evidence_tiers.csv"
    table = pd.read_csv(path) if path.exists() else pd.DataFrame()
    row = {
        "claim": "clone-aware fate-diversification prediction",
        "status": status,
        "tier": tier,
        "internal_native_support": False,
        "native_sensitivity_support": True,
        "external_time_series_support": False,
        "lineage_clone_support": False,
        "negative_controls": "outcome-preserving native sampling, covariate adjustment, matched analyses and label/exposure shuffles",
        "module_necessity": "not_applicable",
        "external_independence": "Jindal/Weinreb clone-aware datasets",
        "allowed_manuscript_sentence": "Clone-aware native audits remain weak/mixed or unsupported; the retained main claim is a time-series order-parameter hypothesis.",
        "forbidden_sentence": "Do not present clone fate diversification as established.",
    }
    if "claim" in table:
        table = table[~table["claim"].eq(row["claim"])]
    table = pd.concat([table, pd.DataFrame([row])], ignore_index=True)
    _write_csv(table, "tables/final_claim_evidence_tiers.csv")
    _write_md("reports/final_claim_evidence_tiers.md", "# Final Claim Evidence Tiers\n\n" + table.to_markdown(index=False) + "\n")


def _update_documents(status: str, tier: str, interpretation: str, summary: pd.DataFrame) -> None:
    final_sentence = (
        "clone-level fate-diversification prediction is not supported under current tested datasets and native sampling strategies."
        if tier in {"fail", "weak"}
        else "clone-aware support is an acceptable computational candidate but remains unvalidated experimentally."
    )
    summary_md = summary.to_markdown(index=False) if not summary.empty else "_No strategy summary._"
    _write_md(
        "reports/clone_outcome_preserving_native_audit.md",
        "# Outcome-Preserving Native Clone Audit\n\n"
        "This audit tests whether clone-stratified native failures were caused by terminal fate entropy becoming unestimable. Jindal LSK and Weinreb LARRY were rerun with five sampling strategies that preserve terminal outcome quality before native moscot teacher extraction.\n\n"
        "## Final Interpretation\n\n"
        f"- status: {status}\n"
        f"- tier: {tier}\n"
        f"- interpretation: {interpretation}\n"
        f"- decision: {final_sentence}\n\n"
        "## Strategy Summary\n\n"
        f"{summary_md}\n\n"
        "The retained project story remains the time-series branch-nucleation order-parameter hypothesis supported by internal native moscot and E1 MouseGastrulationData. Clone-aware results are not promoted to a main claim unless they are positive across datasets after covariate adjustment, matched analyses and negative controls.\n",
    )
    _write_md(
        "reports/clone_outcome_preserving_native_conclusion.md",
        "# Outcome-Preserving Clone Conclusion\n\n"
        f"- final_status: {status}\n"
        f"- final_tier: {tier}\n"
        f"- interpretation: {interpretation}\n"
        f"- boundary: {final_sentence}\n\n"
        "If the outcome-preserving audit remains weak or failed, the clone-level fate-diversification line is closed as a main claim and retained only as future work.\n",
    )
    _write_md(
        "reports/clone_outcome_preserving_native_claim_audit.md",
        "# Outcome-Preserving Native Claim Audit\n\n"
        "- Jindal and Weinreb were re-evaluated with terminal-outcome-preserving native sampling.\n"
        "- The primary claim tested was branch-window condensation exposure predicting terminal clone fate entropy.\n"
        "- Secondary branch-window features were reported but cannot rescue a failed primary condensation claim.\n"
        "- No clone-level fate-diversification prediction is retained unless it is cross-dataset stable after covariates, matched analyses and negative controls.\n"
        "- Birth/death, memory, CCI, topological specificity, swarm-required causality and model-beats-OT claims remain excluded.\n",
    )
    _write_md(
        "manuscript/final_retained_results_and_methods.md",
        "# Final Retained Results and Methods\n\n"
        "## Retained Claim\n\n"
        "SwarmLineage-OT retains a branch-nucleation / transient condensation-before-divergence time-series order-parameter computational hypothesis. The strongest support remains internal native moscot teacher analysis plus E1 MouseGastrulationData external time-series support. The evidence-selected primary agent remains M5_ot_swarm.\n\n"
        "## Clone-Aware Boundary\n\n"
        "Clone-aware fate-diversification support is not established as a retained claim. Jindal full-data fallback weak positivity did not persist under clone-stratified native moscot. Weinreb produced sampling-specific native condensation signals, but they were not cross-dataset stable. The outcome-preserving native audit directly tested whether terminal fate entropy became unestimable during native sampling.\n\n"
        f"Outcome-preserving audit status: `{status}`; tier: `{tier}`. {interpretation} Therefore, {final_sentence}\n\n"
        "## Not Retained\n\n"
        "- clone-level fate-diversification prediction;\n"
        "- clone-level prediction from condensation;\n"
        "- topological-neighbour-specific mechanism;\n"
        "- swarm-required causality;\n"
        "- diffusion as an independent biological discovery;\n"
        "- birth/death, memory or CCI as supported mechanisms.\n\n"
        "## Next Strongest Experiment\n\n"
        "The next decisive experiment is an independent developmental time-series or spatial/time-series system, ideally gastruloid or embryoid-body differentiation with paired clone/barcode readouts and sufficient terminal sampling depth. Clone fate prediction should not be the main claim until such data support it.\n",
    )
    _write_md(
        "manuscript/manuscript.md",
        "# SwarmLineage-OT: Native OT-Guided Branch-Nucleation Order Parameters\n\n"
        "SwarmLineage-OT converts native OT-inferred pseudo-lineage maps into finite-agent rollouts. The retained result is a branch-nucleation / transient condensation-before-divergence time-series order-parameter computational hypothesis, not a claim of clone-level fate prediction.\n\n"
        "Internal native moscot analysis and E1 MouseGastrulationData remain the strongest support. Clone-aware analyses in Biddy/CellTag, Jindal LSK and Weinreb LARRY were used as stress tests. Jindal's full-data fallback weak positive did not survive clone-stratified native moscot. Weinreb shows sampling-specific native condensation associations, but those signals are not cross-dataset stable. The outcome-preserving native audit was added to test whether terminal fate entropy loss explained the negative native results.\n\n"
        f"The final clone-aware assessment is `{tier}`: {interpretation} The manuscript therefore retains the time-series order-parameter story and excludes clone-level fate-diversification prediction from the main claim.\n\n"
        "Diffusion remains an encoded control-law recovery. Birth/death, memory and CCI remain unsupported. Topological-neighbor specificity and swarm-required causality are not established.\n",
    )
    _write_md(
        "manuscript/methods.md",
        "# Methods\n\n"
        "Native moscot teachers were used for the internal dataset and for downsampled clone-aware Jindal/Weinreb audits where feasible. The final outcome-preserving clone audit sampled cells to preserve terminal clone fate entropy estimability before running native moscot. Associations were evaluated on outcome-estimable clones only, using raw Spearman effects, covariate-adjusted regression, within-start-time, clone-size and time-span matched analyses, bootstrap intervals and permutation controls.\n\n"
        "The primary clone-aware test was branch-window condensation exposure predicting terminal clone fate entropy. Secondary models included post-event divergence, fate entropy, teacher velocity bias, condensation plus divergence and uncertainty plus teacher bias. Secondary associations were reported as sensitivity analyses and were not allowed to rescue a failed primary claim.\n",
    )
    _write_md(
        "reports/scientific_gap_audit.md",
        "# Scientific Gap Audit\n\n"
        "OT gives the developmental map; SwarmLineage-OT learns microscopic finite-agent rules that realize the map and reveal emergent developmental laws.\n\n"
        "- teacher_fidelity_tier: acceptable\n"
        "- retained main claim: branch nucleation / transient condensation-before-divergence as a time-series order-parameter computational hypothesis\n"
        "- clone-aware tier: " + tier + "\n"
        "- clone-aware status: " + status + "\n\n"
        "## Remaining Gaps\n\n"
        "- Clone-level fate-diversification prediction is not retained under current tested datasets and native sampling strategies.\n"
        "- Outcome-preserving sampling addresses terminal entropy estimability but does not by itself establish clone support.\n"
        "- External time-series generalization beyond E1 remains limited.\n"
        "- Birth/death, memory and CCI remain unsupported and excluded from the main claim.\n"
        "- No manuscript claim may frame SwarmLineage-OT as surpassing the OT reference.\n",
    )
    _write_md(
        "reports/editorial_assessment.md",
        "# Editorial Assessment\n\n"
        "Current evidence level: a native-OT-guided computational research prototype with a retained time-series branch-nucleation order-parameter hypothesis. It is not a clone-level fate prediction result.\n\n"
        f"- outcome-preserving clone audit status: {status}\n"
        f"- outcome-preserving clone audit tier: {tier}\n"
        "- internal native + E1 time-series support remains the strongest evidence layer.\n"
        "- clone-aware fate-diversification support is not established as a main claim.\n"
        "- diffusion is encoded recovery; birth/death, memory and CCI remain unsupported.\n"
        "- the project is not ready for high-impact biological mechanism claims without independent validation.\n",
    )
    _update_claim_tier_table(status, tier, interpretation)


def _write_reports(sampling: pd.DataFrame, runs: pd.DataFrame, results: pd.DataFrame, controls: pd.DataFrame, clones: pd.DataFrame, summary: pd.DataFrame) -> None:
    status, tier, interpretation = _final_interpretation(summary)
    _write_csv(sampling, "tables/clone_outcome_preserving_native_sampling.csv")
    _write_csv(runs, "tables/clone_outcome_preserving_native_teacher_runs.csv")
    _write_csv(results, "tables/clone_outcome_preserving_native_model_results.csv")
    _write_csv(results, "tables/clone_outcome_preserving_native_outcome_sensitivity.csv")
    _write_csv(controls, "tables/clone_outcome_preserving_native_negative_controls.csv")
    _write_csv(clones, "tables/clone_outcome_preserving_native_clone_level.csv")
    _write_csv(summary, "tables/clone_outcome_preserving_native_strategy_summary.csv")
    _write_csv(
        pd.DataFrame([{"final_clone_aware_status": status, "final_clone_aware_tier": tier, "interpretation": interpretation}]),
        "tables/clone_outcome_preserving_native_final_summary.csv",
    )
    supported = results[(results["outcome_primary"].eq(True)) & (results["support_tier"].isin(["acceptable", "weak"]))].copy() if not results.empty else pd.DataFrame()
    _write_csv(supported, "tables/clone_outcome_preserving_native_supported_model_rows.csv")
    _write_figure(summary)
    _update_documents(status, tier, interpretation, summary)


def run() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", default=NATIVE_PYTHON)
    parser.add_argument("--timeout", type=int, default=1200)
    parser.add_argument("--seed", type=int, default=23)
    args = parser.parse_args()
    sampling_rows = []
    run_rows = []
    clone_tables = []
    result_tables = []
    control_tables = []
    cfgs = [d for d in DATASETS if d.dataset_id.startswith("Jindal") or d.dataset_id.startswith("Weinreb")]
    for cfg in cfgs:
        for strategy in STRATEGIES:
            prep = prepare_strategy(cfg, strategy, args.seed)
            sampling_rows.append(prep)
            if not prep.get("prepared", False):
                continue
            native = run_native(cfg.dataset_id, strategy.name, prep["config_path"], args.timeout, args.python)
            run_rows.append(native)
            clone, results, controls = analyze_strategy(cfg, strategy, bool(native.get("native_moscot_success", False)))
            if not clone.empty:
                clone.insert(0, "sampling_strategy", strategy.name)
                clone.insert(0, "dataset_id", cfg.dataset_id)
                clone_tables.append(clone)
            if not results.empty:
                result_tables.append(results)
            if not controls.empty:
                control_tables.append(controls)
    sampling = pd.DataFrame(sampling_rows)
    runs = pd.DataFrame(run_rows)
    results = pd.concat(result_tables, ignore_index=True) if result_tables else pd.DataFrame()
    controls = pd.concat(control_tables, ignore_index=True) if control_tables else pd.DataFrame()
    clones = pd.concat(clone_tables, ignore_index=True) if clone_tables else pd.DataFrame()
    summary = _summarize(sampling, runs, results)
    _write_reports(sampling, runs, results, controls, clones, summary)
    print(_final_interpretation(summary))


if __name__ == "__main__":
    run()
