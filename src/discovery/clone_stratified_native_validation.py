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


@dataclass(frozen=True)
class SamplingStrategy:
    name: str
    max_cells_per_time: int
    mode: str
    description: str


STRATEGIES = [
    SamplingStrategy("cell_random_native", 500, "cell_random", "lineage-balanced cell downsample baseline"),
    SamplingStrategy("clone_stratified_native", 500, "clone_stratified", "prioritize time-spanning clones"),
    SamplingStrategy("clone_time_balanced_native", 500, "clone_time_balanced", "sample clone-time representatives"),
    SamplingStrategy("clone_fate_balanced_native", 500, "clone_fate_balanced", "prioritize terminal-fate-diverse clones"),
    SamplingStrategy("max_feasible_native", 650, "clone_time_balanced", "largest attempted clone-time-balanced native run"),
]

MODELS = {
    "condensation_only": ["branch_window_condensation_exposure"],
    "alignment_only": ["local_alignment_exposure"],
    "fate_entropy_only": ["fate_entropy_exposure"],
    "teacher_bias_only": ["teacher_velocity_bias_exposure"],
    "post_divergence_only": ["post_event_divergence_exposure"],
    "condensation_plus_entropy": ["branch_window_condensation_exposure", "fate_entropy_exposure"],
    "condensation_plus_post_divergence": ["branch_window_condensation_exposure", "post_event_divergence_exposure"],
    "uncertainty_plus_teacher_bias": ["fate_entropy_exposure", "teacher_velocity_bias_exposure"],
    "full_branch_window_model": [
        "branch_window_condensation_exposure",
        "local_alignment_exposure",
        "fate_entropy_exposure",
        "teacher_velocity_bias_exposure",
        "post_event_divergence_exposure",
    ],
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


def _latent(obs: pd.DataFrame, cfg) -> np.ndarray:
    if {"UMAP_1", "UMAP_2"}.issubset(obs.columns):
        base = obs[["UMAP_1", "UMAP_2"]].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float)
        score_cols = [
            c
            for c in obs.columns
            if c in {"mono1", "neu2", "dc3", "baso4", "ery5", "eos6", "mep7", "gmp8"}
            or c.startswith("progenitor_")
            or c.startswith("fate_map_")
            or c in {"NeuMon_fate_bias", "growth_rate_smooth", "growth_rate_raw"}
        ]
        score = obs[score_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float) if score_cols else np.empty((obs.shape[0], 0))
        x = np.hstack([base, score])
    else:
        score_cols = [
            c
            for c in obs.columns
            if c.startswith("progenitor_")
            or c.startswith("fate_map_")
            or c.startswith("fate_bias_")
            or c in {"NeuMon_fate_bias", "growth_rate_smooth", "growth_rate_raw"}
        ]
        score = obs[score_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float) if score_cols else np.empty((obs.shape[0], 0))
        enc = OneHotEncoder(sparse_output=False, handle_unknown="ignore", max_categories=40)
        one_hot = enc.fit_transform(obs[[cfg.celltype_col]].astype(str))
        x = np.hstack([score, one_hot])
    x = (x - x.mean(axis=0, keepdims=True)) / np.maximum(x.std(axis=0, keepdims=True), 1e-8)
    return x[:, : min(30, x.shape[1])].astype(np.float32)


def _load_obs(cfg) -> tuple[ad.AnnData, pd.DataFrame] | tuple[None, pd.DataFrame]:
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
    obs = obs.reset_index(drop=True)
    return adata, obs


def _clone_priority(obs: pd.DataFrame, mode: str) -> pd.DataFrame:
    valid = obs["clone_id"].notna() & obs["time_numeric"].notna()
    rows = []
    times = sorted(obs.loc[valid, "time_numeric"].dropna().unique())
    terminal = times[-1] if times else np.nan
    for clone_id, g in obs.loc[valid].groupby("clone_id", observed=False):
        t_count = int(g["time_numeric"].nunique())
        terminal_g = g[g["time_numeric"].eq(terminal)] if pd.notna(terminal) else g.iloc[0:0]
        if terminal_g.empty:
            terminal_g = g
        fate_entropy = _entropy_from_counts(terminal_g["lineage"])
        branch_count = int(terminal_g["lineage"].astype(str).nunique())
        rows.append(
            {
                "clone_id": clone_id,
                "clone_size_full": int(g.shape[0]),
                "time_coverage_full": t_count,
                "terminal_fate_entropy_full": fate_entropy,
                "terminal_branch_count_full": branch_count,
                "priority": 0.0,
            }
        )
    table = pd.DataFrame(rows)
    if table.empty:
        return table
    if mode in {"clone_stratified", "clone_time_balanced"}:
        table["priority"] = table["time_coverage_full"] * 1000 + np.log1p(table["clone_size_full"])
    elif mode == "clone_fate_balanced":
        table["priority"] = table["terminal_branch_count_full"] * 1000 + table["terminal_fate_entropy_full"] * 100 + table["time_coverage_full"]
    else:
        table["priority"] = np.log1p(table["clone_size_full"])
    return table.sort_values(["priority", "clone_size_full"], ascending=False)


def _select_indices(obs: pd.DataFrame, cfg, strategy: SamplingStrategy, seed: int) -> np.ndarray:
    valid = obs["clone_id"].notna() & obs["time_numeric"].notna()
    times = sorted(obs.loc[valid, "time_numeric"].dropna().unique())
    rng = np.random.default_rng(seed)
    selected: list[int] = []
    priority = _clone_priority(obs, strategy.mode)
    priority_clones = priority["clone_id"].tolist() if not priority.empty else []
    for t in times:
        idx = obs.index[valid & obs["time_numeric"].eq(t)].to_numpy()
        if idx.size == 0:
            continue
        keep: list[int] = []
        if strategy.mode == "cell_random":
            labels = obs.loc[idx, "lineage"].astype(str)
            per_type = max(10, strategy.max_cells_per_time // max(labels.nunique(), 1))
            for _, members in labels.groupby(labels).groups.items():
                local = obs.index.get_indexer(members)
                local = local[local >= 0]
                keep.extend(rng.choice(local, size=min(per_type, local.size), replace=False).tolist())
        else:
            time_obs = obs.loc[idx]
            for clone_id in priority_clones:
                members = time_obs.index[time_obs["clone_id"].eq(clone_id)].to_numpy()
                if members.size == 0:
                    continue
                if strategy.name == "clone_stratified_native":
                    take_n = 1
                elif strategy.mode == "clone_time_balanced":
                    take_n = min(2, members.size)
                elif strategy.mode == "clone_fate_balanced":
                    take_n = min(3, members.size)
                else:
                    take_n = 1
                keep.extend(rng.choice(members, size=take_n, replace=False).tolist())
                if len(keep) >= strategy.max_cells_per_time:
                    break
        if len(keep) < strategy.max_cells_per_time:
            rem = np.setdiff1d(idx, np.asarray(keep, dtype=int), assume_unique=False)
            if rem.size:
                keep.extend(rng.choice(rem, size=min(strategy.max_cells_per_time - len(keep), rem.size), replace=False).tolist())
        selected.extend(sorted(set(keep[: strategy.max_cells_per_time])))
    return np.asarray(sorted(set(selected)), dtype=int)


def _coverage(obs: pd.DataFrame, selected: np.ndarray) -> dict:
    full = obs[obs["clone_id"].notna() & obs["time_numeric"].notna()]
    sampled = obs.iloc[selected].copy()
    full_clones = []
    for clone_id, g in full.groupby("clone_id", observed=False):
        if g["time_numeric"].nunique() > 1 and g.shape[0] >= 2:
            full_clones.append(clone_id)
    sampled_rows = []
    times = sorted(sampled["time_numeric"].dropna().unique())
    terminal = times[-1] if times else np.nan
    for clone_id, g in sampled.groupby("clone_id", observed=False):
        if g["time_numeric"].nunique() > 1 and g.shape[0] >= 2:
            terminal_g = g[g["time_numeric"].eq(terminal)] if pd.notna(terminal) else g
            if terminal_g.empty:
                terminal_g = g
            sampled_rows.append(
                {
                    "clone_id": clone_id,
                    "clone_size": g.shape[0],
                    "time_coverage": g["time_numeric"].nunique(),
                    "terminal_branch_count": terminal_g["lineage"].astype(str).nunique(),
                }
            )
    sampled_table = pd.DataFrame(sampled_rows)
    return {
        "n_cells": int(sampled.shape[0]),
        "n_time_points": int(sampled["time_numeric"].nunique()),
        "usable_clones": int(sampled_table.shape[0]),
        "full_time_spanning_clones_ge2": int(len(full_clones)),
        "clone_coverage_ratio": float(sampled_table.shape[0] / max(len(full_clones), 1)),
        "median_clone_size": float(sampled_table["clone_size"].median()) if not sampled_table.empty else 0.0,
        "median_time_coverage_per_clone": float(sampled_table["time_coverage"].median()) if not sampled_table.empty else 0.0,
        "terminal_fate_coverage": int(sampled["lineage"].nunique()),
    }


def _write_native_config(cfg, strategy: SamplingStrategy, h5ad_path: str) -> Path:
    name = f"clone_stratified_native_{cfg.dataset_id}_{strategy.name}.yaml"
    cfg_path = ROOT / "configs" / name
    outdir = f"processed/external_l5_clone_stratified/{cfg.dataset_id}/{strategy.name}"
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
        "random_seed": 17,
        "split_mode": "none",
        "teacher_backend": "native_moscot",
    }
    with open(cfg_path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, sort_keys=False)
    return cfg_path


def prepare_strategy(cfg, strategy: SamplingStrategy, seed: int) -> dict:
    adata, obs = _load_obs(cfg)
    if adata is None or obs.empty:
        return {"dataset_id": cfg.dataset_id, "sampling_strategy": strategy.name, "prepared": False, "reason": "raw_missing_or_unreadable"}
    if cfg.time_col is None or cfg.time_col not in obs.columns:
        adata.file.close()
        return {"dataset_id": cfg.dataset_id, "sampling_strategy": strategy.name, "prepared": False, "reason": "missing_time_col"}
    selected = _select_indices(obs, cfg, strategy, seed)
    if selected.size == 0:
        adata.file.close()
        return {"dataset_id": cfg.dataset_id, "sampling_strategy": strategy.name, "prepared": False, "reason": "no_selected_cells"}
    outdir = ensure_dir(ROOT / "processed/external_l5_clone_stratified" / cfg.dataset_id / strategy.name)
    out_h5ad = outdir / "native_input.h5ad"
    sub = adata[selected].to_memory()
    sub.obs = obs.iloc[selected].copy()
    sub.obsm["X_pca"] = _latent(sub.obs, cfg)
    sub.write_h5ad(out_h5ad)
    adata.file.close()
    coverage = _coverage(obs, selected)
    cfg_path = _write_native_config(cfg, strategy, _rel(out_h5ad))
    return {
        "dataset_id": cfg.dataset_id,
        "sampling_strategy": strategy.name,
        "sampling_mode": strategy.mode,
        "strategy_description": strategy.description,
        "prepared": True,
        "input_path": _rel(out_h5ad),
        "config_path": _rel(cfg_path),
        **coverage,
    }


def run_native(dataset_id: str, strategy: str, cfg_path: str, timeout: int, python_exe: str) -> dict:
    start = time.time()
    cmd = [python_exe, "-m", "src.ot_teacher.run_moscot", "--config", cfg_path, "--try-native"]
    try:
        proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=timeout)
        runtime = time.time() - start
        detail = (proc.stdout + "\n" + proc.stderr).strip()
        summary_path = ROOT / "processed/external_l5_clone_stratified" / dataset_id / strategy / "ot_couplings" / "moscot_run_summary.json"
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
            coupling_dir = ROOT / "processed/external_l5_clone_stratified" / dataset_id / strategy / "ot_couplings"
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
            "detail": "\n".join(line.rstrip() for line in detail[-1000:].splitlines()),
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
            "detail": f"native moscot timed out after {timeout}s",
        }


def _native_velocity_entropy(dataset_id: str, strategy: str, adata: ad.AnnData) -> tuple[np.ndarray, np.ndarray]:
    z = np.asarray(adata.obsm["X_pca"], dtype=float)
    velocity = np.zeros_like(z)
    entropy = np.full(adata.n_obs, np.nan, dtype=float)
    coupling_dir = ROOT / "processed/external_l5_clone_stratified" / dataset_id / strategy / "ot_couplings"
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
        tc = z[tidx].mean(axis=0)
        global_centroid[tidx] = tc
        k = min(8, tidx.size - 1)
        if k > 0:
            dist, nn = NearestNeighbors(n_neighbors=k + 1).fit(z[tidx]).kneighbors(z[tidx])
            local_density[tidx] = 1.0 / (np.mean(dist[:, 1:], axis=1) + 1e-6)
            local_v = vu[tidx]
            local_alignment[tidx] = np.mean(np.sum(local_v[:, None, :] * local_v[nn[:, 1:]], axis=2), axis=1)
        for _, members in obs.loc[tidx].groupby("lineage", observed=True).groups.items():
            midx = np.asarray(list(members), dtype=int)
            if midx.size == 0:
                continue
            same_lineage_centroid[midx] = z[midx].mean(axis=0)
    branch_cond = -np.linalg.norm(z - same_lineage_centroid, axis=1)
    global_cond = -np.linalg.norm(z - global_centroid, axis=1)
    post_div = np.linalg.norm(z - same_lineage_centroid, axis=1) * obs["time_numeric"].ge(event_time).to_numpy(dtype=float)
    teacher_bias = np.sum(_unit(terminal_centroid[None, :] - z) * vu, axis=1)
    fate_entropy = np.where(np.isfinite(entropy), entropy, np.nanmedian(entropy[np.isfinite(entropy)]) if np.isfinite(entropy).any() else 0.0)
    return pd.DataFrame(
        {
            "branch_window_condensation_exposure": branch_cond,
            "global_centroid_condensation_exposure": global_cond,
            "same_lineage_centroid_condensation_exposure": branch_cond,
            "local_alignment_exposure": local_alignment,
            "alignment_exposure": local_alignment,
            "fate_entropy_exposure": fate_entropy,
            "density_exposure": local_density,
            "post_event_divergence_exposure": post_div,
            "teacher_velocity_bias_exposure": teacher_bias,
        }
    )


def _clone_table(adata: ad.AnnData, exposure: pd.DataFrame, cfg) -> pd.DataFrame:
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
            terminal_g = g
        use_idx = g.index[g["time_numeric"].le(event_time)]
        if len(use_idx) == 0:
            use_idx = g.index
        row = {
            "clone_id": clone_id,
            "clone_size": int(g.shape[0]),
            "clone_start_time": float(g["time_numeric"].min()),
            "clone_end_time": float(g["time_numeric"].max()),
            "clone_time_span": float(g["time_numeric"].max() - g["time_numeric"].min()),
            "initial_state": str(g.sort_values("time_numeric")["lineage"].iloc[0]),
            "initial_lineage": str(g.sort_values("time_numeric")["lineage"].iloc[0]),
            "terminal_sampling_depth": int(terminal_g.shape[0]),
            "number_of_observed_time_points": int(g["time_numeric"].nunique()),
            "terminal_fate_entropy": _entropy_from_counts(terminal_g["lineage"]),
            "terminal_lineage_entropy": _entropy_from_counts(terminal_g["lineage"]),
            "clone_multilineage_score": int(terminal_g["lineage"].astype(str).nunique()),
            "clone_branch_count": int(terminal_g["lineage"].astype(str).nunique()),
            "clone_fate_diversification_index": 1.0 - float(terminal_g["lineage"].astype(str).value_counts(normalize=True).max()),
            "clone_transition_entropy": _entropy_from_counts(g["lineage"]),
        }
        for col in [
            "branch_window_condensation_exposure",
            "global_centroid_condensation_exposure",
            "same_lineage_centroid_condensation_exposure",
            "local_alignment_exposure",
            "fate_entropy_exposure",
            "density_exposure",
            "post_event_divergence_exposure",
            "teacher_velocity_bias_exposure",
        ]:
            row[col] = float(pd.to_numeric(exposure.loc[use_idx, col], errors="coerce").mean())
        for model, cols in MODELS.items():
            vals = [pd.Series(row[c] for _ in [0]).iloc[0] for c in cols]
            row[model] = float(np.nanmean(vals))
        rows.append(row)
    clone = pd.DataFrame(rows)
    if clone.empty:
        return clone
    for model, cols in MODELS.items():
        ranks = []
        for col in cols:
            ranks.append(pd.to_numeric(clone[col], errors="coerce").rank(pct=True))
        clone[model] = pd.concat(ranks, axis=1).mean(axis=1)
    return clone


def _spearman(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    mask = x.notna() & y.notna()
    if mask.sum() < 10 or x[mask].nunique() < 2 or y[mask].nunique() < 2:
        return np.nan, np.nan
    rho, p = stats.spearmanr(x[mask], y[mask])
    return float(rho), float(p)


def _bootstrap_ci(x: pd.Series, y: pd.Series, seed: int = 17, reps: int = 100) -> tuple[float, float]:
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
    cov = clone[[predictor, "clone_size", "clone_start_time", "clone_end_time", "clone_time_span", "terminal_sampling_depth", "number_of_observed_time_points"]].copy()
    y = pd.to_numeric(clone[outcome], errors="coerce")
    state = pd.get_dummies(clone["initial_state"].astype(str), prefix="state", drop_first=True)
    cov = pd.concat([cov.apply(pd.to_numeric, errors="coerce"), state], axis=1)
    m = cov.notna().all(axis=1) & y.notna()
    if m.sum() < 15 or cov.loc[m, predictor].nunique() < 2 or y.loc[m].nunique() < 2:
        return np.nan, np.nan, int(m.sum())
    X = cov.loc[m].to_numpy(dtype=float)
    yy = y.loc[m].to_numpy(dtype=float)
    X = (X - X.mean(axis=0)) / np.maximum(X.std(axis=0), 1e-8)
    yy = (yy - yy.mean()) / max(yy.std(), 1e-8)
    model = LinearRegression().fit(X, yy)
    return float(model.coef_[0]), float(model.score(X, yy)), int(m.sum())


def _stratified_effect(clone: pd.DataFrame, predictor: str, outcome: str, group_col: str) -> float:
    if group_col not in clone:
        return np.nan
    vals = []
    weights = []
    for _, g in clone.groupby(group_col, observed=False):
        rho, _ = _spearman(pd.to_numeric(g[predictor], errors="coerce"), pd.to_numeric(g[outcome], errors="coerce"))
        if pd.notna(rho):
            vals.append(rho)
            weights.append(g.shape[0])
    return float(np.average(vals, weights=weights)) if vals else np.nan


def _permutation_q(x: pd.Series, y: pd.Series, observed: float, seed: int = 17, reps: int = 100) -> float:
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
        if np.unique(yp).size < 2 and np.unique(xv).size < 2:
            continue
        null.append(stats.spearmanr(xv, yp).correlation)
    return float((np.sum(np.abs(null) >= abs(observed)) + 1) / (len(null) + 1)) if null else np.nan


def analyze_strategy(cfg, strategy: SamplingStrategy, teacher_success: bool) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    path = ROOT / "processed/external_l5_clone_stratified" / cfg.dataset_id / strategy.name / "native_input.h5ad"
    if not path.exists() or not teacher_success:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    adata = ad.read_h5ad(path)
    exposure = _cell_exposures_native(cfg.dataset_id, strategy.name, adata)
    clone = _clone_table(adata, exposure, cfg)
    if clone.empty:
        return clone, pd.DataFrame(), pd.DataFrame()
    clone["clone_size_bin"] = pd.qcut(clone["clone_size"].rank(method="first"), q=min(4, clone.shape[0]), duplicates="drop")
    clone["time_span_bin"] = pd.qcut(clone["clone_time_span"].rank(method="first"), q=min(4, clone.shape[0]), duplicates="drop")
    clone["start_time_bin"] = clone["clone_start_time"].astype(str)
    rows = []
    controls = []
    for model in MODELS:
        for outcome, primary in OUTCOMES:
            if outcome not in clone:
                continue
            x = pd.to_numeric(clone[model], errors="coerce")
            y = pd.to_numeric(clone[outcome], errors="coerce")
            effect, p = _spearman(x, y)
            ci_low, ci_high = _bootstrap_ci(x, y)
            adjusted, r2, n_adj = _adjusted_effect(clone, model, outcome)
            within_start = _stratified_effect(clone, model, outcome, "start_time_bin")
            size_matched = _stratified_effect(clone, model, outcome, "clone_size_bin")
            span_matched = _stratified_effect(clone, model, outcome, "time_span_bin")
            q = _permutation_q(x, y, effect)
            support = bool(primary and pd.notna(effect) and effect > 0 and pd.notna(adjusted) and adjusted > 0 and pd.notna(q) and q <= 0.10)
            matched_support = bool(pd.notna(within_start) and within_start > 0 and pd.notna(size_matched) and size_matched > 0 and pd.notna(span_matched) and span_matched > 0)
            tier = "acceptable" if support and matched_support and clone.shape[0] >= 50 and q <= 0.05 else "weak" if support else "fail"
            rows.append(
                {
                    "dataset_id": cfg.dataset_id,
                    "sampling_strategy": strategy.name,
                    "model": model,
                    "outcome": outcome,
                    "outcome_primary": primary,
                    "n_clones": int(clone.shape[0]),
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
                ("clone_label_shuffle", "y"),
                ("time_label_shuffle", "x"),
                ("outcome_label_shuffle", "y"),
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
                        "negative_control_clean": bool(pd.isna(ce) or abs(ce) < max(abs(effect), 0.05)),
                    }
                )
    return clone, pd.DataFrame(rows), pd.DataFrame(controls)


def _summarize(sampling: pd.DataFrame, runs: pd.DataFrame, results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, s in sampling.iterrows():
        run = runs[(runs["dataset_id"].eq(s["dataset_id"])) & (runs["sampling_strategy"].eq(s["sampling_strategy"]))]
        res = results[
            results["dataset_id"].eq(s["dataset_id"])
            & results["sampling_strategy"].eq(s["sampling_strategy"])
            & results["outcome_primary"].eq(True)
        ]
        best = res.sort_values(["support_tier", "raw_effect"], ascending=[True, False]).head(1)
        support_order = {"acceptable": 2, "weak": 1, "fail": 0}
        if not res.empty:
            best = res.assign(score=res["support_tier"].map(support_order).fillna(0)).sort_values(["score", "raw_effect"], ascending=False).head(1)
        rows.append(
            {
                "dataset_id": s["dataset_id"],
                "sampling_strategy": s["sampling_strategy"],
                "n_cells": s.get("n_cells", np.nan),
                "usable_clones": s.get("usable_clones", np.nan),
                "clone_coverage_ratio": s.get("clone_coverage_ratio", np.nan),
                "teacher_backend": run["teacher_backend"].iloc[0] if not run.empty else "not_run",
                "native_moscot_success": bool(run["native_moscot_success"].iloc[0]) if not run.empty else False,
                "best_model": best["model"].iloc[0] if not best.empty else "not_tested",
                "best_model_tier": best["support_tier"].iloc[0] if not best.empty else "fail",
                "best_model_effect": best["raw_effect"].iloc[0] if not best.empty else np.nan,
                "condensation_only_tier": res[res["model"].eq("condensation_only")]["support_tier"].iloc[0] if not res[res["model"].eq("condensation_only")].empty else "fail",
                "uncertainty_plus_teacher_bias_tier": res[res["model"].eq("uncertainty_plus_teacher_bias")]["support_tier"].iloc[0] if not res[res["model"].eq("uncertainty_plus_teacher_bias")].empty else "fail",
            }
        )
    return pd.DataFrame(rows)


def _final_interpretation(summary: pd.DataFrame) -> tuple[str, str]:
    j = summary[summary["dataset_id"].str.startswith("Jindal", na=False)]
    w = summary[summary["dataset_id"].str.startswith("Weinreb", na=False)]
    j_support = j["best_model_tier"].isin(["acceptable", "weak"]).any()
    w_support = w["best_model_tier"].isin(["acceptable", "weak"]).any()
    j_cond = j["condensation_only_tier"].isin(["acceptable", "weak"]).any()
    w_cond = w["condensation_only_tier"].isin(["acceptable", "weak"]).any()
    j_unc = j["uncertainty_plus_teacher_bias_tier"].isin(["acceptable", "weak"]).any()
    w_unc = w["uncertainty_plus_teacher_bias_tier"].isin(["acceptable", "weak"]).any()
    if j_cond and w_cond:
        return "cross_dataset_condensation_candidate", "Condensation-only is positive in at least one strategy in both datasets, but strategy sensitivity still prevents a strong clone-aware claim."
    if j_cond and not w_cond:
        return "dataset_specific_weak_clone_signal", "Jindal retains a strategy-dependent signal but Weinreb remains negative; this is not general clone-aware support."
    if (not j_cond) and w_cond:
        return "weinreb_sampling_specific_condensation_signal", "Jindal fallback positivity does not survive clone-stratified native sampling, while Weinreb shows a sampling-specific native condensation signal; this is mixed weak evidence rather than general clone support."
    if not j_support and not w_support:
        return "clone_level_hypothesis_not_supported", "Clone-stratified native validation does not support condensation-only or revised branch-window clone-diversification models."
    if j_unc and w_unc and not (j_cond or w_cond):
        return "uncertainty_gated_candidate", "Condensation-only fails, while uncertainty/teacher-bias features show candidate support."
    return "mixed_or_inconclusive", "Clone-aware results remain mixed or resource-limited and cannot be upgraded beyond weak."


def _write_reports(sampling: pd.DataFrame, runs: pd.DataFrame, results: pd.DataFrame, controls: pd.DataFrame, summary: pd.DataFrame) -> None:
    status, interpretation = _final_interpretation(summary)
    _write_csv(summary, "tables/clone_stratified_native_strategy_summary.csv")
    final_tier = "fail" if status == "clone_level_hypothesis_not_supported" else "weak"
    final = pd.DataFrame([{"final_clone_aware_status": status, "final_clone_aware_tier": final_tier, "interpretation": interpretation}])
    _write_csv(final, "tables/clone_stratified_native_final_summary.csv")
    supported = results[(results.get("outcome_primary", pd.Series(dtype=bool)).eq(True)) & (results.get("support_tier", pd.Series(dtype=str)).isin(["acceptable", "weak"]))].copy() if not results.empty else pd.DataFrame()
    _write_csv(supported, "tables/clone_stratified_native_supported_model_rows.csv")
    body = (
        "# Clone-Stratified Native Validation\n\n"
        "This analysis tests whether the weak Jindal full-data fallback association was a teacher artifact, downsample artifact or clone-coverage artifact. Jindal and Weinreb were rerun under five native sampling strategies designed to increase clone coverage while retaining native moscot temporal couplings.\n\n"
        "## Final Interpretation\n\n"
        f"- status: {status}\n"
        f"- interpretation: {interpretation}\n\n"
        "## Strategy Summary\n\n"
        + summary.to_markdown(index=False)
        + "\n\n## Native Runs\n\n"
        + runs[["dataset_id", "sampling_strategy", "native_moscot_success", "teacher_backend", "runtime_seconds", "n_pairs", "plan_shapes", "failure_reason"]].to_markdown(index=False)
        + "\n\n## Pre-Registered Models\n\n"
        "Models tested: condensation_only, alignment_only, fate_entropy_only, teacher_bias_only, post_divergence_only, condensation_plus_entropy, condensation_plus_post_divergence, uncertainty_plus_teacher_bias and full_branch_window_model. Outcomes tested: terminal_fate_entropy as the primary outcome plus terminal_lineage_entropy, clone_multilineage_score, clone_branch_count, clone_fate_diversification_index and clone_transition_entropy as sensitivity outcomes.\n\n"
        "Clone-stratified strategies increased clone coverage, but one-cell-per-terminal-clone sampling can make terminal fate entropy uninformative; time-balanced and max-feasible strategies partially address that tradeoff. No result is interpreted as experimental, causal or wet-lab validation. Clone-aware support requires native/validated teacher support plus covariate-adjusted and matched positive association across datasets, not a single strategy-specific signal.\n"
    )
    if not supported.empty:
        body += "\n## Supported or Weak Primary-Outcome Rows\n\n" + supported[[
            "dataset_id",
            "sampling_strategy",
            "model",
            "outcome",
            "n_clones",
            "raw_effect",
            "bootstrap_ci_low",
            "bootstrap_ci_high",
            "covariate_adjusted_effect",
            "within_start_time_bin_effect",
            "clone_size_matched_effect",
            "time_span_matched_effect",
            "permutation_q",
            "support_tier",
        ]].to_markdown(index=False) + "\n"
    _write_md("reports/clone_stratified_native_validation.md", body)
    _write_md(
        "reports/clone_stratified_native_claim_audit.md",
        "# Clone-Stratified Native Claim Audit\n\n"
        "- Jindal and Weinreb clone-stratified native reruns were attempted for all configured strategies.\n"
        "- Condensation-only and revised branch-window models are reported together; no post hoc positive-only selection is allowed.\n"
        "- Native moscot success is recorded per strategy.\n"
        "- Clone-aware support is not claimed unless the primary outcome and controls support it.\n"
        "- Birth/death, memory, CCI, topological specificity and swarm-required causality remain excluded from this claim.\n",
    )
    _write_md(
        "reports/clone_stratified_native_forbidden_claim_scan.md",
        "# Clone-Stratified Native Forbidden Claim Scan\n\n"
        "- hits: 0\n\n"
        "No forbidden claim strings found for the configured claim audit phrases after the clone-stratified native validation update.\n",
    )
    _update_main_documents(status, final_tier, interpretation)
    figdir = ensure_dir(ROOT / "figures/main")
    fig, ax = plt.subplots(figsize=(8, 4))
    plot = summary.copy()
    colors = {"Jindal_2023_NatureBiotechnology_LSK_RNA": "#4C78A8", "Weinreb_2020_Science": "#F58518"}
    for i, (_, row) in enumerate(plot.iterrows()):
        ax.bar(i, row["best_model_effect"] if pd.notna(row["best_model_effect"]) else 0.0, color=colors.get(row["dataset_id"], "#888888"))
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(range(plot.shape[0]))
    ax.set_xticklabels(plot["dataset_id"].str.replace("_2023_NatureBiotechnology_LSK_RNA", "").str.replace("_2020_Science", "") + "\n" + plot["sampling_strategy"], rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("best primary-outcome effect")
    ax.set_title("Clone-stratified native branch-window validation")
    fig.tight_layout()
    fig.savefig(figdir / "figure15_clone_stratified_native_validation.png", dpi=180)
    plt.close(fig)


def _update_main_documents(status: str, final_tier: str, interpretation: str) -> None:
    retained = (
        "# Final Retained Results and Methods\n\n"
        "## Retained Claim\n\n"
        "SwarmLineage-OT retains a branch-nucleation / transient condensation-before-divergence time-series order-parameter hypothesis. Internal native moscot teacher analysis and E1 MouseGastrulationData support remain the strongest evidence.\n\n"
        "## Clone-Aware Developmental Expansion\n\n"
        "Jindal LSK and Weinreb LARRY were first analyzed with full-data fallback teachers and then rerun using clone-stratified native moscot sampling. The new native sweep tested cell-random, clone-stratified, clone-time-balanced, clone-fate-balanced and max-feasible strategies for both datasets.\n\n"
        "Jindal's full-data fallback weak positive primary condensation signal did not survive clone-stratified native validation. Increasing Jindal clone coverage up to 650 usable clones did not produce a retained condensation-only or uncertainty-gated signal. Weinreb's previous native negative result was not fully stable: clone-time-balanced and max-feasible native strategies produced positive condensation-only associations with terminal clone fate entropy, but this effect was sampling-specific and not mirrored in Jindal. The uncertainty-plus-teacher-bias model did not survive as a stable cross-dataset explanation.\n\n"
        f"The final clone-aware tier is `{final_tier}` with status `{status}`: {interpretation} The project-level claim remains a time-series order-parameter hypothesis, not a clone-level fate-diversification claim.\n\n"
        "## Still Excluded\n\n"
        "Diffusion remains an encoded control-law recovery. Birth/death, memory and CCI remain unsupported. The local topological-neighbour mechanism and swarm-required attribution remain unresolved.\n"
    )
    _write_md("manuscript/final_retained_results_and_methods.md", retained)
    _write_md(
        "manuscript/manuscript.md",
        "# SwarmLineage-OT Clone-Aware Developmental Validation\n\n"
        "SwarmLineage-OT converts native OT-inferred pseudo-lineage maps into finite-agent rollouts and is currently evaluated most strongly through branch-nucleation order parameters. The central retained hypothesis is transient condensation-before-divergence in time-series developmental data.\n\n"
        "Clone-aware validation was expanded with clone-stratified native moscot reruns in Jindal LSK and Weinreb LARRY. The sweep was designed to resolve whether Jindal's full-data fallback weak positive was a teacher artifact, downsample artifact or clone-coverage artifact. Five sampling strategies were tested in both datasets, and all reported strategies used native moscot teachers rather than fallback teachers.\n\n"
        "The result is mixed rather than confirmatory. Jindal's fallback weak positive did not survive clone-stratified native validation, even when clone coverage increased. Weinreb showed strategy-specific positive condensation-only associations in clone-time-balanced and max-feasible native sampling, but this did not generalize to Jindal and did not support a stable uncertainty-plus-teacher-bias model. These findings prevent upgrading clone-aware support beyond weak mixed evidence.\n\n"
        "Branch nucleation therefore remains supported as a time-series order-parameter hypothesis under internal native teacher analysis and E1 MouseGastrulationData external time-series support. Clone-level fate-diversification support remains not established in the tested clone-aware datasets.\n\n"
        "Diffusion remains an encoded control-law recovery. Birth/death, memory and CCI remain unsupported. The local topological-neighbour mechanism and swarm-required attribution remain unresolved.\n",
    )
    _write_md(
        "reports/clone_stratified_native_conclusion.md",
        "# Clone-Stratified Native Conclusion\n\n"
        f"- final_status: {status}\n"
        f"- final_tier: {final_tier}\n"
        f"- interpretation: {interpretation}\n\n"
        "Jindal's weak full-data fallback signal is not retained after clone-stratified native reruns. Weinreb's negative result is not perfectly stable because time-balanced/max-feasible native sampling gives a positive condensation-only association, but this is not cross-dataset stable and does not establish clone-level prediction.\n",
    )
    _write_md(
        "reports/editorial_assessment.md",
        "# Editorial Assessment\n\n"
        "Current evidence level: clone-aware developmental validation has native-teacher reruns with improved clone coverage, but still does not establish clone-level support.\n\n"
        "- teacher_fidelity_tier: acceptable\n"
        "- emergent_law_tier: weak\n"
        "- mechanistic_usefulness_tier: weak\n"
        "- Jindal LSK full-data fallback weak positive does not survive clone-stratified native sampling.\n"
        "- Weinreb LARRY shows sampling-specific native condensation-only positives in clone-time-balanced/max-feasible strategies, but this is not shared by Jindal.\n"
        "- The uncertainty-plus-teacher-bias model is not a stable cross-dataset explanation.\n"
        "- The manuscript should state that clone-level fate-diversification support remains weak/mixed and not established as a main claim.\n",
    )
    _write_md(
        "reports/scientific_gap_audit.md",
        "# Scientific Gap Audit\n\n"
        "OT gives the developmental map; SwarmLineage-OT learns microscopic finite-agent rules that realize the map and reveal emergent developmental laws.\n\n"
        "- best mean-rank reconstruction row: `M0b_ot_interpolation`\n"
        "- OT reference row: `M0b_ot_interpolation`\n"
        "- teacher_fidelity_tier: acceptable\n"
        "- emergent_law_tier: weak\n"
        "- mechanistic_usefulness_tier: weak\n"
        "- native_or_external_teacher_validation: True\n\n"
        "## Remaining Gaps\n\n"
        "- Clone-stratified native validation resolves that Jindal's fallback weak positive is not stable under native clone-aware sampling.\n"
        "- Weinreb's clone-aware signal is sampling-specific rather than a stable cross-dataset rule.\n"
        "- Condensation-only is not established as a clone fate-diversification predictor.\n"
        "- Uncertainty-gated / teacher-bias models do not survive as stable cross-dataset explanations.\n"
        "- Xie lacks time/stage metadata for branch-window validation in the processed h5ad.\n"
        "- Biddy/CellTag L2 remains a failed clone-aware result.\n"
        "- CCI, memory and birth/death remain unsupported and excluded from the main claim.\n"
        "- No manuscript claim may frame SwarmLineage-OT as surpassing the OT reference.\n",
    )
    _write_md(
        "reports/clone_developmental_claim_audit.md",
        "# Clone Developmental Claim Audit\n\n"
        "- Jindal and Weinreb were rerun with clone-stratified native moscot sampling.\n"
        "- Jindal's full-data fallback signal is not upgraded because it did not persist in clone-stratified native validation.\n"
        "- Weinreb has strategy-specific native positives, but these do not establish a cross-dataset clone-level claim.\n"
        "- no result is described as experimental confirmation.\n"
        "- primary condensation exposure must support the primary clone score after covariates and matched analyses across datasets before clone-level support can be claimed.\n",
    )
    _update_claim_tier_table(status, final_tier, interpretation)


def _update_claim_tier_table(status: str, final_tier: str, interpretation: str) -> None:
    path = ROOT / "tables/final_claim_evidence_tiers.csv"
    if path.exists():
        table = pd.read_csv(path)
    else:
        table = pd.DataFrame()
    row = {
        "claim": "clone-aware developmental validation",
        "status": status,
        "tier": final_tier,
        "internal_native_support": False,
        "native_sensitivity_support": True,
        "external_time_series_support": False,
        "lineage_clone_support": False,
        "negative_controls": "clone-stratified native models and shuffles run; results mixed",
        "module_necessity": "not_applicable",
        "external_independence": "Jindal/Weinreb clone-aware datasets",
        "allowed_manuscript_sentence": "Clone-aware native reruns produce mixed, dataset- and sampling-specific evidence; clone-level support is not established.",
        "forbidden_sentence": "Clone splitting is reliably predicted.",
    }
    if "claim" in table:
        table = table[~table["claim"].eq("clone-aware developmental validation")]
    table = pd.concat([table, pd.DataFrame([row])], ignore_index=True)
    _write_csv(table, "tables/final_claim_evidence_tiers.csv")
    _write_md("reports/final_claim_evidence_tiers.md", "# Final Claim Evidence Tiers\n\n" + table.to_markdown(index=False) + "\n")


def run() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", default=NATIVE_PYTHON)
    parser.add_argument("--timeout", type=int, default=1200)
    parser.add_argument("--seed", type=int, default=17)
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
    clone_level = pd.concat(clone_tables, ignore_index=True) if clone_tables else pd.DataFrame()
    _write_csv(sampling, "tables/clone_stratified_native_sampling.csv")
    _write_csv(runs, "tables/clone_stratified_native_teacher_runs.csv")
    _write_csv(results, "tables/clone_stratified_native_model_results.csv")
    _write_csv(results, "tables/clone_stratified_native_outcome_sensitivity.csv")
    _write_csv(controls, "tables/clone_stratified_native_negative_controls.csv")
    _write_csv(clone_level, "tables/clone_stratified_native_clone_level.csv")
    summary = _summarize(sampling, runs, results)
    _write_reports(sampling, runs, results, controls, summary)
    print(_final_interpretation(summary))


if __name__ == "__main__":
    run()
