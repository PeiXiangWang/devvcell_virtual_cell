from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances

from src.utils.config import ensure_dir, load_config, write_json


def _native_status(timeout: int) -> dict:
    code = "import moscot, json; print(json.dumps({'version': getattr(moscot, '__version__', 'unknown')}))"
    try:
        out = subprocess.check_output([sys.executable, "-c", code], text=True, stderr=subprocess.STDOUT, timeout=timeout)
        return {"available": True, "detail": out.strip(), "timeout_seconds": timeout}
    except subprocess.TimeoutExpired:
        return {"available": False, "detail": f"import timed out after {timeout}s", "timeout_seconds": timeout}
    except Exception as exc:
        return {"available": False, "detail": f"{type(exc).__name__}: {exc}", "timeout_seconds": timeout}


def _scaled_cost(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    cost = pairwise_distances(x, y, metric="sqeuclidean")
    positive = cost[cost > 0]
    scale = float(np.median(positive)) if positive.size else 1.0
    return cost / max(scale, 1e-8)


def _sinkhorn_numpy(a: np.ndarray, b: np.ndarray, cost: np.ndarray, reg: float, max_iter: int = 800, tol: float = 1e-7) -> np.ndarray:
    kernel = np.exp(-cost / max(reg, 1e-8))
    kernel = np.maximum(kernel, 1e-300)
    u = np.ones_like(a)
    v = np.ones_like(b)
    for _ in range(max_iter):
        u_prev = u.copy()
        u = a / np.maximum(kernel @ v, 1e-300)
        v = b / np.maximum(kernel.T @ u, 1e-300)
        if np.max(np.abs(u - u_prev)) < tol:
            break
    plan = (u[:, None] * kernel) * v[None, :]
    return plan / np.maximum(plan.sum(), 1e-300)


def _stratified_idx(obs: pd.DataFrame, max_n: int, seed: int, key: str) -> np.ndarray:
    if obs.shape[0] <= max_n:
        return np.arange(obs.shape[0])
    rng = np.random.default_rng(seed)
    labels = obs[key].astype(str) if key in obs else pd.Series("all", index=obs.index)
    selected: list[int] = []
    per = max(10, max_n // max(1, labels.nunique()))
    for _, idx in labels.groupby(labels).groups.items():
        pos = obs.index.get_indexer(idx)
        selected.extend(rng.choice(pos, size=min(per, len(pos)), replace=False).tolist())
    if len(selected) < max_n:
        rem = np.setdiff1d(np.arange(obs.shape[0]), np.asarray(selected), assume_unique=False)
        selected.extend(rng.choice(rem, size=min(max_n - len(selected), len(rem)), replace=False).tolist())
    if len(selected) > max_n:
        selected = rng.choice(np.asarray(selected), size=max_n, replace=False).tolist()
    return np.asarray(sorted(selected), dtype=int)


def compute_pair_coupling(
    z: np.ndarray,
    obs: pd.DataFrame,
    t0: float,
    t1: float,
    time_key: str,
    cell_type_key: str,
    max_cells: int,
    epsilon: float,
    seed: int,
) -> dict:
    idx0_all = np.where(obs[time_key].to_numpy(dtype=float) == float(t0))[0]
    idx1_all = np.where(obs[time_key].to_numpy(dtype=float) == float(t1))[0]
    sub0 = obs.iloc[idx0_all]
    sub1 = obs.iloc[idx1_all]
    idx0 = idx0_all[_stratified_idx(sub0, max_cells, seed, cell_type_key)]
    idx1 = idx1_all[_stratified_idx(sub1, max_cells, seed + 1, cell_type_key)]
    x = z[idx0]
    y = z[idx1]
    cost = _scaled_cost(x, y)
    a = np.full(x.shape[0], 1.0 / x.shape[0])
    b = np.full(y.shape[0], 1.0 / y.shape[0])
    method = "numpy_sinkhorn"
    try:
        plan = _sinkhorn_numpy(a, b, cost, reg=epsilon)
        if not np.all(np.isfinite(plan)) or plan.sum() <= 0:
            raise FloatingPointError("non-finite plan")
    except Exception:
        method = "softmax_fallback"
        logits = -cost / max(epsilon, 1e-6)
        logits -= logits.max(axis=1, keepdims=True)
        row = np.exp(logits)
        row /= np.maximum(row.sum(axis=1, keepdims=True), 1e-12)
        plan = row / row.sum()
    row_sum = np.maximum(plan.sum(axis=1, keepdims=True), 1e-12)
    transition = plan / row_sum
    bary = transition @ y
    entropy = -np.sum(np.clip(transition, 1e-12, 1.0) * np.log(np.clip(transition, 1e-12, 1.0)), axis=1)
    entropy = entropy / max(np.log(transition.shape[1]), 1.0)
    full_counts0 = obs.iloc[idx0_all][cell_type_key].astype(str).value_counts()
    full_counts1 = obs.iloc[idx1_all][cell_type_key].astype(str).value_counts()
    source_types = obs.iloc[idx0][cell_type_key].astype(str).to_numpy()
    growth = np.array([
        (full_counts1.get(ct, 0.5) + 0.5) / (full_counts0.get(ct, 0.5) + 0.5) for ct in source_types
    ])
    target_types = obs.iloc[idx1][cell_type_key].astype(str).to_numpy()
    return {
        "source_indices": idx0.astype(np.int64),
        "target_indices": idx1.astype(np.int64),
        "plan": plan.astype(np.float32),
        "transition": transition.astype(np.float32),
        "barycentric": bary.astype(np.float32),
        "entropy": entropy.astype(np.float32),
        "growth": growth.astype(np.float32),
        "source_time": float(t0),
        "target_time": float(t1),
        "source_types": source_types,
        "target_types": target_types,
        "method": method,
        "mean_cost": float(np.sum(plan * cost)),
    }


def run_ot(cfg: dict, label: str = "moscot") -> dict:
    ensure_dir(cfg.get("couplings_dir", "processed/ot_couplings"))
    adata = ad.read_h5ad(cfg["adata_path"])
    z = np.asarray(adata.obsm[cfg.get("latent_key", "X_pca")], dtype=float)
    obs = adata.obs.copy()
    time_key = cfg.get("time_key", "time_numeric")
    cell_type_key = cfg.get("cell_type_key", "lineage")
    times = sorted(pd.to_numeric(obs[time_key], errors="coerce").dropna().unique())
    rows = []
    for pair_i, (t0, t1) in enumerate(zip(times[:-1], times[1:])):
        result = compute_pair_coupling(
            z=z,
            obs=obs,
            t0=float(t0),
            t1=float(t1),
            time_key=time_key,
            cell_type_key=cell_type_key,
            max_cells=int(cfg.get("max_cells_per_time", 650)),
            epsilon=float(cfg.get("epsilon", 0.08)),
            seed=int(cfg.get("random_seed", 17)) + pair_i,
        )
        out = Path(cfg.get("couplings_dir", "processed/ot_couplings")) / f"{label}_t{t0:g}_to_t{t1:g}.npz"
        np.savez_compressed(
            out,
            source_indices=result["source_indices"],
            target_indices=result["target_indices"],
            plan=result["plan"],
            transition=result["transition"],
            barycentric=result["barycentric"],
            entropy=result["entropy"],
            growth=result["growth"],
            source_types=result["source_types"],
            target_types=result["target_types"],
            source_time=np.array(result["source_time"]),
            target_time=np.array(result["target_time"]),
            method=np.array(result["method"]),
        )
        rows.append(
            {
                "method_label": label,
                "source_time": t0,
                "target_time": t1,
                "file": str(out),
                "n_source": int(result["source_indices"].size),
                "n_target": int(result["target_indices"].size),
                "solver": result["method"],
                "transport_cost": result["mean_cost"],
                "mean_entropy": float(np.mean(result["entropy"])),
                "mean_growth": float(np.mean(result["growth"])),
            }
        )
    index = pd.DataFrame(rows)
    index_path = Path(cfg.get("couplings_dir", "processed/ot_couplings")) / f"{label}_coupling_index.csv"
    index.to_csv(index_path, index=False)
    return {"index_path": str(index_path), "pairs": rows}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/ot_teacher.yaml")
    args = parser.parse_args()
    cfg = load_config(args.config)
    ensure_dir("logs")
    native = _native_status(int(cfg.get("native_moscot_timeout_seconds", 45)))
    summary = run_ot(cfg, label="moscot")
    summary["native_moscot_status"] = native
    summary["native_moscot_used"] = False
    summary["note"] = "Native moscot is recorded but this prototype uses auditable POT/SciPy fallback couplings for reproducible quick execution."
    write_json(Path(cfg.get("couplings_dir", "processed/ot_couplings")) / "moscot_run_summary.json", summary)
    print(json.dumps({"moscot_fallback_pairs": len(summary["pairs"]), "native": native}, indent=2))


if __name__ == "__main__":
    main()
