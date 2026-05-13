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
    allowed_mask: np.ndarray | None = None,
) -> dict:
    idx0_all = np.where(obs[time_key].to_numpy(dtype=float) == float(t0))[0]
    idx1_all = np.where(obs[time_key].to_numpy(dtype=float) == float(t1))[0]
    if allowed_mask is not None:
        idx0_all = idx0_all[allowed_mask[idx0_all]]
        idx1_all = idx1_all[allowed_mask[idx1_all]]
    if idx0_all.size == 0 or idx1_all.size == 0:
        raise ValueError(f"Cannot build teacher edge {t0}->{t1}: one side has no allowed cells.")
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


def _skip_edge_for_holdout(t0: float, t1: float, cfg: dict) -> bool:
    split_mode = str(cfg.get("split_mode", "none"))
    holdout = cfg.get("holdout_time")
    if holdout is None:
        return False
    holdout = float(holdout)
    if split_mode == "teacher_edge_holdout" and float(t0) < holdout < float(t1):
        return True
    if split_mode == "teacher_edge_holdout" and (np.isclose(float(t0), holdout) or np.isclose(float(t1), holdout)):
        return True
    if split_mode == "strict_time_holdout" and (np.isclose(float(t0), holdout) or np.isclose(float(t1), holdout)):
        return True
    return False


def run_toy_sinkhorn_teacher(cfg: dict, label: str = "teacher") -> dict:
    ensure_dir(cfg.get("couplings_dir", "processed/ot_couplings"))
    adata = ad.read_h5ad(cfg["adata_path"])
    z = np.asarray(adata.obsm[cfg.get("latent_key", "X_pca")], dtype=float)
    obs = adata.obs.copy()
    time_key = cfg.get("time_key", "time_numeric")
    cell_type_key = cfg.get("cell_type_key", "lineage")
    allowed_mask = np.ones(adata.n_obs, dtype=bool)
    if "split_role" in obs:
        allowed_mask &= obs["split_role"].astype(str).to_numpy() != "eval_holdout"
    times = sorted(pd.to_numeric(obs.loc[allowed_mask, time_key], errors="coerce").dropna().unique())
    rows = []
    for pair_i, (t0, t1) in enumerate(zip(times[:-1], times[1:])):
        if _skip_edge_for_holdout(float(t0), float(t1), cfg):
            rows.append(
                {
                    "method_label": label,
                    "source_time": t0,
                    "target_time": t1,
                    "file": "",
                    "n_source": 0,
                    "n_target": 0,
                    "solver": "heldout_edge_skipped",
                    "teacher_backend": "toy_sinkhorn_fallback",
                    "transport_cost": np.nan,
                    "mean_entropy": np.nan,
                    "mean_growth": np.nan,
                }
            )
            continue
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
            allowed_mask=allowed_mask,
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
                "teacher_backend": "toy_sinkhorn_fallback",
                "transport_cost": result["mean_cost"],
                "mean_entropy": float(np.mean(result["entropy"])),
                "mean_growth": float(np.mean(result["growth"])),
            }
        )
    index = pd.DataFrame(rows)
    index_path = Path(cfg.get("teacher_index_path", Path(cfg.get("couplings_dir", "processed/ot_couplings")) / f"{label}_coupling_index.csv"))
    index.to_csv(index_path, index=False)
    return {"index_path": str(index_path), "pairs": rows, "teacher_backend": "toy_sinkhorn_fallback"}


def _native_pair_from_subproblem(
    prob,
    full_obs: pd.DataFrame,
    full_obs_names: pd.Index,
    z_full: np.ndarray,
    time_key: str,
    cell_type_key: str,
) -> dict:
    solution = prob.solution
    if solution is None:
        raise RuntimeError("subproblem has no solution")
    plan = np.asarray(solution.transport_matrix, dtype=float)
    if plan.ndim != 2 or plan.size == 0 or not np.isfinite(plan).all() or plan.sum() <= 0:
        raise RuntimeError("native moscot returned invalid transport matrix")
    plan = plan / max(plan.sum(), 1e-300)
    src_names = prob.adata_src.obs_names
    tgt_names = prob.adata_tgt.obs_names
    idx0 = full_obs_names.get_indexer(src_names)
    idx1 = full_obs_names.get_indexer(tgt_names)
    if (idx0 < 0).any() or (idx1 < 0).any():
        raise RuntimeError("could not map native moscot subproblem obs_names to original AnnData indices")
    row_sum = np.maximum(plan.sum(axis=1, keepdims=True), 1e-12)
    transition = plan / row_sum
    y = z_full[idx1]
    bary = transition @ y
    entropy = -np.sum(np.clip(transition, 1e-12, 1.0) * np.log(np.clip(transition, 1e-12, 1.0)), axis=1)
    entropy = entropy / max(np.log(transition.shape[1]), 1.0)
    source_types = full_obs.iloc[idx0][cell_type_key].astype(str).to_numpy()
    target_types = full_obs.iloc[idx1][cell_type_key].astype(str).to_numpy()
    t0 = float(pd.to_numeric(full_obs.iloc[idx0][time_key], errors="coerce").iloc[0])
    t1 = float(pd.to_numeric(full_obs.iloc[idx1][time_key], errors="coerce").iloc[0])
    allowed0 = full_obs[time_key].to_numpy(dtype=float) == t0
    allowed1 = full_obs[time_key].to_numpy(dtype=float) == t1
    full_counts0 = full_obs.loc[allowed0, cell_type_key].astype(str).value_counts()
    full_counts1 = full_obs.loc[allowed1, cell_type_key].astype(str).value_counts()
    growth = np.array([(full_counts1.get(ct, 0.5) + 0.5) / (full_counts0.get(ct, 0.5) + 0.5) for ct in source_types])
    return {
        "source_indices": idx0.astype(np.int64),
        "target_indices": idx1.astype(np.int64),
        "plan": plan.astype(np.float32),
        "transition": transition.astype(np.float32),
        "barycentric": bary.astype(np.float32),
        "entropy": entropy.astype(np.float32),
        "growth": growth.astype(np.float32),
        "source_time": t0,
        "target_time": t1,
        "source_types": source_types,
        "target_types": target_types,
        "method": "native_moscot_temporalproblem",
        "mean_cost": float(getattr(solution, "cost", np.nan)),
        "converged": bool(getattr(solution, "converged", False)),
    }


def run_native_moscot_teacher(cfg: dict, native_status: dict, label: str = "teacher") -> dict:
    """Run native moscot TemporalProblem and extract solved transport matrices.

    This function returns success only after transport matrices have been
    materialized and written to the same NPZ contract used by the fallback
    teacher. Any import, solve or extraction failure is returned as a closed
    failure; callers must then use the toy fallback instead.
    """
    if not bool(cfg.get("use_native_moscot", False)):
        return {"attempted": False, "success": False, "reason": "use_native_moscot=false"}
    if not bool(native_status.get("available", False)):
        return {"attempted": False, "success": False, "reason": f"native import unavailable: {native_status.get('detail', 'unknown')}"}
    try:
        from moscot.problems.time import TemporalProblem

        ensure_dir(cfg.get("couplings_dir", "processed/ot_couplings"))
        adata_full = ad.read_h5ad(cfg["adata_path"])
        z_full = np.asarray(adata_full.obsm[cfg.get("latent_key", "X_pca")], dtype=float)
        obs_full = adata_full.obs.copy()
        time_key = cfg.get("time_key", "time_numeric")
        cell_type_key = cfg.get("cell_type_key", "lineage")
        allowed_mask = np.ones(adata_full.n_obs, dtype=bool)
        if "split_role" in obs_full:
            allowed_mask &= obs_full["split_role"].astype(str).to_numpy() != "eval_holdout"
        adata = adata_full[allowed_mask].copy()
        obs = adata.obs.copy()
        max_cells = int(cfg.get("native_max_cells_per_time", cfg.get("max_cells_per_time", 650)))
        selected = []
        for i, time_value in enumerate(sorted(pd.to_numeric(obs[time_key], errors="coerce").dropna().unique())):
            idx_all = np.where(obs[time_key].to_numpy(dtype=float) == float(time_value))[0]
            idx = idx_all[_stratified_idx(obs.iloc[idx_all], max_cells, int(cfg.get("random_seed", 17)) + i, cell_type_key)]
            selected.extend(idx.tolist())
        selected = np.asarray(sorted(set(selected)), dtype=int)
        adata_native = adata[selected].copy()
        subset = []
        times = sorted(pd.to_numeric(adata_native.obs[time_key], errors="coerce").dropna().unique())
        for t0, t1 in zip(times[:-1], times[1:]):
            if not _skip_edge_for_holdout(float(t0), float(t1), cfg):
                subset.append((float(t0), float(t1)))
        if not subset:
            raise RuntimeError("no native moscot edges remain after holdout filtering")
        problem = TemporalProblem(adata_native)
        problem = problem.prepare(time_key=time_key, joint_attr=cfg.get("latent_key", "X_pca"), policy="explicit", subset=subset)
        problem = problem.solve(
            epsilon=float(cfg.get("epsilon", 0.08)),
            tau_a=float(cfg.get("tau_a", 1.0)),
            tau_b=float(cfg.get("tau_b", 1.0)),
            max_iterations=int(cfg.get("native_max_iterations", 800)),
            batch_size=cfg.get("native_batch_size", None),
            jit=bool(cfg.get("native_jit", False)),
            device=cfg.get("native_device", "cpu"),
        )
        rows = []
        for pair_i, (key, prob) in enumerate(problem.problems.items()):
            result = _native_pair_from_subproblem(prob, obs_full, pd.Index(adata_full.obs_names), z_full, time_key, cell_type_key)
            out = Path(cfg.get("couplings_dir", "processed/ot_couplings")) / f"{label}_native_moscot_t{result['source_time']:g}_to_t{result['target_time']:g}.npz"
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
                    "source_time": result["source_time"],
                    "target_time": result["target_time"],
                    "file": str(out),
                    "n_source": int(result["source_indices"].size),
                    "n_target": int(result["target_indices"].size),
                    "solver": result["method"],
                    "teacher_backend": "native_moscot",
                    "transport_cost": result["mean_cost"],
                    "mean_entropy": float(np.mean(result["entropy"])),
                    "mean_growth": float(np.mean(result["growth"])),
                    "converged": result["converged"],
                    "problem_key": str(key),
                }
            )
        index = pd.DataFrame(rows)
        index_path = Path(cfg.get("teacher_index_path", Path(cfg.get("couplings_dir", "processed/ot_couplings")) / f"{label}_coupling_index.csv"))
        index.to_csv(index_path, index=False)
        return {"attempted": True, "success": True, "index_path": str(index_path), "pairs": rows, "teacher_backend": "native_moscot"}
    except Exception as exc:
        return {"attempted": True, "success": False, "reason": f"{type(exc).__name__}: {exc}"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/ot_teacher.yaml")
    parser.add_argument("--quick-fixture", action="store_true")
    parser.add_argument("--try-native", action="store_true", help="Attempt native moscot extraction for this run; fail closed to toy fallback on any error.")
    args = parser.parse_args()
    cfg = load_config(args.config)
    if args.try_native:
        cfg = dict(cfg)
        cfg["use_native_moscot"] = True
    if args.quick_fixture:
        cfg = dict(cfg)
        cfg["adata_path"] = "processed/quick_fixture/swarmlineage_input.h5ad"
        cfg["teacher_path"] = "processed/quick_fixture/ot_teacher.h5ad"
        cfg["couplings_dir"] = "processed/quick_fixture/ot_couplings"
        cfg["fate_probabilities_path"] = "processed/quick_fixture/ot_fate_probabilities.parquet"
        cfg["teacher_index_path"] = "processed/quick_fixture/ot_couplings/teacher_coupling_index.csv"
        cfg["holdout_time"] = 14.0
        cfg["max_cells_per_time"] = 80
    ensure_dir("logs")
    native = _native_status(int(cfg.get("native_moscot_timeout_seconds", 45)))
    native_result = run_native_moscot_teacher(cfg, native, label="teacher")
    if native_result.get("success", False):
        summary = dict(native_result)
        summary["note"] = "Native moscot TemporalProblem transport matrices were extracted and written to teacher NPZ files."
    else:
        summary = run_toy_sinkhorn_teacher(cfg, label="teacher")
        summary["note"] = "Fallback couplings are explicitly labelled toy_sinkhorn_fallback and must not be reported as moscot results."
    summary["native_moscot_status"] = native
    summary["native_moscot_execution"] = {k: v for k, v in native_result.items() if k != "pairs"}
    summary["native_moscot_used"] = bool(native_result.get("success", False))
    write_json(Path(cfg.get("couplings_dir", "processed/ot_couplings")) / "moscot_run_summary.json", summary)
    print(json.dumps({"teacher_pairs": len(summary["pairs"]), "teacher_backend": summary["teacher_backend"], "native": native}, indent=2))


if __name__ == "__main__":
    main()
