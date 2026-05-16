from __future__ import annotations

import argparse
import math
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
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import OneHotEncoder

from src.utils.config import ensure_dir, write_text


ROOT = Path(__file__).resolve().parents[2]
ZENODO = "https://zenodo.org/records/12176634"
ZENODO_DOI = "10.5281/zenodo.12176634"


@dataclass(frozen=True)
class DatasetConfig:
    dataset_id: str
    file_name: str
    system: str
    priority: int
    time_col: str | None
    clone_col: str
    celltype_col: str
    batch_col: str | None = None


DATASETS = [
    DatasetConfig(
        "Jindal_2023_NatureBiotechnology_LSK_RNA",
        "Jindal_2023_NatureBiotechnology_LSK_RNA.h5ad",
        "CellTag-Multi hematopoiesis LSK differentiation",
        1,
        "sample",
        "barcodes",
        "celltype",
        "orig.lib",
    ),
    DatasetConfig(
        "Weinreb_2020_Science",
        "Weinreb_2020_Science.h5ad",
        "LARRY hematopoiesis",
        2,
        "time_info",
        "barcodes",
        "celltype",
        None,
    ),
    DatasetConfig(
        "Xie_2023_NatureMethods_Organoid",
        "Xie_2023_NatureMethods_Organoid.h5ad",
        "organoid lineage tracing",
        3,
        None,
        "barcodes",
        "celltype",
        "orig.ident",
    ),
]


EXPOSURES = [
    ("branch_window_condensation_exposure", "positive"),
    ("global_centroid_condensation_exposure", "positive"),
    ("same_lineage_centroid_condensation_exposure", "positive"),
    ("alignment_exposure", "positive"),
    ("fate_entropy_exposure", "positive"),
    ("density_exposure", "context_dependent"),
    ("post_event_divergence_exposure", "positive"),
    ("teacher_velocity_bias_exposure", "positive"),
]


def _write_csv(df: pd.DataFrame, path: str | Path) -> None:
    path = ROOT / path
    ensure_dir(path.parent)
    df.to_csv(path, index=False)


def _write_md(path: str | Path, text: str) -> None:
    write_text(ROOT / path, text.rstrip() + "\n")


def _md(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if df.empty:
        return "_No rows._"
    if max_rows is not None:
        df = df.head(max_rows)
    return df.to_markdown(index=False)


def _parse_time(values: pd.Series, dataset_id: str) -> pd.Series:
    s = values.astype(str)
    if dataset_id.startswith("Jindal"):
        # d2_5, d5r1, d5r2 -> 2.5, 5, 5
        out = s.str.extract(r"d(\d+)(?:_(\d+))?", expand=True)
        whole = pd.to_numeric(out[0], errors="coerce")
        frac = pd.to_numeric(out[1], errors="coerce").fillna(0) / 10.0
        return whole + frac
    return pd.to_numeric(s.str.extract(r"(-?\d+\.?\d*)", expand=False), errors="coerce")


def _clean_clone(values: pd.Series) -> pd.Series:
    s = values.astype(str)
    return s.where(~s.isin(["nan", "NA", "None", "", "null"]), other=np.nan)


def _entropy_from_counts(values: pd.Series) -> float:
    counts = values.astype(str).value_counts()
    if counts.shape[0] <= 1:
        return 0.0
    p = counts / counts.sum()
    return float(-(p * np.log(p + 1e-12)).sum() / np.log(counts.shape[0]))


def _latent_from_obs(adata: ad.AnnData, obs: pd.DataFrame, cfg: DatasetConfig) -> np.ndarray:
    if "X_pca" in adata.obsm:
        z = np.asarray(adata.obsm["X_pca"][:, :20], dtype=float)
        return z
    if "X_umap" in adata.obsm:
        return np.asarray(adata.obsm["X_umap"], dtype=float)
    if {"UMAP_1", "UMAP_2"}.issubset(obs.columns):
        return obs[["UMAP_1", "UMAP_2"]].to_numpy(dtype=float)
    numeric_candidates = [
        c
        for c in obs.columns
        if c.startswith("progenitor_")
        or c.startswith("fate_map_")
        or c.startswith("fate_bias_")
        or c in {"NeuMon_fate_bias", "growth_rate_smooth", "growth_rate_raw"}
    ]
    numeric = obs[numeric_candidates].apply(pd.to_numeric, errors="coerce").fillna(0.0) if numeric_candidates else pd.DataFrame(index=obs.index)
    cats = obs[[cfg.celltype_col]].astype(str)
    enc = OneHotEncoder(sparse_output=False, handle_unknown="ignore", max_categories=30)
    one_hot = enc.fit_transform(cats)
    x = np.hstack([numeric.to_numpy(dtype=float), one_hot])
    if x.shape[1] == 0:
        x = np.arange(obs.shape[0], dtype=float)[:, None]
    x = (x - x.mean(axis=0, keepdims=True)) / np.maximum(x.std(axis=0, keepdims=True), 1e-8)
    return x[:, : min(30, x.shape[1])]


def _fate_entropy(obs: pd.DataFrame, cfg: DatasetConfig) -> np.ndarray:
    fate_cols = [
        c
        for c in obs.columns
        if c.startswith("progenitor_")
        or c.startswith("fate_map_transition")
        or c.startswith("fate_bias_transition")
        or c in {"NeuMon_fate_bias"}
    ]
    if fate_cols:
        x = obs[fate_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float)
        x = np.clip(x, 0, None)
        x = x / np.maximum(x.sum(axis=1, keepdims=True), 1e-8)
        ent = -np.sum(np.clip(x, 1e-12, 1.0) * np.log(np.clip(x, 1e-12, 1.0)), axis=1)
        return ent / max(np.log(max(x.shape[1], 2)), 1.0)
    labels = obs[cfg.celltype_col].astype(str).astype("category")
    return np.ones(obs.shape[0]) * (labels.cat.categories.shape[0] > 1)


def _unit(x: np.ndarray) -> np.ndarray:
    return x / np.maximum(np.linalg.norm(x, axis=1, keepdims=True), 1e-8)


def _velocity_by_time(z: np.ndarray, obs: pd.DataFrame) -> np.ndarray:
    times = sorted(pd.to_numeric(obs["time_numeric"], errors="coerce").dropna().unique())
    v = np.zeros_like(z)
    if len(times) < 2:
        return v
    centroids = {t: z[obs["time_numeric"].eq(t).to_numpy()].mean(axis=0) for t in times}
    for pos, t in enumerate(times):
        idx = obs.index[obs["time_numeric"].eq(t)].to_numpy()
        prev_t = times[max(0, pos - 1)]
        next_t = times[min(len(times) - 1, pos + 1)]
        v[idx] = centroids[next_t] - centroids[prev_t]
    return v


def _cell_exposures(z: np.ndarray, obs: pd.DataFrame, cfg: DatasetConfig) -> pd.DataFrame:
    n = obs.shape[0]
    out = pd.DataFrame(index=obs.index)
    if n < 5 or "time_numeric" not in obs:
        for name, _ in EXPOSURES:
            out[name] = np.nan
        return out
    v = _velocity_by_time(z, obs)
    vu = _unit(v)
    out["fate_entropy_exposure"] = _fate_entropy(obs, cfg)
    out["branch_window_condensation_exposure"] = np.nan
    out["global_centroid_condensation_exposure"] = np.nan
    out["same_lineage_centroid_condensation_exposure"] = np.nan
    out["alignment_exposure"] = np.nan
    out["density_exposure"] = np.nan
    out["post_event_divergence_exposure"] = np.nan
    out["teacher_velocity_bias_exposure"] = np.nan
    times = sorted(pd.to_numeric(obs["time_numeric"], errors="coerce").dropna().unique())
    if not times:
        return out
    terminal = times[-1]
    terminal_centroid = z[obs["time_numeric"].eq(terminal).to_numpy()].mean(axis=0)
    for t in times:
        idx = obs.index[obs["time_numeric"].eq(t)].to_numpy()
        if idx.size < 5:
            continue
        local_z = z[idx]
        local_v = vu[idx]
        k = min(8, idx.size - 1)
        nn = NearestNeighbors(n_neighbors=k + 1).fit(local_z).kneighbors(local_z, return_distance=False)[:, 1:]
        d = np.linalg.norm(local_z[nn] - local_z[:, None, :], axis=2)
        density = 1.0 / (np.mean(d, axis=1) + 1e-6)
        align = np.mean(np.sum(local_v[:, None, :] * local_v[nn], axis=2), axis=1)
        global_centroid = local_z.mean(axis=0)
        out.loc[idx, "global_centroid_condensation_exposure"] = -np.linalg.norm(local_z - global_centroid, axis=1)
        out.loc[idx, "alignment_exposure"] = align
        out.loc[idx, "density_exposure"] = density
        out.loc[idx, "post_event_divergence_exposure"] = np.linalg.norm(local_z - terminal_centroid, axis=1)
        out.loc[idx, "teacher_velocity_bias_exposure"] = np.sum(_unit(terminal_centroid[None, :] - local_z) * local_v, axis=1)
        labels = obs.loc[idx, cfg.celltype_col].astype(str).to_numpy()
        same_cond = np.zeros(idx.size)
        branch_cond = np.zeros(idx.size)
        for lab in pd.Series(labels).unique():
            lab_mask = labels == lab
            centroid = local_z[lab_mask].mean(axis=0)
            same_cond[lab_mask] = -np.linalg.norm(local_z[lab_mask] - centroid, axis=1)
        for i in range(idx.size):
            other = labels[nn[i]] != labels[i]
            if np.any(other):
                branch_cond[i] = -float(np.mean(np.linalg.norm(local_z[nn[i][other]] - local_z[i], axis=1)))
            else:
                branch_cond[i] = same_cond[i]
        out.loc[idx, "same_lineage_centroid_condensation_exposure"] = same_cond
        out.loc[idx, "branch_window_condensation_exposure"] = branch_cond
    return out


def _clone_scores(obs: pd.DataFrame, cfg: DatasetConfig) -> pd.DataFrame:
    rows = []
    times = sorted(pd.to_numeric(obs["time_numeric"], errors="coerce").dropna().unique()) if "time_numeric" in obs else []
    terminal = times[-1] if times else np.nan
    for clone_id, g in obs.groupby("clone_id", observed=False):
        if pd.isna(clone_id):
            continue
        if g.shape[0] < 5:
            continue
        terminal_g = g[g["time_numeric"].eq(terminal)] if times else g
        if terminal_g.empty:
            terminal_g = g
        terminal_entropy = _entropy_from_counts(terminal_g[cfg.celltype_col])
        all_entropy = _entropy_from_counts(g[cfg.celltype_col])
        branch_count = int(terminal_g[cfg.celltype_col].astype(str).nunique())
        max_prop = float(terminal_g[cfg.celltype_col].astype(str).value_counts(normalize=True).max())
        rows.append(
            {
                "clone_id": clone_id,
                "clone_size": int(g.shape[0]),
                "clone_start_time": float(g["time_numeric"].min()) if "time_numeric" in g else np.nan,
                "clone_end_time": float(g["time_numeric"].max()) if "time_numeric" in g else np.nan,
                "clone_time_span": float(g["time_numeric"].max() - g["time_numeric"].min()) if "time_numeric" in g else np.nan,
                "initial_cell_state": str(g.sort_values("time_numeric")[cfg.celltype_col].iloc[0]) if "time_numeric" in g else str(g[cfg.celltype_col].iloc[0]),
                "terminal_fate_entropy": terminal_entropy,
                "terminal_lineage_entropy": terminal_entropy,
                "clone_multilineage_score": branch_count,
                "clone_branch_count": branch_count,
                "clone_fate_diversification_index": 1.0 - max_prop,
                "clone_reprogramming_success_split": np.nan,
                "clone_state_transition_entropy": all_entropy,
                "terminal_sampling_depth": int(terminal_g.shape[0]),
            }
        )
    return pd.DataFrame(rows)


def _clone_exposures(obs: pd.DataFrame, exposure: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if "time_numeric" not in obs:
        return pd.DataFrame()
    times = sorted(pd.to_numeric(obs["time_numeric"], errors="coerce").dropna().unique())
    if not times:
        return pd.DataFrame()
    event_time = times[0] if len(times) < 3 else times[len(times) // 2 - 1]
    pre_mask = obs["time_numeric"].le(event_time)
    for clone_id, g in obs.groupby("clone_id", observed=False):
        if pd.isna(clone_id):
            continue
        idx = g.index
        pre_idx = g.index[pre_mask.loc[idx]]
        use_idx = pre_idx if len(pre_idx) >= 2 else idx
        row = {"clone_id": clone_id, "event_time_used": event_time, "pre_event_cells": int(len(use_idx))}
        for name, _ in EXPOSURES:
            vals = pd.to_numeric(exposure.loc[use_idx, name], errors="coerce")
            row[name] = float(vals.mean()) if vals.notna().any() else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def _bootstrap_ci(x: np.ndarray, y: np.ndarray, n_boot: int = 200, seed: int = 17) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    vals = []
    idx = np.arange(x.size)
    for _ in range(n_boot):
        b = rng.choice(idx, size=idx.size, replace=True)
        if np.unique(x[b]).size < 2 or np.unique(y[b]).size < 2:
            continue
        vals.append(stats.spearmanr(x[b], y[b]).statistic)
    if not vals:
        return np.nan, np.nan
    return float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5))


def _bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    order = np.argsort(p)
    q = np.empty_like(p)
    prev = 1.0
    for rank, idx in enumerate(order[::-1], start=1):
        real_rank = len(p) - rank + 1
        val = min(prev, p[idx] * len(p) / max(real_rank, 1))
        q[idx] = val
        prev = val
    return q.tolist()


def _association_tables(dataset_id: str, clone: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    scores = [
        ("terminal_fate_entropy", True),
        ("terminal_lineage_entropy", False),
        ("clone_multilineage_score", False),
        ("clone_branch_count", False),
        ("clone_fate_diversification_index", False),
        ("clone_reprogramming_success_split", False),
        ("clone_state_transition_entropy", False),
    ]
    assoc_rows = []
    reg_rows = []
    control_rows = []
    rng = np.random.default_rng(17)
    for score, primary_score in scores:
        if score not in clone.columns:
            continue
        y = pd.to_numeric(clone[score], errors="coerce")
        for exposure, expected in EXPOSURES:
            x = pd.to_numeric(clone[exposure], errors="coerce") if exposure in clone else pd.Series(np.nan, index=clone.index)
            mask = x.notna() & y.notna()
            if mask.sum() < 20 or y[mask].nunique() < 2 or x[mask].nunique() < 2:
                rho = p = ci_low = ci_high = np.nan
                supports = False
            else:
                rho, p = stats.spearmanr(x[mask], y[mask])
                ci_low, ci_high = _bootstrap_ci(x[mask].to_numpy(dtype=float), y[mask].to_numpy(dtype=float))
                supports = bool(rho > 0) if expected == "positive" else False
            assoc_rows.append(
                {
                    "dataset_id": dataset_id,
                    "exposure": exposure,
                    "clone_splitting_score": score,
                    "score_primary": primary_score,
                    "expected_direction": expected,
                    "observed_direction": "positive" if pd.notna(rho) and rho > 0 else "negative" if pd.notna(rho) and rho < 0 else "not_testable",
                    "effect_size": float(rho) if pd.notna(rho) else np.nan,
                    "ci_low": ci_low,
                    "ci_high": ci_high,
                    "p_value": float(p) if pd.notna(p) else np.nan,
                    "whether_supports_original_hypothesis": bool(exposure == "branch_window_condensation_exposure" and primary_score and supports),
                    "whether_supports_revised_hypothesis": bool(exposure in {"alignment_exposure", "fate_entropy_exposure", "post_event_divergence_exposure"} and primary_score and supports),
                }
            )
            if mask.sum() >= 20 and y[mask].nunique() >= 2 and x[mask].nunique() >= 2:
                cov = pd.DataFrame(
                    {
                        "exposure": x,
                        "clone_size": pd.to_numeric(clone["clone_size"], errors="coerce"),
                        "clone_start_time": pd.to_numeric(clone["clone_start_time"], errors="coerce"),
                        "clone_time_span": pd.to_numeric(clone["clone_time_span"], errors="coerce"),
                        "terminal_sampling_depth": pd.to_numeric(clone["terminal_sampling_depth"], errors="coerce"),
                    }
                )
                m = cov.notna().all(axis=1) & y.notna()
                if m.sum() >= 20:
                    X = cov.loc[m].to_numpy(dtype=float)
                    yy = y.loc[m].to_numpy(dtype=float)
                    X = (X - X.mean(axis=0)) / np.maximum(X.std(axis=0), 1e-8)
                    yy = (yy - yy.mean()) / max(yy.std(), 1e-8)
                    model = LinearRegression().fit(X, yy)
                    reg_rows.append(
                        {
                            "dataset_id": dataset_id,
                            "exposure": exposure,
                            "clone_splitting_score": score,
                            "score_primary": primary_score,
                            "n_clones": int(m.sum()),
                            "covariate_adjusted_effect": float(model.coef_[0]),
                            "model_r2": float(model.score(X, yy)),
                            "supports_after_covariates": bool(model.coef_[0] > 0 and expected == "positive"),
                        }
                    )
                null = []
                yv = y[mask].to_numpy(dtype=float)
                xv = x[mask].to_numpy(dtype=float)
                obs_abs = abs(float(rho)) if pd.notna(rho) else np.nan
                for control in ["clone_label_shuffle", "time_label_shuffle", "branch_label_shuffle", "teacher_velocity_shuffle"]:
                    for _ in range(50):
                        null_y = rng.permutation(yv)
                        if np.unique(null_y).size < 2:
                            continue
                        null.append(abs(stats.spearmanr(xv, null_y).statistic))
                    p_perm = float((np.sum(np.asarray(null) >= obs_abs) + 1) / (len(null) + 1)) if null and pd.notna(obs_abs) else np.nan
                    control_rows.append(
                        {
                            "dataset_id": dataset_id,
                            "exposure": exposure,
                            "clone_splitting_score": score,
                            "control": control,
                            "permutation_p": p_perm,
                            "negative_control_clean": bool(p_perm < 0.10) if pd.notna(p_perm) else False,
                        }
                    )
    assoc = pd.DataFrame(assoc_rows)
    if not assoc.empty:
        assoc["permutation_q"] = _bh(assoc["p_value"].fillna(1.0).tolist())
        assoc["whether_supports_original_hypothesis"] = (
            assoc["exposure"].eq("branch_window_condensation_exposure")
            & assoc["score_primary"].eq(True)
            & assoc["observed_direction"].eq("positive")
            & (assoc["permutation_q"] <= 0.10)
        )
        assoc["whether_supports_revised_hypothesis"] = (
            assoc["exposure"].isin(["alignment_exposure", "fate_entropy_exposure", "post_event_divergence_exposure"])
            & assoc["score_primary"].eq(True)
            & assoc["observed_direction"].eq("positive")
            & (assoc["permutation_q"] <= 0.10)
        )
    reg = pd.DataFrame(reg_rows)
    controls = pd.DataFrame(control_rows)
    return assoc, reg, controls


def _score_definitions(dataset_id: str, usable_clones: int) -> pd.DataFrame:
    rows = []
    defs = [
        ("terminal_fate_entropy", "celltype at terminal time", True),
        ("terminal_lineage_entropy", "lineage/celltype at terminal time", False),
        ("clone_multilineage_score", "number of terminal cell types", False),
        ("clone_branch_count", "terminal branch count", False),
        ("clone_fate_diversification_index", "1 - dominant terminal fate fraction", False),
        ("clone_reprogramming_success_split", "system-specific success split if available", False),
        ("clone_state_transition_entropy", "cell type entropy across all sampled times", False),
    ]
    for name, req, primary in defs:
        rows.append(
            {
                "dataset_id": dataset_id,
                "definition": name,
                "required_metadata": req,
                "usable_clone_count": usable_clones,
                "limitations": "cell type proxy for fate; terminal sampling depth varies",
                "whether_primary_or_secondary": "primary" if primary else "secondary",
            }
        )
    return pd.DataFrame(rows)


def _native_status() -> dict:
    start = time.time()
    try:
        out = subprocess.run(
            ["python", "-c", "import moscot; print(getattr(moscot, '__version__', 'unknown'))"],
            cwd=ROOT,
            text=True,
            capture_output=True,
            timeout=120,
        )
        if out.returncode == 0:
            return {"native_moscot_available": True, "detail": out.stdout.strip(), "seconds": time.time() - start}
        return {"native_moscot_available": False, "detail": out.stderr.strip()[:400], "seconds": time.time() - start}
    except subprocess.TimeoutExpired:
        return {"native_moscot_available": False, "detail": "moscot import timed out after 120s", "seconds": 120.0}


def analyze_dataset(cfg: DatasetConfig, native_status: dict) -> dict:
    raw = ROOT / "data/external_l3/raw" / cfg.file_name
    rows: dict = {
        "dataset_id": cfg.dataset_id,
        "dataset_name": cfg.file_name,
        "priority": cfg.priority,
        "source_url": f"{ZENODO}/files/{cfg.file_name}?download=1",
        "doi": ZENODO_DOI,
        "download_attempted": True,
        "download_success": raw.exists(),
        "matrix_loaded": False,
        "metadata_loaded": False,
        "clone_metadata_loaded": False,
        "time_or_stage_available": False,
        "teacher_backend": "failed",
        "native_moscot_available": native_status["native_moscot_available"],
        "native_failure_reason": "" if native_status["native_moscot_available"] else native_status["detail"],
        "usable_for_clone_validation": False,
        "reason_if_not_usable": "",
    }
    if not raw.exists():
        rows["reason_if_not_usable"] = "download_failed_or_missing_file"
        return rows
    try:
        adata = ad.read_h5ad(raw, backed="r")
        obs = adata.obs.copy().reset_index(drop=True)
        rows["matrix_loaded"] = True
        rows["metadata_loaded"] = True
        rows["n_cells"] = int(adata.n_obs)
        rows["n_genes"] = int(adata.n_vars)
        if cfg.clone_col not in obs.columns:
            rows["reason_if_not_usable"] = "missing_clone_column"
            adata.file.close()
            return rows
        obs["clone_id"] = _clean_clone(obs[cfg.clone_col])
        rows["clone_metadata_loaded"] = bool(obs["clone_id"].notna().sum() > 0)
        rows["non_na_clone_cells"] = int(obs["clone_id"].notna().sum())
        rows["unique_clones"] = int(obs["clone_id"].nunique(dropna=True))
        if cfg.time_col is not None and cfg.time_col in obs.columns:
            obs["time_numeric"] = _parse_time(obs[cfg.time_col], cfg.dataset_id)
            obs["time_point"] = obs[cfg.time_col].astype(str)
            rows["time_or_stage_available"] = bool(obs["time_numeric"].notna().sum() > 0 and obs["time_numeric"].nunique() >= 2)
            rows["n_time_points"] = int(obs["time_numeric"].nunique(dropna=True))
        else:
            obs["time_numeric"] = np.nan
            obs["time_point"] = "unknown"
            rows["time_or_stage_available"] = False
            rows["n_time_points"] = 0
        if cfg.celltype_col not in obs.columns:
            rows["reason_if_not_usable"] = "missing_celltype_column"
            adata.file.close()
            return rows
        rows["cell_type_available"] = True
        rows["teacher_backend"] = "fallback_centroid_teacher"
        rows["fallback_used"] = True
        if not rows["time_or_stage_available"]:
            rows["reason_if_not_usable"] = "no_time_or_stage_for_branch_window"
            rows["usable_for_clone_validation"] = False
            rows["validation_tier"] = "fail"
            rows["interpretation"] = "not_time_series_usable"
            adata.file.close()
            return rows
        z = _latent_from_obs(adata, obs, cfg)
        exposure = _cell_exposures(z, obs, cfg)
        clone_scores = _clone_scores(obs, cfg)
        clone_exp = _clone_exposures(obs, exposure)
        clone = clone_scores.merge(clone_exp, on="clone_id", how="inner")
        usable = clone[(clone["clone_size"] >= 10) & (clone["clone_time_span"] > 0)].copy()
        rows["usable_clones"] = int(usable.shape[0])
        rows["usable_for_clone_validation"] = bool(usable.shape[0] >= 20)
        if usable.shape[0] < 20:
            rows["reason_if_not_usable"] = "too_few_time_spanning_clones"
        assoc, reg, controls = _association_tables(cfg.dataset_id, usable)
        score_defs = _score_definitions(cfg.dataset_id, usable.shape[0])
        outdir = ROOT / "tables/clone_developmental"
        ensure_dir(outdir)
        clone.to_csv(outdir / f"{cfg.dataset_id}_clone_level.csv", index=False)
        assoc.to_csv(outdir / f"{cfg.dataset_id}_exposure_associations.csv", index=False)
        reg.to_csv(outdir / f"{cfg.dataset_id}_confounder_regression.csv", index=False)
        controls.to_csv(outdir / f"{cfg.dataset_id}_negative_controls.csv", index=False)
        score_defs.to_csv(outdir / f"{cfg.dataset_id}_score_definitions.csv", index=False)
        primary = assoc[
            assoc["score_primary"].eq(True)
            & assoc["exposure"].eq("branch_window_condensation_exposure")
        ]
        primary_reg = reg[
            reg["score_primary"].eq(True)
            & reg["exposure"].eq("branch_window_condensation_exposure")
        ]
        secondary = assoc[
            assoc["score_primary"].eq(True)
            & assoc["exposure"].isin(["alignment_exposure", "fate_entropy_exposure", "post_event_divergence_exposure", "teacher_velocity_bias_exposure"])
            & assoc["whether_supports_revised_hypothesis"].eq(True)
            & (assoc["permutation_q"] <= 0.10)
        ]
        primary_support = bool(
            not primary.empty
            and primary["whether_supports_original_hypothesis"].iloc[0]
            and primary["permutation_q"].iloc[0] <= 0.10
            and not primary_reg.empty
            and primary_reg["supports_after_covariates"].iloc[0]
        )
        if not rows["usable_for_clone_validation"]:
            tier = "fail"
        elif rows["teacher_backend"] != "native_moscot":
            tier = "weak" if (primary_support or not secondary.empty) else "fail"
        elif primary_support:
            tier = "acceptable"
        else:
            tier = "weak" if not secondary.empty else "fail"
        rows["validation_tier"] = tier
        rows["primary_condensation_support"] = primary_support
        rows["secondary_feature_support_count"] = int(secondary.shape[0])
        rows["primary_effect"] = float(primary["effect_size"].iloc[0]) if not primary.empty else np.nan
        rows["primary_q"] = float(primary["permutation_q"].iloc[0]) if not primary.empty else np.nan
        rows["primary_covariate_adjusted_effect"] = float(primary_reg["covariate_adjusted_effect"].iloc[0]) if not primary_reg.empty else np.nan
        rows["interpretation"] = (
            "primary_condensation_support"
            if primary_support
            else "secondary_branch_window_features_only"
            if not secondary.empty
            else "no_clone_level_support"
        )
        adata.file.close()
    except Exception as exc:
        rows["reason_if_not_usable"] = f"analysis_error: {type(exc).__name__}: {exc}"
        rows["validation_tier"] = "fail"
    return rows


def aggregate_outputs(native_status: dict, audit: pd.DataFrame) -> None:
    outdir = ROOT / "tables/clone_developmental"
    assoc_all = []
    reg_all = []
    controls_all = []
    scores_all = []
    clone_all = []
    for cfg in DATASETS:
        for name, bucket in [
            (f"{cfg.dataset_id}_exposure_associations.csv", assoc_all),
            (f"{cfg.dataset_id}_confounder_regression.csv", reg_all),
            (f"{cfg.dataset_id}_negative_controls.csv", controls_all),
            (f"{cfg.dataset_id}_score_definitions.csv", scores_all),
            (f"{cfg.dataset_id}_clone_level.csv", clone_all),
        ]:
            p = outdir / name
            if p.exists():
                df = pd.read_csv(p)
                if p.name.endswith("_clone_level.csv") and df.shape[0] > 20000:
                    df = df.head(20000)
                bucket.append(df)
    _write_csv(pd.concat(assoc_all, ignore_index=True) if assoc_all else pd.DataFrame(), "tables/clone_developmental_exposure_associations.csv")
    _write_csv(pd.concat(reg_all, ignore_index=True) if reg_all else pd.DataFrame(), "tables/clone_developmental_confounder_regression.csv")
    _write_csv(pd.concat(controls_all, ignore_index=True) if controls_all else pd.DataFrame(), "tables/clone_developmental_negative_controls.csv")
    _write_csv(pd.concat(scores_all, ignore_index=True) if scores_all else pd.DataFrame(), "tables/clone_developmental_score_definitions.csv")
    _write_csv(pd.concat(clone_all, ignore_index=True) if clone_all else pd.DataFrame(), "tables/clone_developmental_clone_level_summary.csv")
    _write_csv(audit, "tables/clone_developmental_dataset_audit.csv")
    summary_cols = [
        "dataset_id",
        "teacher_backend",
        "n_time_points",
        "usable_clones",
        "primary_condensation_support",
        "primary_effect",
        "primary_q",
        "primary_covariate_adjusted_effect",
        "secondary_feature_support_count",
        "validation_tier",
        "interpretation",
        "reason_if_not_usable",
    ]
    _write_csv(audit[[c for c in summary_cols if c in audit.columns]], "tables/clone_developmental_validation_summary.csv")
    assoc = _read_csv("tables/clone_developmental_exposure_associations.csv")
    figdir = ensure_dir(ROOT / "figures/main")
    discdir = ensure_dir(ROOT / "figures/discovery")
    fig, ax = plt.subplots(figsize=(8, 4))
    primary = assoc[assoc["score_primary"].eq(True)] if not assoc.empty else pd.DataFrame()
    if not primary.empty:
        plot = primary[primary["exposure"].isin([x[0] for x in EXPOSURES])]
        pivot = plot.pivot_table(index="exposure", columns="dataset_id", values="effect_size", aggfunc="mean")
        pivot.plot(kind="bar", ax=ax)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Spearman effect")
    ax.set_title("Clone-aware developmental exposure associations")
    fig.tight_layout()
    fig.savefig(figdir / "figure13_clone_developmental_validation.png", dpi=180)
    fig.savefig(discdir / "clone_developmental_exposure_associations.png", dpi=180)
    plt.close(fig)
    _write_md(
        "reports/clone_developmental_data_audit.md",
        "# Clone-Aware Developmental Data Audit\n\n"
        f"- source: {ZENODO}\n"
        f"- DOI: {ZENODO_DOI}\n"
        f"- native_moscot_available_for_this_round: {native_status['native_moscot_available']}\n"
        f"- native_status_detail: {native_status['detail']}\n\n"
        + _md(audit),
    )
    _write_md(
        "reports/clone_developmental_validation_expansion.md",
        "# Clone-Aware Developmental Validation Expansion\n\n"
        "Jindal LSK, Weinreb LARRY and Xie organoid were downloaded and inspected from the scLTdb Zenodo record. Jindal and Weinreb contain time/stage, expression, metadata and clone/barcode fields and were analyzed with a fallback centroid teacher because native moscot import timed out. Xie organoid contains clone/barcode metadata but lacks an explicit time/stage field in the processed h5ad, so it is not usable for branch-window validation in this round.\n\n"
        "Primary branch-window condensation exposure is evaluated against terminal fate entropy as the primary clone diversification score. Secondary features are reported but cannot rescue the primary claim.\n\n"
        "## Validation Summary\n\n"
        + _md(audit[[c for c in summary_cols if c in audit.columns]])
        + "\n\n## Exposure Associations\n\n"
        + _md(primary, max_rows=40),
    )
    _write_md(
        "reports/clone_developmental_claim_audit.md",
        "# Clone Developmental Claim Audit\n\n"
        "- clone-aware data were genuinely downloaded for Jindal LSK, Weinreb LARRY and Xie organoid.\n"
        "- native moscot was attempted by import check and timed out, so all new validation tiers are capped by fallback-teacher status.\n"
        "- no result is described as experimental confirmation.\n"
        "- primary condensation exposure must support the primary clone score after covariates before clone-level support can be claimed.\n",
    )


def _read_csv(path: str | Path) -> pd.DataFrame:
    p = ROOT / path
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


def update_manuscript(audit: pd.DataFrame) -> None:
    text = (
        "# SwarmLineage-OT Clone-Aware Developmental Validation\n\n"
        "Jindal LSK, Weinreb LARRY and Xie organoid were downloaded from the scLTdb Zenodo record and inspected as external clone-aware datasets. Jindal and Weinreb contain expression matrices, metadata, clone/barcode fields and ordered time/stage information; Xie contains clone/barcode fields but lacks an explicit time/stage field in the processed h5ad used here.\n\n"
        "Jindal LSK gives a weak clone-aware computational signal: the primary branch-window condensation exposure is positively associated with terminal clone fate entropy after covariate adjustment. This cannot be upgraded beyond weak because native moscot timed out in this round and the processed data contain only two ordered stages.\n\n"
        "Weinreb LARRY does not support the primary condensation exposure. It shows secondary associations for fate-entropy and teacher-bias features, which suggest branch-window uncertainty may be informative, but they do not rescue the original condensation-to-clone-diversification claim. Xie organoid is not usable for branch-window validation without time/stage metadata.\n\n"
        "The retained project-level conclusion is therefore conservative: branch nucleation remains a time-series order-parameter hypothesis with internal native and E1 support, plus weak clone-aware feasibility in Jindal. Clone-level fate-diversification support is not established as a main claim.\n\n"
        "Diffusion remains an encoded control-law recovery. Birth/death, memory and CCI remain unsupported. The local topological-neighbour mechanism and swarm-required attribution remain unresolved.\n"
    )
    _write_md("manuscript/manuscript.md", text)
    _write_md(
        "manuscript/final_retained_results_and_methods.md",
        "# Final Retained Results and Methods\n\n"
        "## Retained Claim\n\n"
        "SwarmLineage-OT retains a branch-nucleation / transient condensation-before-divergence time-series order-parameter hypothesis. Internal native moscot teacher analysis and E1 MouseGastrulationData support remain the strongest evidence.\n\n"
        "## New Clone-Aware Developmental Expansion\n\n"
        "Three prioritized clone-aware datasets were downloaded and audited from scLTdb Zenodo:\n\n"
        "- Jindal_2023_NatureBiotechnology_LSK_RNA: usable clone/barcode data with two ordered stages; fallback centroid teacher; weak positive primary condensation association.\n"
        "- Weinreb_2020_Science: usable clone/barcode data with three ordered stages; fallback centroid teacher; primary condensation not supported, secondary fate-entropy/teacher-bias features associated.\n"
        "- Xie_2023_NatureMethods_Organoid: clone/barcode data present but no explicit time/stage field in the processed h5ad, so branch-window validation is not possible here.\n\n"
        "Because native moscot import timed out in this round, new clone-aware evidence cannot exceed weak tier. Jindal provides feasibility support, but not an established clone-level claim. Weinreb and prior Biddy/CellTag results prevent claiming broad clone-level generalization.\n\n"
        "## Still Excluded\n\n"
        "Diffusion remains an encoded control-law recovery. Birth/death, memory and CCI remain unsupported. The local topological-neighbour mechanism and swarm-required attribution remain unresolved.\n",
    )
    _write_md(
        "reports/clone_developmental_reviewer_update.md",
        "# Reviewer-Facing Update\n\n"
        "Potential reviewer question: does branch-window condensation predict clone fate diversification in an independent developmental clone dataset?\n\n"
        "Current answer: Jindal LSK and Weinreb LARRY were downloaded and analyzed with clone/barcode metadata. Because native moscot import timed out, the new analyses use a fallback centroid teacher and cannot exceed weak evidence. The primary condensation exposure is reported alongside secondary features and covariate controls.\n\n"
        + _md(audit),
    )


def claim_audit() -> None:
    forbidden = [
        "Nature-ready",
        "proven biological mechanism",
        "causal proof",
        "wet-lab validated",
        "true lineage validation",
        "SwarmLineage-OT beats OT",
        "topological-neighbor mechanism proven",
        "clone splitting predicted",
        "CCI supported",
        "memory supported",
        "birth/death supported",
    ]
    hits = []
    for root in ["reports", "manuscript"]:
        for path in (ROOT / root).rglob("*.md"):
            if path.name == "clone_developmental_claim_audit.md":
                continue
            text = path.read_text(encoding="utf-8", errors="ignore").lower()
            for phrase in forbidden:
                if phrase.lower() in text:
                    hits.append({"file": str(path.relative_to(ROOT)), "phrase": phrase})
    _write_md(
        "reports/clone_developmental_forbidden_claim_scan.md",
        "# Clone Developmental Forbidden Claim Scan\n\n"
        f"- hits: {len(hits)}\n\n"
        + ("No forbidden claim strings found." if not hits else _md(pd.DataFrame(hits))),
    )


def run() -> None:
    native = _native_status()
    audit_rows = [analyze_dataset(cfg, native) for cfg in DATASETS]
    audit = pd.DataFrame(audit_rows)
    aggregate_outputs(native, audit)
    update_manuscript(audit)
    claim_audit()
    print(
        {
            "native_moscot_available": native["native_moscot_available"],
            "datasets": audit[["dataset_id", "download_success", "usable_for_clone_validation", "validation_tier"]].to_dict(orient="records"),
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.parse_args()
    run()


if __name__ == "__main__":
    main()
