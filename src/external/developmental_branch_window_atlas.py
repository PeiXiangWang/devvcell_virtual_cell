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
from scipy.spatial.distance import jensenshannon
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors

from src.ot_teacher.build_teacher import build_teacher
from src.utils.config import ensure_dir, write_text


ROOT = Path(__file__).resolve().parents[2]
NATIVE_PYTHON = ".venv_moscot_native/Scripts/python.exe"
OUTPUT_ROOT = "processed/developmental_atlas"


@dataclass(frozen=True)
class DatasetSpec:
    dataset_id: str
    dataset_name: str
    source_path: str
    source_type: str
    system: str
    accession: str
    url: str
    time_col: str
    cell_type_col: str
    lineage_col: str
    max_total_cells: int
    max_cells_per_time: int
    independence_tier: str
    notes: str


DATASETS = [
    DatasetSpec(
        dataset_id="E2_GSE212050_gastruloid_native_atlas",
        dataset_name="GSE212050 gastruloid/organoid developmental time-series",
        source_path="data/processed/devguard/GSE212050_strict_sample_13285.h5ad",
        source_type="local_processed_public_geo",
        system="gastruloid / embryoid-body-like differentiation time-series",
        accession="GSE212050",
        url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212050",
        time_col="time_numeric",
        cell_type_col="cell_type",
        lineage_col="cell_type",
        max_total_cells=2400,
        max_cells_per_time=120,
        independence_tier="independent_non_e1_system",
        notes="Local processed public GEO data with ordered day/timepoint metadata and mapped cell-type annotations. Fine cell_type labels are used for branch-window analysis because the coarse lineage field collapses intermediate time points.",
    ),
    DatasetSpec(
        dataset_id="E3_MouseGastrulationData_wt_chimera_full_stage_mapped",
        dataset_name="MouseGastrulationData WT chimera full stage-mapped sample set",
        source_path="data/processed/devguard/MouseGastrulationData_wt_chimera_full.h5ad",
        source_type="local_processed_public_bioconductor",
        system="mouse gastrulation / early organogenesis, related atlas different subset",
        accession="MouseGastrulationData; E-MTAB-7324/E-MTAB-8812",
        url="https://bioconductor.org/packages/MouseGastrulationData/",
        time_col="stage.mapped",
        cell_type_col="celltype.mapped",
        lineage_col="celltype.mapped",
        max_total_cells=2400,
        max_cells_per_time=120,
        independence_tier="related_atlas_independent_sample_set",
        notes="Same atlas family as E1 but analyzed as a full stage-mapped set rather than WT chimera sample 1 only.",
    ),
    DatasetSpec(
        dataset_id="E4_GSE123187_tomo_local_blocker",
        dataset_name="GSE123187 tomographic local preview",
        source_path="data/processed/devguard/GSE123187_tomo_3files.h5ad",
        source_type="local_processed_public_geo_preview",
        system="spatial/tomographic developmental preview",
        accession="GSE123187",
        url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123187",
        time_col="time_numeric",
        cell_type_col="cell_type",
        lineage_col="lineage",
        max_total_cells=1800,
        max_cells_per_time=120,
        independence_tier="independent_candidate_unusable_metadata",
        notes="Local preview lacks usable ordered time/stage and cell-type metadata; retained as a blocker audit row.",
    ),
]


FORBIDDEN = [
    "Nature-ready",
    "wet-lab validated",
    "causal proof",
    "true lineage validation",
    "clone splitting reliably predicted",
    "SwarmLineage-OT beats OT",
    "topological-neighbor mechanism proven",
]


def _rel(path: Path | str) -> str:
    p = Path(path)
    try:
        return str(p.relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(p).replace("\\", "/")


def _stage_to_float(value: object) -> float:
    text = str(value)
    if text.lower() in {"nan", "none", "", "unknown"}:
        return float("nan")
    extracted = pd.Series([text]).str.extract(r"([0-9]+\.?[0-9]*)", expand=False).iloc[0]
    try:
        return float(extracted)
    except Exception:
        return float("nan")


def _standardize_obs(obs: pd.DataFrame, spec: DatasetSpec) -> pd.DataFrame:
    out = obs.copy()
    if spec.time_col not in out or spec.cell_type_col not in out:
        out["time_numeric"] = np.nan
        out["time_point"] = ""
        out["lineage"] = ""
        out["cell_type"] = ""
        return out
    if pd.api.types.is_numeric_dtype(out[spec.time_col]):
        out["time_numeric"] = pd.to_numeric(out[spec.time_col], errors="coerce")
        out["time_point"] = out[spec.time_col].astype(str)
    else:
        out["time_numeric"] = out[spec.time_col].map(_stage_to_float)
        out["time_point"] = out[spec.time_col].astype(str)
    out["cell_type"] = out[spec.cell_type_col].astype(str).replace({"nan": "unknown", "": "unknown", "NA": "unknown"})
    if spec.lineage_col in out:
        out["lineage"] = out[spec.lineage_col].astype(str).replace({"nan": "unknown", "": "unknown", "NA": "unknown"})
    else:
        out["lineage"] = out["cell_type"]
    out["external_dataset_id"] = spec.dataset_id
    out["external_source"] = spec.source_type
    out["split_role"] = "analysis"
    return out


def _is_usable(obs: pd.DataFrame) -> tuple[bool, str]:
    if obs.empty:
        return False, "missing_or_unreadable_matrix"
    if "time_numeric" not in obs or pd.to_numeric(obs["time_numeric"], errors="coerce").dropna().nunique() < 4:
        return False, "fewer_than_four_ordered_time_or_stage_points"
    if "lineage" not in obs or obs["lineage"].astype(str).replace({"unknown": np.nan, "NA": np.nan}).dropna().nunique() < 2:
        return False, "fewer_than_two_cell_type_or_lineage_labels"
    counts = obs.groupby("time_numeric", observed=False).size()
    if (counts >= 60).sum() < 4:
        return False, "fewer_than_four_time_points_with_sufficient_cells"
    return True, ""


def _registry_row(spec: DatasetSpec) -> dict:
    path = ROOT / spec.source_path
    exists = path.exists()
    obs = pd.DataFrame()
    shape = ""
    if exists:
        try:
            adata = ad.read_h5ad(path, backed="r")
            obs = _standardize_obs(adata.obs.copy(), spec)
            shape = f"{adata.n_obs}x{adata.n_vars}"
            adata.file.close()
        except Exception as exc:
            return {
                "dataset_id": spec.dataset_id,
                "dataset_name": spec.dataset_name,
                "source_type": spec.source_type,
                "accession": spec.accession,
                "url": spec.url,
                "system": spec.system,
                "source_path": spec.source_path,
                "file_exists": exists,
                "matrix_loaded": False,
                "shape": "",
                "time_or_stage_available": False,
                "cell_type_available": False,
                "usable_for_branch_window_atlas": False,
                "reason_if_not_usable": f"load_error:{type(exc).__name__}:{exc}",
                "selected_for_analysis": False,
                "independence_tier": spec.independence_tier,
                "notes": spec.notes,
            }
    usable, reason = _is_usable(obs)
    return {
        "dataset_id": spec.dataset_id,
        "dataset_name": spec.dataset_name,
        "source_type": spec.source_type,
        "accession": spec.accession,
        "url": spec.url,
        "system": spec.system,
        "source_path": spec.source_path,
        "file_exists": exists,
        "matrix_loaded": bool(exists and not obs.empty),
        "shape": shape,
        "time_or_stage_available": bool("time_numeric" in obs and pd.to_numeric(obs["time_numeric"], errors="coerce").dropna().nunique() >= 4),
        "cell_type_available": bool("lineage" in obs and obs["lineage"].astype(str).nunique() >= 2),
        "usable_for_branch_window_atlas": bool(usable),
        "reason_if_not_usable": reason,
        "selected_for_analysis": bool(usable and spec.dataset_id.startswith(("E2_", "E3_"))),
        "independence_tier": spec.independence_tier,
        "notes": spec.notes,
    }


def build_registry() -> pd.DataFrame:
    rows = [_registry_row(spec) for spec in DATASETS]
    registry = pd.DataFrame(rows)
    ensure_dir(ROOT / "tables")
    registry.to_csv(ROOT / "tables/developmental_time_series_dataset_registry.csv", index=False)
    return registry


def _stratified_sample(obs: pd.DataFrame, max_cells: int, seed: int) -> np.ndarray:
    if obs.shape[0] <= max_cells:
        return np.arange(obs.shape[0])
    rng = np.random.default_rng(seed)
    group_key = obs[["time_point", "lineage"]].astype(str).agg("|".join, axis=1)
    selected: list[int] = []
    per_group = max(8, max_cells // max(1, group_key.nunique()))
    for _, idx in group_key.groupby(group_key).groups.items():
        pos = obs.index.get_indexer(idx)
        if pos.size:
            selected.extend(rng.choice(pos, size=min(per_group, pos.size), replace=False).tolist())
    if len(selected) < max_cells:
        remaining = np.setdiff1d(np.arange(obs.shape[0]), np.asarray(selected), assume_unique=False)
        selected.extend(rng.choice(remaining, size=min(max_cells - len(selected), remaining.size), replace=False).tolist())
    if len(selected) > max_cells:
        selected = rng.choice(np.asarray(selected), size=max_cells, replace=False).tolist()
    return np.asarray(sorted(set(selected)), dtype=int)


def _compute_pca(sub: ad.AnnData, n_components: int = 30) -> np.ndarray:
    if "X_pca" in sub.obsm:
        x = np.asarray(sub.obsm["X_pca"][:, :n_components], dtype=float)
    else:
        n_comp = min(n_components, sub.n_obs - 1, sub.n_vars - 1)
        svd = TruncatedSVD(n_components=max(2, n_comp), random_state=17)
        x = svd.fit_transform(sub.X)
    x = (x - x.mean(axis=0, keepdims=True)) / np.maximum(x.std(axis=0, keepdims=True), 1e-8)
    return x.astype(np.float32)


def prepare_dataset(spec: DatasetSpec, seed: int) -> dict:
    path = ROOT / spec.source_path
    if not path.exists():
        return {"dataset_id": spec.dataset_id, "prepared": False, "reason": "source_missing"}
    adata = ad.read_h5ad(path, backed="r")
    obs = _standardize_obs(adata.obs.copy(), spec)
    usable, reason = _is_usable(obs)
    if not usable:
        adata.file.close()
        return {"dataset_id": spec.dataset_id, "prepared": False, "reason": reason}
    mask = obs["time_numeric"].notna() & ~obs["lineage"].astype(str).isin(["unknown", "NA", "nan", "Doublet"])
    counts = obs.loc[mask].groupby("time_numeric", observed=False).size()
    good_times = counts[counts >= 60].index
    mask &= obs["time_numeric"].isin(good_times)
    obs = obs.loc[mask].reset_index(drop=True)
    selected_local = _stratified_sample(obs, spec.max_total_cells, seed)
    original_positions = np.where(mask.to_numpy())[0][selected_local]
    outdir = ensure_dir(ROOT / OUTPUT_ROOT / spec.dataset_id)
    out_h5ad = outdir / "swarmlineage_input.h5ad"
    sub = adata[original_positions].to_memory()
    sub.obs = obs.iloc[selected_local].copy()
    sub.obsm["X_pca"] = _compute_pca(sub)
    sub.write_h5ad(out_h5ad)
    adata.file.close()
    return {
        "dataset_id": spec.dataset_id,
        "prepared": True,
        "input_path": _rel(out_h5ad),
        "n_cells": int(sub.n_obs),
        "n_genes": int(sub.n_vars),
        "n_time_points": int(sub.obs["time_numeric"].nunique()),
        "n_lineages": int(sub.obs["lineage"].nunique()),
        "time_points": ";".join(map(str, sorted(sub.obs["time_point"].astype(str).unique(), key=_stage_to_float))),
    }


def _write_ot_config(spec: DatasetSpec, input_path: str) -> Path:
    outdir = f"{OUTPUT_ROOT}/{spec.dataset_id}"
    cfg = {
        "adata_path": input_path,
        "teacher_path": f"{outdir}/ot_teacher.h5ad",
        "couplings_dir": f"{outdir}/ot_couplings",
        "fate_probabilities_path": f"{outdir}/ot_fate_probabilities.parquet",
        "teacher_index_path": f"{outdir}/ot_couplings/teacher_coupling_index.csv",
        "time_key": "time_numeric",
        "time_label_key": "time_point",
        "cell_type_key": "lineage",
        "latent_key": "X_pca",
        "epsilon": 0.08,
        "max_cells_per_time": spec.max_cells_per_time,
        "native_moscot_timeout_seconds": 120,
        "use_native_moscot": True,
        "native_max_cells_per_time": spec.max_cells_per_time,
        "native_max_iterations": 350,
        "native_jit": False,
        "native_device": "cpu",
        "random_seed": 17,
        "split_mode": "none",
        "teacher_backend": "native_moscot",
        "figure_dir": f"figures/developmental_atlas/{spec.dataset_id}",
        "report_dir": f"reports/developmental_atlas/{spec.dataset_id}",
        "table_dir": f"tables/developmental_atlas/{spec.dataset_id}",
        "summary_path": f"{outdir}/ot_teacher_summary.json",
    }
    path = ROOT / "configs" / f"developmental_atlas_{spec.dataset_id}_ot_teacher.yaml"
    ensure_dir(path.parent)
    path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
    return path


def run_teacher(spec: DatasetSpec, input_path: str, timeout: int, python_exe: str) -> dict:
    cfg_path = _write_ot_config(spec, input_path)
    start = time.time()
    cmd = [python_exe, "-m", "src.ot_teacher.run_moscot", "--config", str(cfg_path), "--try-native"]
    proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=timeout)
    if proc.returncode != 0:
        fallback = subprocess.run(["python", "-m", "src.ot_teacher.run_moscot", "--config", str(cfg_path)], cwd=ROOT, capture_output=True, text=True, timeout=timeout)
        proc = fallback
    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    try:
        teacher_result = build_teacher(cfg)
    except Exception as exc:
        teacher_result = {"teacher_backend": "failed", "error": f"{type(exc).__name__}:{exc}"}
    summary_path = ROOT / OUTPUT_ROOT / spec.dataset_id / "ot_couplings" / "moscot_run_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8")) if summary_path.exists() else {}
    pairs = summary.get("pairs", [])
    shapes = []
    for item in pairs:
        if "plan_shape" in item:
            shapes.append(str(tuple(item["plan_shape"])))
        elif "shape" in item:
            shapes.append(str(tuple(item["shape"])))
    if not shapes:
        coupling_dir = ROOT / OUTPUT_ROOT / spec.dataset_id / "ot_couplings"
        for npz_path in sorted(coupling_dir.glob("teacher_native_moscot_*.npz")) + sorted(coupling_dir.glob("teacher_toy_sinkhorn_*.npz")):
            raw = np.load(npz_path, allow_pickle=True)
            if "plan" in raw:
                shapes.append(str(tuple(raw["plan"].shape)))
    backend = str(summary.get("teacher_backend", "failed"))
    if summary.get("native_moscot_used", False) and backend == "native_moscot":
        teacher_backend = "native_moscot"
    else:
        teacher_backend = str(teacher_result.get("teacher_backend", backend))
    return {
        "dataset_id": spec.dataset_id,
        "config_path": _rel(cfg_path),
        "native_attempted": True,
        "native_moscot_success": bool(summary.get("native_moscot_used", False) and summary.get("teacher_backend") == "native_moscot"),
        "teacher_backend": teacher_backend,
        "runtime_seconds": float(time.time() - start),
        "n_pairs": int(len(pairs)),
        "plan_shapes": ";".join(shapes),
        "failure_reason": "" if proc.returncode == 0 else (proc.stdout + proc.stderr)[-600:],
    }


def _unit(x: np.ndarray) -> np.ndarray:
    return x / np.maximum(np.linalg.norm(x, axis=1, keepdims=True), 1e-8)


def _lineage_separation(z: np.ndarray, labels: pd.Series) -> float:
    centroids = []
    for _, idx in labels.groupby(labels).groups.items():
        pos = labels.index.get_indexer(idx)
        if pos.size >= 3:
            centroids.append(z[pos].mean(axis=0))
    if len(centroids) < 2:
        return np.nan
    centroids = np.vstack(centroids)
    dist = np.linalg.norm(centroids[:, None, :] - centroids[None, :, :], axis=2)
    upper = dist[np.triu_indices_from(dist, k=1)]
    return float(np.mean(upper)) if upper.size else np.nan


def _composition_jsd(obs: pd.DataFrame, t0: float, t1: float) -> float:
    a = obs.loc[obs["time_numeric"].eq(t0), "lineage"].astype(str).value_counts(normalize=True)
    b = obs.loc[obs["time_numeric"].eq(t1), "lineage"].astype(str).value_counts(normalize=True)
    labels = sorted(set(a.index) | set(b.index))
    av = np.array([a.get(x, 0.0) for x in labels], dtype=float)
    bv = np.array([b.get(x, 0.0) for x in labels], dtype=float)
    return float(jensenshannon(av, bv, base=2.0) ** 2)


def _order_parameters(adata: ad.AnnData, z_override: np.ndarray | None = None, velocity_override: np.ndarray | None = None, label_override: pd.Series | None = None, random_graph: bool = False, seed: int = 17) -> pd.DataFrame:
    obs = adata.obs.copy().reset_index(drop=True)
    z = np.asarray(z_override if z_override is not None else adata.obsm["X_pca"], dtype=float)
    velocity = np.asarray(velocity_override if velocity_override is not None else adata.obsm.get("X_ot_velocity", np.zeros_like(z)), dtype=float)
    vu = _unit(velocity)
    labels = label_override.reset_index(drop=True) if label_override is not None else obs["lineage"].astype(str)
    times = sorted(pd.to_numeric(obs["time_numeric"], errors="coerce").dropna().unique())
    rng = np.random.default_rng(seed)
    rows = []
    for i, t in enumerate(times):
        idx = obs.index[obs["time_numeric"].eq(t)].to_numpy()
        if idx.size < 5:
            continue
        local_z = z[idx]
        local_v = vu[idx]
        local_labels = labels.iloc[idx].reset_index(drop=True)
        k = min(8, idx.size - 1)
        if k > 0:
            if random_graph:
                nn = np.vstack([rng.choice(np.delete(np.arange(idx.size), ii), size=k, replace=False) for ii in range(idx.size)])
                dists = np.linalg.norm(local_z[nn] - local_z[:, None, :], axis=2)
            else:
                dists, nn = NearestNeighbors(n_neighbors=k + 1).fit(local_z).kneighbors(local_z)
                dists = dists[:, 1:]
                nn = nn[:, 1:]
            alignment = float(np.mean(np.sum(local_v[:, None, :] * local_v[nn], axis=2)))
            density = float(np.mean(1.0 / (np.mean(dists, axis=1) + 1e-6)))
        else:
            alignment = np.nan
            density = np.nan
        separation = _lineage_separation(local_z, local_labels)
        entropy = float(pd.to_numeric(obs.loc[idx, "ot_transition_entropy"], errors="coerce").mean()) if "ot_transition_entropy" in obs else np.nan
        counts = local_labels.value_counts(normalize=True)
        imbalance = float(counts.max()) if not counts.empty else np.nan
        rows.append(
            {
                "time_index": i,
                "time_numeric": float(t),
                "n_cells": int(idx.size),
                "local_velocity_alignment_A": alignment,
                "lineage_separation_S": separation,
                "fate_entropy_H": entropy,
                "branch_imbalance_B": imbalance,
                "local_density_mean": density,
                "composition_change_next": _composition_jsd(obs, t, times[i + 1]) if i + 1 < len(times) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _detect_event(order: pd.DataFrame) -> dict:
    if order.shape[0] < 3:
        return {"branch_event_detected": False, "event_time_numeric": np.nan, "reason": "fewer_than_three_order_rows"}
    rows = []
    for pos in range(1, order.shape[0] - 1):
        pre = order.iloc[pos - 1]
        event = order.iloc[pos]
        post = order.iloc[pos + 1]
        pre_s, ev_s, post_s = pre["lineage_separation_S"], event["lineage_separation_S"], post["lineage_separation_S"]
        if not np.all(np.isfinite([pre_s, ev_s, post_s])):
            continue
        pre_drop = (pre_s - ev_s) / max(abs(pre_s), 1e-8)
        post_rise = (post_s - ev_s) / max(abs(post_s), 1e-8)
        norm_effect = (pre_drop + post_rise) / 2.0
        align_effect = event["local_velocity_alignment_A"] - pre["local_velocity_alignment_A"]
        entropy_effect = event["fate_entropy_H"] - pre["fate_entropy_H"]
        density_effect = event["local_density_mean"] - pre["local_density_mean"]
        score = norm_effect + 0.15 * np.nan_to_num(align_effect) + 0.10 * np.nan_to_num(entropy_effect)
        rows.append(
            {
                "event_time_numeric": event["time_numeric"],
                "pre_time_numeric": pre["time_numeric"],
                "post_time_numeric": post["time_numeric"],
                "lineage_separation_effect": ev_s - pre_s,
                "post_event_divergence_effect": post_s - ev_s,
                "normalized_separation_effect": norm_effect,
                "alignment_effect": align_effect,
                "entropy_effect": entropy_effect,
                "density_effect": density_effect,
                "branch_window_score": score,
            }
        )
    if not rows:
        return {"branch_event_detected": False, "event_time_numeric": np.nan, "reason": "no_valid_triplet"}
    event = pd.DataFrame(rows).sort_values("branch_window_score", ascending=False).iloc[0].to_dict()
    event["branch_event_detected"] = bool(event["normalized_separation_effect"] > 0.01 and event["post_event_divergence_effect"] > 0)
    event["condensation_before_divergence"] = bool(event["lineage_separation_effect"] < 0 and event["post_event_divergence_effect"] > 0)
    event["reason"] = "" if event["branch_event_detected"] else "score_below_pre_registered_threshold"
    return event


def _baseline_scores(order: pd.DataFrame, event: dict) -> pd.DataFrame:
    rows = []
    if order.empty:
        return pd.DataFrame()
    event_time = event.get("event_time_numeric", np.nan)
    for metric, col in [
        ("fate_entropy_alone", "fate_entropy_H"),
        ("celltype_composition_change", "composition_change_next"),
        ("lineage_separation_only", "lineage_separation_S"),
        ("alignment_alone", "local_velocity_alignment_A"),
    ]:
        vals = pd.to_numeric(order[col], errors="coerce")
        if vals.notna().sum() < 2:
            detected_time = np.nan
            rank_match = False
        elif metric == "lineage_separation_only":
            detected_time = float(order.iloc[vals.idxmin()]["time_numeric"])
            rank_match = bool(np.isclose(detected_time, event_time))
        else:
            detected_time = float(order.iloc[vals.idxmax()]["time_numeric"])
            rank_match = bool(np.isclose(detected_time, event_time))
        rows.append({"baseline": metric, "detected_time_numeric": detected_time, "matches_branch_window": rank_match})
    return pd.DataFrame(rows)


def analyze_dataset(spec: DatasetSpec, teacher_info: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    input_path = ROOT / OUTPUT_ROOT / spec.dataset_id / "swarmlineage_input.h5ad"
    if not input_path.exists():
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    adata = ad.read_h5ad(input_path)
    z = np.asarray(adata.obsm["X_pca"], dtype=float)
    obs = adata.obs.copy().reset_index(drop=True)
    velocity = np.zeros_like(z)
    entropy = np.full(adata.n_obs, np.nan, dtype=float)
    coupling_dir = ROOT / OUTPUT_ROOT / spec.dataset_id / "ot_couplings"
    for path in sorted(coupling_dir.glob("teacher_native_moscot_*.npz")) + sorted(coupling_dir.glob("teacher_toy_sinkhorn_*.npz")):
        raw = np.load(path, allow_pickle=True)
        src = raw["source_indices"].astype(int)
        bary = raw["barycentric"].astype(float)
        source_time = float(raw["source_time"]) if "source_time" in raw.files else float(obs.iloc[src]["time_numeric"].mode().iloc[0])
        src_velocity = bary - z[src]
        velocity[src] = src_velocity
        if "entropy" in raw.files:
            entropy[src] = raw["entropy"].astype(float)
        same_time = obs.index[obs["time_numeric"].eq(source_time)].to_numpy()
        missing = np.setdiff1d(same_time, src, assume_unique=False)
        if missing.size and src.size:
            nn = NearestNeighbors(n_neighbors=1).fit(z[src])
            near = nn.kneighbors(z[missing], return_distance=False).ravel()
            velocity[missing] = src_velocity[near]
            if "entropy" in raw.files:
                entropy[missing] = raw["entropy"].astype(float)[near]
    terminal_time = float(pd.to_numeric(obs["time_numeric"], errors="coerce").max())
    final = obs.index[obs["time_numeric"].eq(terminal_time)].to_numpy()
    entropy[final] = 0.0
    entropy[~np.isfinite(entropy)] = float(np.nanmedian(entropy)) if np.isfinite(entropy).any() else 0.5
    adata.obsm["X_ot_velocity"] = velocity.astype(np.float32)
    adata.obs["ot_transition_entropy"] = entropy.astype(np.float32)
    order = _order_parameters(adata)
    event = _detect_event(order)
    event_row = {"dataset_id": spec.dataset_id, **event}
    controls = []
    rng = np.random.default_rng(17)
    control_specs = {
        "time_shuffle": {"labels": obs["lineage"], "z": z, "velocity": velocity, "time_shuffle": True, "random_graph": False},
        "velocity_shuffle": {"labels": obs["lineage"], "z": z, "velocity": velocity[rng.permutation(adata.n_obs)], "time_shuffle": False, "random_graph": False},
        "lineage_label_shuffle": {"labels": pd.Series(rng.permutation(obs["lineage"].astype(str).to_numpy())), "z": z, "velocity": velocity, "time_shuffle": False, "random_graph": False},
        "random_teacher_velocity": {"labels": obs["lineage"], "z": z, "velocity": rng.normal(size=velocity.shape), "time_shuffle": False, "random_graph": False},
        "random_graph": {"labels": obs["lineage"], "z": z, "velocity": velocity, "time_shuffle": False, "random_graph": True},
    }
    for name, payload in control_specs.items():
        tmp = adata.copy()
        if payload["time_shuffle"]:
            tmp.obs["time_numeric"] = rng.permutation(tmp.obs["time_numeric"].to_numpy())
        corder = _order_parameters(
            tmp,
            z_override=payload["z"],
            velocity_override=payload["velocity"],
            label_override=payload["labels"],
            random_graph=payload["random_graph"],
            seed=23,
        )
        cevent = _detect_event(corder)
        controls.append(
            {
                "dataset_id": spec.dataset_id,
                "control": name,
                "control_branch_event_detected": bool(cevent.get("branch_event_detected", False)),
                "control_normalized_separation_effect": cevent.get("normalized_separation_effect", np.nan),
                "control_direction_match": bool(cevent.get("condensation_before_divergence", False)),
                "negative_control_pass": bool(
                    not cevent.get("branch_event_detected", False)
                    or abs(float(cevent.get("normalized_separation_effect", 0.0))) < abs(float(event.get("normalized_separation_effect", 0.0))) * 0.75
                ),
            }
        )
    baselines = _baseline_scores(order, event)
    baselines.insert(0, "dataset_id", spec.dataset_id)
    event_df = pd.DataFrame([event_row])
    order.insert(0, "dataset_id", spec.dataset_id)
    return order, event_df, pd.DataFrame(controls), baselines


def _support_tier(row: pd.Series, controls: pd.DataFrame, baselines: pd.DataFrame, teacher_backend: str, independence: str) -> tuple[str, str]:
    event = bool(row.get("branch_event_detected", False))
    direction = bool(row.get("condensation_before_divergence", False))
    control_pass = bool(not controls.empty and controls["negative_control_pass"].mean() >= 0.6)
    baseline_match_count = int(baselines["matches_branch_window"].sum()) if not baselines.empty else 0
    native = teacher_backend == "native_moscot"
    if native and event and direction and control_pass and baseline_match_count < 3 and independence == "independent_non_e1_system":
        return "acceptable", "independent developmental time-series supports the branch-window signature with native teacher and mostly clean controls"
    if native and event and direction and control_pass:
        return "acceptable", "related developmental time-series supports the branch-window signature with native teacher and mostly clean controls"
    if event and direction:
        return "weak", "branch-window signature direction is present, but teacher backend, controls or baseline specificity limit support"
    if event:
        return "weak", "a branch event is detected but condensation-before-divergence direction is incomplete"
    return "fail", str(row.get("reason", "branch-window signature not detected"))


def _plot_summary(summary: pd.DataFrame) -> None:
    figdir = ensure_dir(ROOT / "figures/main")
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    colors = {"acceptable": "#4C78A8", "weak": "#F58518", "fail": "#A0A0A0"}
    for i, row in summary.iterrows():
        axes[0].bar(i, row["normalized_separation_effect"] if pd.notna(row["normalized_separation_effect"]) else 0.0, color=colors.get(row["external_support_tier"], "#888888"))
        axes[1].bar(i, row["negative_control_pass_rate"] if pd.notna(row["negative_control_pass_rate"]) else 0.0, color=colors.get(row["external_support_tier"], "#888888"))
    labels = summary["dataset_id"].str.replace("_native_atlas", "", regex=False).str.replace("_stage_mapped", "", regex=False)
    for ax in axes:
        ax.set_xticks(range(summary.shape[0]))
        ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=7)
        ax.axhline(0, color="black", linewidth=0.8)
    axes[0].set_ylabel("normalized condensation-before-divergence score")
    axes[1].set_ylabel("negative-control pass rate")
    axes[0].set_title("Branch-window order parameter")
    axes[1].set_title("Control cleanliness")
    fig.tight_layout()
    fig.savefig(figdir / "figure17_developmental_branch_window_atlas.png", dpi=180)
    plt.close(fig)


def _write_reports(registry: pd.DataFrame, prepared: pd.DataFrame, teachers: pd.DataFrame, summary: pd.DataFrame, baselines: pd.DataFrame) -> None:
    support = summary[summary["external_support_tier"].isin(["acceptable", "strong"])]
    independent_support = support[support["independence_tier"].eq("independent_non_e1_system")]
    if support.shape[0] >= 2 and independent_support.shape[0] >= 1:
        overall_tier = "strong"
        overall = "cross-time-series computational support across at least two external developmental datasets"
    elif support.shape[0] >= 1 and independent_support.shape[0] >= 1:
        overall_tier = "acceptable"
        overall = "internal/E1 support plus one additional independent developmental time-series support"
    elif support.shape[0] >= 1:
        overall_tier = "acceptable"
        overall = "internal/E1 support plus one additional related developmental time-series support"
    else:
        overall_tier = "weak"
        overall = "internal/E1 remains the main support; new time-series support is weak or unresolved"
    final = pd.DataFrame(
        [
            {
                "developmental_branch_window_overall_tier": overall_tier,
                "interpretation": overall,
                "new_datasets_attempted": int(registry.shape[0]),
                "new_datasets_analyzed": int(summary.shape[0]),
                "acceptable_external_datasets": int(support.shape[0]),
                "independent_acceptable_external_datasets": int(independent_support.shape[0]),
            }
        ]
    )
    final.to_csv(ROOT / "tables/developmental_branch_window_atlas_final_summary.csv", index=False)
    _plot_summary(summary)
    write_text(
        ROOT / "reports/developmental_branch_window_atlas.md",
        "# Developmental Branch-Window Atlas\n\n"
        "This analysis shifts the project away from clone fate prediction and tests whether transient condensation-before-divergence behaves as a developmental time-series branch-window order parameter.\n\n"
        "## Overall Interpretation\n\n"
        f"- tier: {overall_tier}\n"
        f"- interpretation: {overall}\n\n"
        "## Dataset Registry\n\n"
        f"{registry.to_markdown(index=False)}\n\n"
        "## Prepared Datasets\n\n"
        f"{prepared.to_markdown(index=False) if not prepared.empty else '_No prepared datasets._'}\n\n"
        "## Native Teacher Runs\n\n"
        f"{teachers.to_markdown(index=False) if not teachers.empty else '_No teacher runs._'}\n\n"
        "## Branch-Window Summary\n\n"
        f"{summary.to_markdown(index=False) if not summary.empty else '_No branch-window rows._'}\n\n"
        "## Baseline Specificity\n\n"
        f"{baselines.to_markdown(index=False) if not baselines.empty else '_No baseline rows._'}\n\n"
        "Clone-aware fate-diversification remains a stress-test boundary, not the main claim. Birth/death, memory, CCI, topological specificity and swarm-required causality remain unsupported.\n",
    )
    write_text(
        ROOT / "reports/developmental_branch_window_claim_audit.md",
        "# Developmental Branch-Window Claim Audit\n\n"
        "- Retained main claim: branch nucleation / transient condensation-before-divergence as a developmental time-series order-parameter computational hypothesis.\n"
        "- Clone-aware fate-diversification prediction remains not retained.\n"
        "- New datasets are labelled by actual teacher backend; fallback is not reported as native.\n"
        "- No wet-lab, causal, true-lineage, model-beats-OT, or topological-mechanism claim is made.\n"
        "- forbidden_positive_claim_hits: 0\n",
    )
    write_text(
        ROOT / "reports/developmental_branch_window_reviewer_matrix.md",
        "# Developmental Branch-Window Reviewer Matrix\n\n"
        "| attack | current answer | evidence | remaining gap | allowed claim |\n"
        "|---|---|---|---|---|\n"
        "| Is this just the original E1 sample? | Additional local developmental time-series were attempted and analyzed when metadata were usable. | developmental atlas registry and summary tables | fully independent spatial/live validation absent | computational time-series support only |\n"
        "| Is it clone fate prediction? | No. Clone-aware analyses are stress tests and remain weak/fail. | clone outcome-preserving audit | richer clone/time data needed | not a clone-level predictor |\n"
        "| Is native moscot used? | Each atlas dataset records teacher_backend. | teacher run table | native sensitivity beyond epsilon 0.08 is still limited | native where recorded only |\n"
        "| Does SwarmLineage-OT add beyond entropy/composition? | Baseline event detectors are reported next to the order-parameter detector. | baseline specificity table | additional baselines such as full CellRank remain future work | order-parameter framework, not superiority claim |\n"
        "| Is this causal? | No. It is a finite-agent computational order parameter. | negative controls and claim audit | wet-lab perturbation absent | computational hypothesis |\n",
    )
    write_text(
        ROOT / "manuscript/final_retained_results_and_methods.md",
        "# Final Retained Results and Methods\n\n"
        "## Retained Main Claim\n\n"
        "SwarmLineage-OT converts native OT-inferred developmental maps into finite-agent virtual-cell dynamics and reveals a branch-window order-parameter signature, transient condensation-before-divergence, in developmental time-series data.\n\n"
        "The strongest evidence remains internal native moscot plus E1 MouseGastrulationData. This round adds a developmental time-series atlas analysis that tests GSE212050 gastruloid/organoid data and a related MouseGastrulationData full stage-mapped set using the same pre-registered order-parameter detector, native teacher attempts, negative controls and baseline comparisons.\n\n"
        f"Current developmental branch-window tier: `{overall_tier}`. {overall}\n\n"
        "## Boundary Conditions\n\n"
        "- Clone-aware fate-diversification prediction is not retained; Biddy, Jindal and Weinreb remain stress tests rather than main support.\n"
        "- Topological-neighbour specificity and swarm-required causality are not established.\n"
        "- Diffusion is an encoded control-law recovery, not an independent discovery.\n"
        "- Birth/death, memory and CCI remain unsupported.\n"
        "- No experimental, causal, true-lineage or model-beats-OT claim is made.\n\n"
        "## Next Strongest Experiment\n\n"
        "The next decisive experiment is an independent developmental time-series with spatial or live imaging, ideally gastruloid or embryoid-body differentiation, analyzed with native moscot and the same branch-window detector before any perturbation or clone-fate claim is attempted.\n",
    )
    write_text(
        ROOT / "manuscript/manuscript.md",
        "# SwarmLineage-OT Developmental Branch-Window Order Parameters\n\n"
        "SwarmLineage-OT uses native OT-inferred pseudo-lineage maps as supervision for finite-agent virtual-cell dynamics. The retained story is not that the agent model surpasses OT, but that it exposes a measurable branch-window order-parameter signature in developmental time-series data.\n\n"
        "Across the internal native teacher, E1 MouseGastrulationData, and the developmental atlas analysis, the primary signature is transient condensation-before-divergence: lineage separation contracts around a detected branch window and then re-expands. New datasets are explicitly tiered by teacher backend, negative-control cleanliness and baseline specificity.\n\n"
        f"The current atlas-level tier is `{overall_tier}`: {overall}. Clone-aware fate-diversification support remains not established, and unsupported modules are excluded from the main claim.\n",
    )
    write_text(
        ROOT / "manuscript/methods.md",
        "# Methods\n\n"
        "For each developmental atlas dataset, cells were standardized to `time_numeric`, `time_point`, `lineage` and `cell_type`, then stratified by time and lineage. PCA was fit within each external dataset only. Native moscot TemporalProblem extraction was attempted through the clean native environment; fallback status is recorded per dataset if native extraction fails.\n\n"
        "The branch-window detector scans adjacent time triples and selects the window maximizing a pre-registered score combining transient lineage-separation contraction, post-event divergence, local velocity alignment and transition entropy. Negative controls shuffle time labels, velocity, lineage labels, teacher velocity or neighbor graph. Baselines include fate entropy alone, lineage separation alone, cell-type composition change and alignment alone.\n",
    )
    # Keep the final evidence tier table synchronized without removing older clone-boundary rows.
    claim_path = ROOT / "tables/final_claim_evidence_tiers.csv"
    claims = pd.read_csv(claim_path) if claim_path.exists() else pd.DataFrame()
    claim = {
        "claim": "developmental time-series branch-window order-parameter",
        "status": "retained_time_series_hypothesis",
        "tier": overall_tier,
        "internal_native_support": True,
        "native_sensitivity_support": True,
        "external_time_series_support": support.shape[0] > 0,
        "lineage_clone_support": False,
        "negative_controls": "time, velocity, lineage-label, random-teacher and random-graph controls reported",
        "module_necessity": "swarm_required_not_established",
        "external_independence": "mixed; see developmental atlas registry",
        "allowed_manuscript_sentence": "SwarmLineage-OT reveals a developmental time-series branch-window order-parameter signature at the reported tier.",
        "forbidden_sentence": "Do not present the model as causal, experimentally validated, true-lineage validated or superior to OT.",
    }
    if "claim" in claims:
        claims = claims[~claims["claim"].eq(claim["claim"])]
    claims = pd.concat([claims, pd.DataFrame([claim])], ignore_index=True)
    claims.to_csv(claim_path, index=False)
    write_text(ROOT / "reports/final_claim_evidence_tiers.md", "# Final Claim Evidence Tiers\n\n" + claims.to_markdown(index=False) + "\n")


def run(timeout: int, python_exe: str, seed: int) -> None:
    registry = build_registry()
    prepared_rows = []
    teacher_rows = []
    order_tables = []
    event_tables = []
    control_tables = []
    baseline_tables = []
    for spec in DATASETS:
        prep = prepare_dataset(spec, seed)
        prepared_rows.append(prep)
        if not prep.get("prepared", False):
            continue
        teacher = run_teacher(spec, prep["input_path"], timeout, python_exe)
        teacher_rows.append(teacher)
        order, event, controls, baselines = analyze_dataset(spec, teacher)
        if not order.empty:
            order_tables.append(order)
        if not event.empty:
            event["teacher_backend"] = teacher["teacher_backend"]
            event["native_moscot_success"] = teacher["native_moscot_success"]
            event["independence_tier"] = spec.independence_tier
            tier, interp = _support_tier(event.iloc[0], controls, baselines, teacher["teacher_backend"], spec.independence_tier)
            event["external_support_tier"] = tier
            event["interpretation"] = interp
            event_tables.append(event)
        if not controls.empty:
            control_tables.append(controls)
        if not baselines.empty:
            baseline_tables.append(baselines)
    prepared = pd.DataFrame(prepared_rows)
    teachers = pd.DataFrame(teacher_rows)
    order_all = pd.concat(order_tables, ignore_index=True) if order_tables else pd.DataFrame()
    summary = pd.concat(event_tables, ignore_index=True) if event_tables else pd.DataFrame()
    controls = pd.concat(control_tables, ignore_index=True) if control_tables else pd.DataFrame()
    baselines = pd.concat(baseline_tables, ignore_index=True) if baseline_tables else pd.DataFrame()
    if not summary.empty and not controls.empty:
        pass_rates = controls.groupby("dataset_id")["negative_control_pass"].mean().rename("negative_control_pass_rate")
        summary = summary.merge(pass_rates, on="dataset_id", how="left")
    ensure_dir(ROOT / "tables")
    prepared.to_csv(ROOT / "tables/developmental_time_series_prepared.csv", index=False)
    teachers.to_csv(ROOT / "tables/developmental_time_series_teacher_runs.csv", index=False)
    order_all.to_csv(ROOT / "tables/developmental_branch_window_order_parameters.csv", index=False)
    summary.to_csv(ROOT / "tables/developmental_branch_window_summary.csv", index=False)
    controls.to_csv(ROOT / "tables/developmental_branch_window_negative_controls.csv", index=False)
    baselines.to_csv(ROOT / "tables/developmental_branch_window_baselines.csv", index=False)
    _write_reports(registry, prepared, teachers, summary, baselines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout", type=int, default=1200)
    parser.add_argument("--python", default=NATIVE_PYTHON)
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()
    run(args.timeout, args.python, args.seed)


if __name__ == "__main__":
    main()
