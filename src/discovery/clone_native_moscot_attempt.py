from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from sklearn.preprocessing import OneHotEncoder

from src.discovery.clone_developmental_validation_v14 import DATASETS, _clean_clone, _parse_time
from src.utils.config import ensure_dir


ROOT = Path(__file__).resolve().parents[2]


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
        enc = OneHotEncoder(sparse_output=False, handle_unknown="ignore", max_categories=30)
        one_hot = enc.fit_transform(obs[[cfg.celltype_col]].astype(str))
        x = np.hstack([score, one_hot])
    x = (x - x.mean(axis=0, keepdims=True)) / np.maximum(x.std(axis=0, keepdims=True), 1e-8)
    return x[:, : min(30, x.shape[1])].astype(np.float32)


def prepare_dataset(cfg, max_cells_per_time: int = 120, seed: int = 17) -> dict:
    raw = ROOT / "data/external_l3/raw" / cfg.file_name
    outdir = ensure_dir(ROOT / "processed/external_l4_native" / cfg.dataset_id)
    out_h5ad = outdir / "native_input.h5ad"
    if not raw.exists():
        return {"dataset_id": cfg.dataset_id, "prepared": False, "reason": "raw_file_missing"}
    adata = ad.read_h5ad(raw, backed="r")
    obs = adata.obs.copy()
    if cfg.time_col is None or cfg.time_col not in obs.columns:
        adata.file.close()
        return {"dataset_id": cfg.dataset_id, "prepared": False, "reason": "missing_time_col"}
    obs["time_numeric"] = _parse_time(obs[cfg.time_col], cfg.dataset_id)
    obs["time_point"] = obs[cfg.time_col].astype(str)
    obs["clone_id"] = _clean_clone(obs[cfg.clone_col])
    obs["lineage"] = obs[cfg.celltype_col].astype(str)
    valid = obs["time_numeric"].notna() & obs["clone_id"].notna()
    times = sorted(obs.loc[valid, "time_numeric"].dropna().unique())
    if len(times) < 2:
        adata.file.close()
        return {"dataset_id": cfg.dataset_id, "prepared": False, "reason": "fewer_than_two_times"}
    rng = np.random.default_rng(seed)
    selected = []
    for i, t in enumerate(times):
        idx = np.where(valid.to_numpy() & obs["time_numeric"].eq(t).to_numpy())[0]
        labels = obs.iloc[idx]["lineage"].astype(str)
        per_type = max(10, max_cells_per_time // max(labels.nunique(), 1))
        keep = []
        for _, members in labels.groupby(labels).groups.items():
            local = obs.index.get_indexer(members)
            local = local[local >= 0]
            keep.extend(rng.choice(local, size=min(per_type, local.size), replace=False).tolist())
        if len(keep) < max_cells_per_time:
            rem = np.setdiff1d(idx, np.asarray(keep), assume_unique=False)
            keep.extend(rng.choice(rem, size=min(max_cells_per_time - len(keep), rem.size), replace=False).tolist())
        selected.extend(sorted(keep[:max_cells_per_time]))
    selected = np.asarray(sorted(set(selected)), dtype=int)
    sub = adata[selected].to_memory()
    sub.obs = obs.iloc[selected].copy()
    sub.obsm["X_pca"] = _latent(sub.obs, cfg)
    sub.write_h5ad(out_h5ad)
    adata.file.close()
    return {"dataset_id": cfg.dataset_id, "prepared": True, "path": _rel(out_h5ad), "n_obs": int(sub.n_obs), "times": [float(t) for t in times]}


def write_config(dataset_id: str, h5ad_path: str, max_cells_per_time: int) -> Path:
    cfg_path = ROOT / "configs" / f"clone_native_{dataset_id}.yaml"
    outdir = f"processed/external_l4_native/{dataset_id}"
    cfg = {
        "adata_path": h5ad_path.replace("\\", "/"),
        "teacher_path": f"{outdir}/ot_teacher.h5ad",
        "couplings_dir": f"{outdir}/ot_couplings",
        "teacher_index_path": f"{outdir}/ot_couplings/teacher_coupling_index.csv",
        "summary_path": f"{outdir}/ot_teacher_summary.json",
        "time_key": "time_numeric",
        "time_label_key": "time_point",
        "cell_type_key": "lineage",
        "latent_key": "X_pca",
        "epsilon": 0.08,
        "max_cells_per_time": max_cells_per_time,
        "native_moscot_timeout_seconds": 120,
        "use_native_moscot": True,
        "native_max_cells_per_time": max_cells_per_time,
        "native_max_iterations": 350,
        "native_jit": False,
        "native_device": "cpu",
        "random_seed": 17,
        "split_mode": "none",
        "teacher_backend": "native_moscot",
    }
    with open(cfg_path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(cfg, handle, sort_keys=False)
    return cfg_path


def run_native(dataset_id: str, cfg_path: Path, python_exe: str, timeout: int) -> dict:
    cmd = [python_exe, "-m", "src.ot_teacher.run_moscot", "--config", str(cfg_path), "--try-native"]
    start = pd.Timestamp.utcnow()
    try:
        out = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, timeout=timeout)
        detail = (out.stdout + "\n" + out.stderr).strip()
        summary_path = ROOT / "processed/external_l4_native" / dataset_id / "ot_couplings" / "moscot_run_summary.json"
        backend = "unknown"
        used = False
        pairs = 0
        if summary_path.exists():
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            backend = summary.get("teacher_backend", "unknown")
            used = bool(summary.get("native_moscot_used", False))
            pairs = len(summary.get("pairs", []))
        return {
            "dataset_id": dataset_id,
            "attempted": True,
            "returncode": out.returncode,
            "native_moscot_used": used,
            "teacher_backend": backend,
            "pairs": pairs,
            "detail": detail[-1200:],
            "started_utc": str(start),
        }
    except subprocess.TimeoutExpired:
        return {
            "dataset_id": dataset_id,
            "attempted": True,
            "returncode": "timeout",
            "native_moscot_used": False,
            "teacher_backend": "timeout",
            "pairs": 0,
            "detail": f"native run timed out after {timeout}s",
            "started_utc": str(start),
        }


def run() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", default=".venv_moscot_native/Scripts/python.exe")
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--max-cells-per-time", type=int, default=120)
    args = parser.parse_args()
    rows = []
    for cfg in [d for d in DATASETS if d.dataset_id.startswith("Jindal") or d.dataset_id.startswith("Weinreb")]:
        prep = prepare_dataset(cfg, max_cells_per_time=args.max_cells_per_time)
        row = dict(prep)
        if prep.get("prepared", False):
            cfg_path = write_config(cfg.dataset_id, prep["path"], args.max_cells_per_time)
            row.update(run_native(cfg.dataset_id, cfg_path, args.python, args.timeout))
            row["config_path"] = _rel(cfg_path)
        rows.append(row)
    out = pd.DataFrame(rows)
    ensure_dir(ROOT / "tables")
    out.to_csv(ROOT / "tables/clone_native_moscot_attempts.csv", index=False)
    ensure_dir(ROOT / "reports")
    (ROOT / "reports/clone_native_moscot_attempts.md").write_text(
        "# Clone Native moscot Attempts\n\n"
        + out.to_markdown(index=False)
        + "\n",
        encoding="utf-8",
    )
    print(out.to_json(orient="records"))


if __name__ == "__main__":
    run()
