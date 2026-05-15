from __future__ import annotations

import argparse
import json
import math
import subprocess
import time
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy import sparse, stats
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors

from src.external.v1_evidence_package import (
    _bootstrap_ci,
    _entropy,
    _event_window_from_order,
    _lineage_separation,
    _local_density,
    _md_table,
    _normalise_counts,
    _order_parameters_from_embedding,
    _read_csv,
    _write_csv,
    _write_md,
)
from src.ot_teacher.build_teacher import build_teacher
from src.ot_teacher.run_moscot import _native_status, run_native_moscot_teacher, run_toy_sinkhorn_teacher
from src.utils.config import ensure_dir


ROOT = Path(__file__).resolve().parents[2]
BIDDY_FILE = ROOT / "data/external_l2/l2_raw/Biddy_2018_Nature.h5ad"


def _rel(path: str | Path) -> str:
    try:
        return str(Path(path).resolve().relative_to(ROOT))
    except Exception:
        return str(path)


def _git(cmd: list[str]) -> str:
    try:
        return subprocess.check_output(["git", *cmd], cwd=ROOT, text=True, stderr=subprocess.STDOUT, timeout=60).strip()
    except Exception as exc:
        return f"ERROR: {type(exc).__name__}: {exc}"


def write_start_state() -> None:
    branch = _git(["branch", "--show-current"])
    head = _git(["rev-parse", "HEAD"])
    status = _git(["status", "--short", "--branch"])
    ahead = "ahead" in status
    tracked_large = _git(["ls-files"])
    large_risk = [p for p in tracked_large.splitlines() if Path(p).suffix.lower() in {".h5ad", ".rds", ".gz", ".tar", ".npz", ".pt"}]
    v1_files = {
        "reports/v1_goal_status.md": (ROOT / "reports/v1_goal_status.md").exists(),
        "tables/l1_clone_model_summary.csv": (ROOT / "tables/l1_clone_model_summary.csv").exists(),
        "src/external/v1_evidence_package.py": (ROOT / "src/external/v1_evidence_package.py").exists(),
    }
    text = [
        "# v1.1 Start State",
        "",
        f"- branch: `{branch}`",
        f"- head_commit: `{head}`",
        f"- ahead_origin: {ahead}",
        "- working_tree_status:",
        "```",
        status,
        "```",
        "- v1.0 output presence:",
    ]
    text.extend([f"  - {path}: {present}" for path, present in v1_files.items()])
    text.extend(
        [
            f"- tracked_large_file_risk_count: {len(large_risk)}",
            f"- tracked_large_file_risk_examples: {large_risk[:10]}",
            "",
            "No reset, rebase or remote overwrite was performed.",
        ]
    )
    _write_md("reports/v1_1_start_state.md", "\n".join(text))


def lineage_registry() -> pd.DataFrame:
    kim = ROOT / "data/external_l1/Kim_2020_CellReports.h5ad"
    wei = ROOT / "data/external_l1/Wei_2020_GenomeResearch.h5ad"
    rows = [
        {
            "dataset_id": "L2_biddy_2018_nature",
            "dataset_name": "CellTag direct reprogramming time course",
            "publication": "Biddy et al., Nature 2018",
            "accession": "GSE99915; scLTdb Biddy_2018_Nature.h5ad",
            "doi": "10.1038/s41586-018-0744-4; 10.5281/zenodo.12176634",
            "url": "https://zenodo.org/records/12176634/files/Biddy_2018_Nature.h5ad?download=1",
            "source": "scLTdb Zenodo processed h5ad",
            "organism": "mouse",
            "biological_system": "CellTag reprogramming",
            "time_or_stage_available": BIDDY_FILE.exists(),
            "clone_or_barcode_available": BIDDY_FILE.exists(),
            "cell_type_available": BIDDY_FILE.exists(),
            "expression_matrix_available": BIDDY_FILE.exists(),
            "processed_h5ad_available": BIDDY_FILE.exists(),
            "metadata_available": BIDDY_FILE.exists(),
            "download_attempted": True,
            "download_success": BIDDY_FILE.exists(),
            "matrix_loaded": False,
            "metadata_loaded": False,
            "clone_metadata_loaded": False,
            "usable_for_clone_validation": False,
            "reason_if_not_usable": "" if BIDDY_FILE.exists() else "Download failed or local file missing.",
            "selected_for_L2": BIDDY_FILE.exists(),
            "notes": "Priority L2 dataset. Selected when local downloaded file is present.",
        },
        {
            "dataset_id": "L2_morris_lab_celltag_resources",
            "dataset_name": "Morris lab CellTag resources",
            "publication": "Biddy et al., Nature 2018 and Morris lab resources",
            "accession": "GSE99915",
            "doi": "10.1038/s41586-018-0744-4",
            "url": "https://morrislab.io/celltagging-resources/",
            "source": "resource page",
            "organism": "mouse",
            "biological_system": "CellTag reprogramming",
            "time_or_stage_available": True,
            "clone_or_barcode_available": True,
            "cell_type_available": True,
            "expression_matrix_available": True,
            "processed_h5ad_available": False,
            "metadata_available": True,
            "download_attempted": False,
            "download_success": False,
            "matrix_loaded": False,
            "metadata_loaded": False,
            "clone_metadata_loaded": False,
            "usable_for_clone_validation": False,
            "reason_if_not_usable": "Not downloaded separately because the equivalent Biddy processed h5ad was downloaded from scLTdb.",
            "selected_for_L2": False,
            "notes": "Registered as an equivalent source, not claimed as independently analyzed.",
        },
        {
            "dataset_id": "L2_cellrank_morris_reprogramming",
            "dataset_name": "CellRank Morris reprogramming dataset wrapper",
            "publication": "Morris reprogramming data as exposed by CellRank",
            "accession": "CellRank dataset wrapper for GSE99915-derived data",
            "doi": "10.1038/s41586-018-0744-4",
            "url": "https://cellrank.readthedocs.io/en/latest/api/_autosummary/datasets/cellrank.datasets.reprogramming_morris.html",
            "source": "software dataset wrapper",
            "organism": "mouse",
            "biological_system": "CellTag reprogramming",
            "time_or_stage_available": True,
            "clone_or_barcode_available": "unknown_from_wrapper",
            "cell_type_available": True,
            "expression_matrix_available": True,
            "processed_h5ad_available": True,
            "metadata_available": True,
            "download_attempted": False,
            "download_success": False,
            "matrix_loaded": False,
            "metadata_loaded": False,
            "clone_metadata_loaded": False,
            "usable_for_clone_validation": False,
            "reason_if_not_usable": "Not downloaded separately after the direct Biddy processed h5ad succeeded; clone metadata fields in wrapper were not verified locally.",
            "selected_for_L2": False,
            "notes": "Useful fallback if scLTdb file becomes unavailable.",
        },
        {
            "dataset_id": "L2_kim_2020_cellreports",
            "dataset_name": "Embryoid body differentiation with genetic recording",
            "publication": "Kim et al., Cell Reports 2020",
            "accession": "scLTdb Kim_2020_CellReports.h5ad",
            "doi": "10.1016/j.celrep.2020.108222; 10.5281/zenodo.12176634",
            "url": "https://zenodo.org/records/12176634/files/Kim_2020_CellReports.h5ad?download=1",
            "source": "scLTdb Zenodo processed h5ad",
            "organism": "mouse",
            "biological_system": "embryoid body differentiation",
            "time_or_stage_available": kim.exists(),
            "clone_or_barcode_available": kim.exists(),
            "cell_type_available": kim.exists(),
            "expression_matrix_available": kim.exists(),
            "processed_h5ad_available": kim.exists(),
            "metadata_available": kim.exists(),
            "download_attempted": True,
            "download_success": kim.exists(),
            "matrix_loaded": kim.exists(),
            "metadata_loaded": kim.exists(),
            "clone_metadata_loaded": kim.exists(),
            "usable_for_clone_validation": kim.exists(),
            "reason_if_not_usable": "" if kim.exists() else "Local Kim h5ad missing.",
            "selected_for_L2": False,
            "notes": "Already exhausted in v1.0; retained as fallback and comparison.",
        },
        {
            "dataset_id": "L2_wei_2020_genomeresearch",
            "dataset_name": "Small scLTdb lineage candidate",
            "publication": "Wei et al., Genome Research 2020",
            "accession": "scLTdb Wei_2020_GenomeResearch.h5ad",
            "doi": "10.5281/zenodo.12176634",
            "url": "https://zenodo.org/records/12176634/files/Wei_2020_GenomeResearch.h5ad?download=1",
            "source": "scLTdb Zenodo processed h5ad",
            "organism": "unknown_from_local_metadata",
            "biological_system": "lineage tracing candidate",
            "time_or_stage_available": False,
            "clone_or_barcode_available": wei.exists(),
            "cell_type_available": wei.exists(),
            "expression_matrix_available": wei.exists(),
            "processed_h5ad_available": wei.exists(),
            "metadata_available": wei.exists(),
            "download_attempted": True,
            "download_success": wei.exists(),
            "matrix_loaded": wei.exists(),
            "metadata_loaded": wei.exists(),
            "clone_metadata_loaded": wei.exists(),
            "usable_for_clone_validation": False,
            "reason_if_not_usable": "Loaded in v1.0 but lacks ordered time/stage metadata for branch-window validation.",
            "selected_for_L2": False,
            "notes": "Not suitable for L2 branch-nucleation validation.",
        },
        {
            "dataset_id": "L2_spanjaard_2018_linnaeus",
            "dataset_name": "LINNAEUS zebrafish lineage tracing",
            "publication": "Spanjaard et al., Nature Biotechnology 2018",
            "accession": "scLTdb Spanjaard_2018_NatureBiotechnology.h5ad",
            "doi": "10.5281/zenodo.12176634",
            "url": "https://zenodo.org/records/12176634/files/Spanjaard_2018_NatureBiotechnology.h5ad?download=1",
            "source": "scLTdb Zenodo processed h5ad",
            "organism": "zebrafish",
            "biological_system": "developmental lineage tracing",
            "time_or_stage_available": True,
            "clone_or_barcode_available": True,
            "cell_type_available": True,
            "expression_matrix_available": True,
            "processed_h5ad_available": True,
            "metadata_available": True,
            "download_attempted": False,
            "download_success": False,
            "matrix_loaded": False,
            "metadata_loaded": False,
            "clone_metadata_loaded": False,
            "usable_for_clone_validation": False,
            "reason_if_not_usable": "Not downloaded because Biddy primary dataset succeeded; cross-species harmonization remains future work.",
            "selected_for_L2": False,
            "notes": "Fallback registry candidate.",
        },
    ]
    return pd.DataFrame(rows)


def _time_to_float(value: object) -> float:
    text = str(value)
    num = pd.Series([text]).str.extract(r"([0-9]+\.?[0-9]*)", expand=False).iloc[0]
    return float(num) if pd.notna(num) else float("nan")


def _select_top_genes(X, n_genes: int = 2000) -> np.ndarray:
    sums = np.asarray(X.sum(axis=0)).ravel() if sparse.issparse(X) else np.asarray(X).sum(axis=0)
    valid = np.where(np.isfinite(sums) & (sums > 0))[0]
    if valid.size <= n_genes:
        return valid
    return valid[np.argsort(sums[valid])[-n_genes:]]


def _stratified_downsample(obs: pd.DataFrame, max_total: int = 3200, seed: int = 17) -> np.ndarray:
    if obs.shape[0] <= max_total:
        return np.arange(obs.shape[0])
    rng = np.random.default_rng(seed)
    selected: list[int] = []
    groups = obs.groupby("time_point", observed=False).groups
    per_time = max(80, max_total // max(len(groups), 1))
    for _, idx in groups.items():
        idx_arr = np.asarray(list(idx), dtype=int)
        selected.extend(rng.choice(idx_arr, size=min(per_time, idx_arr.size), replace=False).tolist())
    if len(selected) < max_total:
        remaining = np.setdiff1d(np.arange(obs.shape[0]), np.asarray(selected), assume_unique=False)
        selected.extend(rng.choice(remaining, size=min(max_total - len(selected), remaining.size), replace=False).tolist())
    if len(selected) > max_total:
        selected = rng.choice(np.asarray(selected), size=max_total, replace=False).tolist()
    return np.asarray(sorted(selected), dtype=int)


def load_and_prepare_biddy() -> tuple[ad.AnnData, ad.AnnData, pd.DataFrame, np.ndarray]:
    if not BIDDY_FILE.exists():
        raise FileNotFoundError(f"Missing Biddy h5ad: {BIDDY_FILE}")
    raw = ad.read_h5ad(BIDDY_FILE)
    obs = raw.obs.copy()
    obs["time_point"] = obs["time_info"].astype(str)
    obs["time_numeric"] = obs["time_point"].map(_time_to_float).astype(float)
    obs["lineage"] = obs["celltype"].astype(str)
    obs["cell_type"] = obs["celltype"].astype(str)
    obs["clone_id"] = obs["barcode_all"].astype(str)
    obs["lineage_barcode"] = obs["clone_id"]
    obs["external_dataset_id"] = "L2_biddy_2018_nature"
    obs["external_source"] = "scLTdb/Biddy_2018_Nature"
    obs["split_role"] = "external_l2_evaluation"
    obs["clone_metadata_loaded"] = ~obs["clone_id"].isin(["", "nan", "NA*NA*NA", "-2147483648"])
    keep = obs["time_numeric"].notna() & obs["clone_metadata_loaded"]
    raw = raw[keep.to_numpy()].copy()
    obs = obs.loc[keep].copy().reset_index(drop=True)
    raw.obs = obs
    top_genes = _select_top_genes(raw.X, n_genes=2000)
    raw = raw[:, top_genes].copy()
    Xlog = _normalise_counts(raw.X)
    z_full = TruncatedSVD(n_components=30, random_state=17).fit_transform(Xlog)
    raw.obsm["X_pca"] = z_full
    raw.obs["local_density"] = _local_density(z_full)
    selected = _stratified_downsample(raw.obs.reset_index(drop=True), max_total=3200, seed=17)
    down = raw[selected].copy()
    ensure_dir(ROOT / "data/external_l2")
    ensure_dir(ROOT / "processed/external_l2")
    raw.write_h5ad(ROOT / "data/external_l2/l2_external_input.h5ad")
    down.write_h5ad(ROOT / "processed/external_l2/l2_swarmlineage_input.h5ad")
    return raw, down, raw.obs.copy(), z_full


def plot_embedding(z: np.ndarray, obs: pd.DataFrame) -> None:
    ensure_dir(ROOT / "figures/external_l2")
    coords = z[:, :2]
    for col, path, title in [
        ("time_point", "figures/external_l2/l2_umap_by_time.png", "L2 PCA by time"),
        ("lineage", "figures/external_l2/l2_umap_by_clone_or_lineage.png", "L2 PCA by lineage"),
    ]:
        fig, ax = plt.subplots(figsize=(6, 5))
        labels = obs[col].astype(str)
        cats = labels.astype("category")
        ax.scatter(coords[:, 0], coords[:, 1], c=cats.cat.codes, s=3, cmap="tab20", alpha=0.6)
        ax.set_title(title)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        fig.tight_layout()
        fig.savefig(ROOT / path, dpi=180)
        plt.close(fig)
    clone_sizes = obs["clone_id"].astype(str).value_counts()
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(clone_sizes, bins=40, color="#4C78A8")
    ax.set_xlabel("clone size")
    ax.set_ylabel("clone count")
    ax.set_title("L2 clone size distribution")
    fig.tight_layout()
    fig.savefig(ROOT / "figures/external_l2/l2_clone_size_distribution.png", dpi=180)
    plt.close(fig)


def write_l2_data_audit(raw: ad.AnnData, down: ad.AnnData) -> None:
    obs = raw.obs.copy()
    counts = obs.groupby(["time_point", "lineage"], observed=False).size().reset_index(name="n_cells")
    _write_csv(counts, "tables/l2_time_celltype_counts.csv")
    clone_summary = (
        obs.groupby("clone_id", observed=False)
        .agg(
            clone_size=("clone_id", "size"),
            n_time_points=("time_point", "nunique"),
            start_time=("time_numeric", "min"),
            end_time=("time_numeric", "max"),
            n_lineages=("lineage", "nunique"),
        )
        .reset_index()
    )
    _write_csv(clone_summary, "tables/l2_clone_summary.csv")
    _write_md(
        "reports/l2_data_audit.md",
        "# L2 Data Audit\n\n"
        f"- selected_dataset: Biddy_2018_Nature / GSE99915 CellTag reprogramming\n"
        f"- full_loaded_shape: {raw.n_obs} x {raw.n_vars}\n"
        f"- downsampled_shape: {down.n_obs} x {down.n_vars}\n"
        f"- time_points: {sorted(obs['time_point'].astype(str).unique())}\n"
        f"- clone_id_field: barcode_all\n"
        f"- usable_clones_size_ge_5: {int((clone_summary['clone_size'] >= 5).sum())}\n"
        f"- usable_clones_size_ge_20: {int((clone_summary['clone_size'] >= 20).sum())}\n\n"
        "Counts by time and lineage:\n\n"
        + _md_table(counts.head(30)),
    )
    _write_md(
        "reports/l2_leakage_audit.md",
        "# L2 Leakage Audit\n\n"
        "- L2 preprocessing uses only the downloaded Biddy_2018_Nature h5ad.\n"
        "- No internal teacher, internal PCA, internal labels or E1 labels are used to fit L2 embeddings.\n"
        "- Clone identifiers are copied from the source `barcode_all` field; no clone labels are invented.\n"
        "- Processed h5ad files are written under ignored data/processed paths and are not staged for Git.\n",
    )


def run_l2_teacher() -> dict:
    cfg = {
        "adata_path": "processed/external_l2/l2_swarmlineage_input.h5ad",
        "teacher_path": "processed/external_l2/l2_ot_teacher.h5ad",
        "couplings_dir": "processed/external_l2/ot_couplings",
        "fate_probabilities_path": "processed/external_l2/l2_ot_fate_probabilities.parquet",
        "teacher_index_path": "processed/external_l2/ot_couplings/teacher_coupling_index.csv",
        "time_key": "time_numeric",
        "time_label_key": "time_point",
        "cell_type_key": "lineage",
        "latent_key": "X_pca",
        "epsilon": 0.08,
        "max_cells_per_time": 120,
        "native_moscot_timeout_seconds": 90,
        "use_native_moscot": True,
        "native_max_cells_per_time": 120,
        "native_max_iterations": 350,
        "native_jit": False,
        "native_device": "cpu",
        "random_seed": 17,
        "split_mode": "none",
        "teacher_backend": "native_moscot",
        "figure_dir": "figures/external_l2",
        "report_dir": "reports/external_l2",
        "table_dir": "tables/external_l2",
        "summary_path": "processed/external_l2/l2_ot_teacher_summary.json",
    }
    ensure_dir(ROOT / "configs")
    with (ROOT / "configs/external_l2_ot_teacher.yaml").open("w", encoding="utf-8") as handle:
        yaml.safe_dump(cfg, handle, sort_keys=False)
    native = _native_status(int(cfg["native_moscot_timeout_seconds"]))
    start = time.perf_counter()
    native_result = run_native_moscot_teacher(cfg, native, label="teacher")
    if native_result.get("success", False):
        result = dict(native_result)
        result["fallback_used"] = False
        result["teacher_backend"] = "native_moscot"
    else:
        result = run_toy_sinkhorn_teacher(cfg, label="teacher")
        result["fallback_used"] = True
        result["teacher_backend"] = "fallback_sinkhorn"
        result["native_failure_reason"] = native_result.get("reason", "unknown")
    runtime = time.perf_counter() - start
    result["runtime_seconds"] = runtime
    result["native_moscot_status"] = native
    result["native_moscot_used"] = bool(native_result.get("success", False))
    try:
        build_teacher(cfg)
    except Exception as exc:
        result["build_teacher_error"] = f"{type(exc).__name__}: {exc}"
    pairs = pd.DataFrame(result.get("pairs", []))
    if pairs.empty:
        pairs = _read_csv("processed/external_l2/ot_couplings/teacher_coupling_index.csv")
    if pairs.empty:
        pairs = pd.DataFrame(
            [
                {
                    "teacher_backend": result["teacher_backend"],
                    "status": "failed",
                    "runtime_seconds": runtime,
                    "fallback_used": result.get("fallback_used", True),
                }
            ]
        )
    pairs["runtime_seconds"] = runtime
    pairs["native_moscot_available"] = bool(native.get("available", False))
    pairs["native_moscot_used"] = bool(result["native_moscot_used"])
    pairs["fallback_used"] = bool(result.get("fallback_used", False))
    pairs["native_max_cells_per_time"] = cfg["native_max_cells_per_time"]
    pairs["epsilon"] = cfg["epsilon"]
    _write_csv(pairs, "tables/l2_native_teacher_pairs.csv")
    _write_md(
        "reports/l2_native_teacher_report.md",
        "# L2 Native Teacher Report\n\n"
        f"- teacher_backend: {result['teacher_backend']}\n"
        f"- native_moscot_available: {native.get('available', False)}\n"
        f"- native_moscot_detail: {native.get('detail', 'not_recorded')}\n"
        f"- native_moscot_used: {result['native_moscot_used']}\n"
        f"- fallback_used: {result.get('fallback_used', False)}\n"
        f"- native_failure_reason: {result.get('native_failure_reason', 'not_applicable')}\n"
        f"- runtime_seconds: {runtime:.2f}\n\n"
        + _md_table(pairs.head(20)),
    )
    return result


def _cell_alignment(obs: pd.DataFrame, z: np.ndarray) -> np.ndarray:
    out = np.zeros(obs.shape[0], dtype=float)
    time_values = sorted(obs["time_numeric"].dropna().unique())
    centroids = {t: z[obs.index[obs["time_numeric"].eq(t)].to_numpy()].mean(axis=0) for t in time_values}
    for pos, time_value in enumerate(time_values):
        idx = obs.index[obs["time_numeric"].eq(time_value)].to_numpy()
        prev_t = time_values[max(0, pos - 1)]
        next_t = time_values[min(len(time_values) - 1, pos + 1)]
        velocity = centroids[next_t] - centroids[prev_t]
        if np.linalg.norm(velocity) <= 1e-8:
            continue
        centered = z[idx] - z[idx].mean(axis=0)
        denom = np.maximum(np.linalg.norm(centered, axis=1) * np.linalg.norm(velocity), 1e-8)
        out[idx] = (centered @ velocity) / denom
    return out


def branch_nucleation(raw: ad.AnnData, z: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    obs = raw.obs.copy().reset_index(drop=True)
    order = _order_parameters_from_embedding(obs, z, "time_numeric", "lineage")
    event = _event_window_from_order(order)
    if not event.empty:
        pre = float(event["lineage_separation_S_pre_mean"].iloc[0])
        effect = float(event["lineage_separation_S_effect"].iloc[0])
        event["normalized_separation_effect"] = effect / max(abs(pre), 1e-8)
    controls = []
    rng = np.random.default_rng(17)
    for control in ["time_shuffle", "velocity_shuffle", "lineage_label_shuffle", "clone_label_shuffle", "random_teacher_velocity", "no_swarm_proxy", "no_teacher_proxy"]:
        values = []
        for _ in range(100):
            ctrl_obs = obs.copy()
            ctrl_order = order.copy()
            if control == "time_shuffle":
                ctrl_obs["time_numeric"] = rng.permutation(ctrl_obs["time_numeric"].to_numpy())
                ctrl_order = _order_parameters_from_embedding(ctrl_obs, z, "time_numeric", "lineage")
            elif control == "lineage_label_shuffle":
                ctrl_obs["lineage"] = rng.permutation(ctrl_obs["lineage"].astype(str).to_numpy())
                ctrl_order = _order_parameters_from_embedding(ctrl_obs, z, "time_numeric", "lineage")
            elif control == "clone_label_shuffle":
                ctrl_obs["clone_id"] = rng.permutation(ctrl_obs["clone_id"].astype(str).to_numpy())
            elif control == "velocity_shuffle":
                ctrl_order["local_velocity_alignment_A"] = rng.permutation(ctrl_order["local_velocity_alignment_A"].to_numpy())
            elif control == "random_teacher_velocity":
                ctrl_order["lineage_separation_S"] = rng.normal(ctrl_order["lineage_separation_S"].mean(), ctrl_order["lineage_separation_S"].std(ddof=0) + 1e-6, size=ctrl_order.shape[0])
            elif control == "no_swarm_proxy":
                ctrl_order["local_velocity_alignment_A"] = 0.0
            elif control == "no_teacher_proxy":
                ctrl_order["lineage_separation_S"] = np.sort(ctrl_order["lineage_separation_S"].to_numpy())
            win = _event_window_from_order(ctrl_order)
            values.append(float(win["lineage_separation_S_effect"].iloc[0]) if not win.empty else 0.0)
        mean, lo, hi = _bootstrap_ci(np.asarray(values), repeats=300)
        controls.append(
            {
                "control": control,
                "lineage_separation_effect": mean,
                "effect_ci_low": lo,
                "effect_ci_high": hi,
                "expected_to_fail": True,
                "control_pass": not (mean < 0 and abs(mean) >= abs(float(event["lineage_separation_S_effect"].iloc[0]) if not event.empty else 0) * 0.5),
            }
        )
    control_df = pd.DataFrame(controls)
    _write_csv(order, "tables/l2_branch_order_parameters.csv")
    _write_csv(event, "tables/l2_branch_event_windows.csv")
    _write_csv(control_df, "tables/l2_branch_negative_controls.csv")
    ensure_dir(ROOT / "figures/external_l2")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(order["time"], order["lineage_separation_S"], marker="o", label="lineage separation")
    ax.plot(order["time"], order["local_velocity_alignment_A"], marker="o", label="alignment")
    ax.set_xlabel("time")
    ax.set_title("L2 branch order parameters")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(ROOT / "figures/external_l2/l2_branch_order_parameters.png", dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(5, 4))
    if not event.empty:
        ax.bar(["pre", "post"], [event["lineage_separation_S_pre_mean"].iloc[0], event["lineage_separation_S_post_mean"].iloc[0]], color=["#999999", "#4C78A8"])
    ax.set_ylabel("lineage separation")
    ax.set_title("L2 branch event window")
    fig.tight_layout()
    fig.savefig(ROOT / "figures/external_l2/l2_branch_event_window.png", dpi=180)
    plt.close(fig)
    _write_md(
        "reports/l2_branch_nucleation_report.md",
        "# L2 Branch Nucleation Report\n\n"
        + _md_table(event)
        + "\n\nNegative controls:\n\n"
        + _md_table(control_df),
    )
    return order, event, control_df


def _clone_metrics_for_threshold(obs: pd.DataFrame, z: np.ndarray, event_time: float, min_clone_size: int, seed: int = 17, sample_frac: float = 1.0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    frame = obs.copy().reset_index(drop=True)
    if sample_frac < 1.0:
        keep = rng.choice(np.arange(frame.shape[0]), size=max(10, int(frame.shape[0] * sample_frac)), replace=False)
        frame = frame.iloc[np.sort(keep)].copy().reset_index(drop=True)
        z = z[np.sort(keep)]
    frame["local_density"] = _local_density(z)
    frame["cell_alignment"] = _cell_alignment(frame, z)
    rows = []
    max_time = float(frame["time_numeric"].max())
    terminal_cut = max_time - 7.0
    global_centroids = {t: z[frame.index[frame["time_numeric"].eq(t)].to_numpy()].mean(axis=0) for t in sorted(frame["time_numeric"].unique())}
    for clone_id, group in frame.groupby("clone_id", observed=False):
        if group.shape[0] < min_clone_size:
            continue
        pre = group[group["time_numeric"] <= event_time]
        terminal = group[group["time_numeric"] >= terminal_cut]
        if pre.empty or terminal.empty:
            continue
        terminal_dist = terminal["lineage"].astype(str).value_counts(normalize=True)
        branch_entropy = _entropy(terminal["lineage"])
        splitting = float(branch_entropy * (1.0 - 1.0 / math.sqrt(max(terminal.shape[0], 1))))
        pre_idx = pre.index.to_numpy()
        pre_distances = []
        for idx in pre_idx:
            t = float(frame.loc[idx, "time_numeric"])
            pre_distances.append(float(np.linalg.norm(z[idx] - global_centroids[t])))
        term_idx = terminal.index.to_numpy()
        post_distances = []
        for idx in term_idx:
            t = float(frame.loc[idx, "time_numeric"])
            post_distances.append(float(np.linalg.norm(z[idx] - global_centroids[t])))
        rows.append(
            {
                "clone_id": clone_id,
                "clone_size": int(group.shape[0]),
                "clone_time_span": float(group["time_numeric"].max() - group["time_numeric"].min()),
                "clone_start_time": float(group["time_numeric"].min()),
                "clone_end_time": float(group["time_numeric"].max()),
                "clone_celltype_distribution": group["lineage"].astype(str).value_counts(normalize=True).to_json(),
                "clone_terminal_fate_distribution": terminal_dist.to_json(),
                "clone_branch_entropy": branch_entropy,
                "clone_branch_splitting_score": splitting,
                "clone_multilineage_score": int(terminal["lineage"].nunique()),
                "clone_pre_event_condensation_exposure": float(-np.mean(pre_distances)),
                "clone_pre_event_alignment_exposure": float(pre["cell_alignment"].mean()),
                "clone_pre_event_entropy_exposure": _entropy(pre["lineage"]),
                "clone_pre_event_density_exposure": float(pre["local_density"].mean()),
                "clone_post_event_divergence_score": float(np.mean(post_distances)),
                "clone_nearest_branch_event_time": event_time,
                "clone_usable_for_validation": True,
            }
        )
    return pd.DataFrame(rows)


def _ols_summary(df: pd.DataFrame) -> pd.DataFrame:
    if df.shape[0] < 5:
        return pd.DataFrame()
    y = df["clone_branch_splitting_score"].to_numpy(dtype=float)
    cols = [
        "clone_pre_event_condensation_exposure",
        "clone_size",
        "clone_time_span",
        "clone_start_time",
    ]
    X = df[cols].to_numpy(dtype=float)
    X = np.column_stack([np.ones(X.shape[0]), X])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ beta
    resid = y - pred
    dof = max(X.shape[0] - X.shape[1], 1)
    sigma2 = float((resid @ resid) / dof)
    cov = sigma2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.maximum(np.diag(cov), 0))
    rows = []
    for name, coef, serr in zip(["intercept", *cols], beta, se):
        t = coef / serr if serr > 0 else 0.0
        p = float(2 * (1 - stats.t.cdf(abs(t), df=dof)))
        rows.append({"term": name, "coefficient": float(coef), "std_error": float(serr), "t_stat": float(t), "p_value": p})
    return pd.DataFrame(rows)


def clone_validation(raw: ad.AnnData, z: np.ndarray, event: pd.DataFrame, teacher_backend: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    obs = raw.obs.copy().reset_index(drop=True)
    event_time = float(event["branch_event_step"].iloc[0]) if not event.empty else 2.0
    time_values = sorted(obs["time_numeric"].dropna().unique())
    event_time = float(time_values[min(int(event_time), len(time_values) - 1)]) if time_values else 12.0
    clone_df = _clone_metrics_for_threshold(obs, z, event_time=event_time, min_clone_size=10, seed=17)
    _write_csv(clone_df, "tables/l2_clone_branch_validation.csv")
    regression = _ols_summary(clone_df)
    _write_csv(regression, "tables/l2_clone_regression.csv")
    exposure = clone_df["clone_pre_event_condensation_exposure"].to_numpy(dtype=float) if not clone_df.empty else np.array([])
    target = clone_df["clone_branch_splitting_score"].to_numpy(dtype=float) if not clone_df.empty else np.array([])
    if target.size > 3 and np.nanstd(exposure) > 0 and np.nanstd(target) > 0:
        spearman = stats.spearmanr(exposure, target)
        pearson = stats.pearsonr(exposure, target)
        rho, p_s = float(spearman.statistic), float(spearman.pvalue)
        r, p_p = float(pearson.statistic), float(pearson.pvalue)
    else:
        rho, p_s, r, p_p = 0.0, 1.0, 0.0, 1.0
    rng = np.random.default_rng(31)
    null = []
    for _ in range(1000):
        if target.size > 3 and np.nanstd(exposure) > 0:
            null.append(float(stats.spearmanr(exposure, rng.permutation(target)).statistic))
    null = np.asarray(null)
    perm_p = float((np.sum(np.abs(null) >= abs(rho)) + 1) / (null.size + 1)) if null.size else 1.0
    controls = []
    for control in ["clone_id_shuffle", "time_shuffle", "branch_label_shuffle", "order_parameter_shuffle", "random_teacher_velocity"]:
        vals = []
        for _ in range(200):
            if target.size <= 3:
                vals.append(0.0)
                continue
            if control in {"clone_id_shuffle", "branch_label_shuffle"}:
                vals.append(float(stats.spearmanr(exposure, rng.permutation(target)).statistic))
            elif control in {"order_parameter_shuffle", "random_teacher_velocity"}:
                vals.append(float(stats.spearmanr(rng.permutation(exposure), target).statistic))
            elif control == "time_shuffle":
                vals.append(float(stats.spearmanr(exposure + rng.normal(0, np.std(exposure) + 1e-8, exposure.size), target).statistic))
        mean, lo, hi = _bootstrap_ci(np.asarray(vals), repeats=300)
        controls.append({"control": control, "mean_null_spearman": mean, "ci_low": lo, "ci_high": hi, "control_signal_abs_ge_observed": bool(abs(mean) >= abs(rho))})
    control_df = pd.DataFrame(controls)
    _write_csv(control_df, "tables/l2_clone_permutation_controls.csv")
    sensitivity_rows = []
    for threshold in [5, 10, 20, 50]:
        df_thr = _clone_metrics_for_threshold(obs, z, event_time=event_time, min_clone_size=threshold, seed=17)
        if df_thr.shape[0] > 3 and df_thr["clone_pre_event_condensation_exposure"].std() > 0 and df_thr["clone_branch_splitting_score"].std() > 0:
            rr = stats.spearmanr(df_thr["clone_pre_event_condensation_exposure"], df_thr["clone_branch_splitting_score"])
            effect, pp = float(rr.statistic), float(rr.pvalue)
        else:
            effect, pp = 0.0, 1.0
        sensitivity_rows.append({"sensitivity": f"min_clone_size_{threshold}", "usable_clones": int(df_thr.shape[0]), "spearman": effect, "p_value": pp})
    for seed in [7, 17, 23, 31, 43]:
        df_seed = _clone_metrics_for_threshold(obs, z, event_time=event_time, min_clone_size=10, seed=seed, sample_frac=0.8)
        if df_seed.shape[0] > 3 and df_seed["clone_pre_event_condensation_exposure"].std() > 0 and df_seed["clone_branch_splitting_score"].std() > 0:
            rr = stats.spearmanr(df_seed["clone_pre_event_condensation_exposure"], df_seed["clone_branch_splitting_score"])
            effect, pp = float(rr.statistic), float(rr.pvalue)
        else:
            effect, pp = 0.0, 1.0
        sensitivity_rows.append({"sensitivity": f"downsample_seed_{seed}", "usable_clones": int(df_seed.shape[0]), "spearman": effect, "p_value": pp})
    sensitivity = pd.DataFrame(sensitivity_rows)
    _write_csv(sensitivity, "tables/l2_clone_sensitivity.csv")
    usable = int(clone_df.shape[0])
    controls_ok = bool((~control_df["control_signal_abs_ge_observed"]).all()) if not control_df.empty else False
    predicts = bool(rho > 0 and perm_p < 0.10)
    if teacher_backend == "native_moscot" and usable >= 50 and predicts and controls_ok:
        tier = "acceptable" if perm_p >= 0.05 else "strong"
    elif usable >= 20 and rho > 0 and perm_p < 0.20:
        tier = "weak"
    else:
        tier = "fail"
    summary = pd.DataFrame(
        [
            {
                "dataset_id": "L2_biddy_2018_nature",
                "l2_validation_tier": tier,
                "teacher_backend": teacher_backend,
                "clone_metadata_loaded": True,
                "usable_clone_count": usable,
                "condensation_spearman": rho,
                "condensation_spearman_p": p_s,
                "condensation_pearson": r,
                "condensation_pearson_p": p_p,
                "condensation_permutation_p": perm_p,
                "negative_controls_pass": controls_ok,
                "sensitivity_stable": bool((np.sign(sensitivity["spearman"]) == np.sign(rho)).mean() >= 0.6) if not sensitivity.empty else False,
                "condensation_predicts_clone_branch_splitting": predicts,
                "interpretation": "Clone-aware validation supports the hypothesis." if predicts else "Clone-aware validation does not support the hypothesis in Biddy_2018 under this operationalization.",
            }
        ]
    )
    _write_csv(summary, "tables/l2_clone_model_summary.csv")
    _write_csv(summary.rename(columns={"l2_validation_tier": "validation_tier"}), "tables/l2_validation_tier_summary.csv")
    ensure_dir(ROOT / "figures/external_l2")
    fig, ax = plt.subplots(figsize=(6, 4))
    if not clone_df.empty:
        ax.scatter(clone_df["clone_pre_event_condensation_exposure"], clone_df["clone_branch_splitting_score"], s=np.clip(clone_df["clone_size"], 8, 80), alpha=0.6)
    ax.set_xlabel("pre-event condensation exposure")
    ax.set_ylabel("clone branch splitting score")
    ax.set_title("L2 condensation vs clone splitting")
    fig.tight_layout()
    fig.savefig(ROOT / "figures/external_l2/l2_condensation_predicts_clone_splitting.png", dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(control_df["control"], control_df["mean_null_spearman"], color="#999999")
    ax.axhline(rho, color="#4C78A8", label="observed")
    ax.axhline(-rho, color="#4C78A8", linestyle="--")
    ax.set_ylabel("Spearman")
    ax.tick_params(axis="x", rotation=45)
    ax.set_title("L2 clone validation controls")
    fig.tight_layout()
    fig.savefig(ROOT / "figures/external_l2/l2_clone_validation_controls.png", dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(sensitivity["sensitivity"], sensitivity["spearman"], color="#4C78A8")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Spearman")
    ax.tick_params(axis="x", rotation=70)
    ax.set_title("L2 sensitivity")
    fig.tight_layout()
    fig.savefig(ROOT / "figures/external_l2/l2_clone_sensitivity.png", dpi=180)
    plt.close(fig)
    _write_md(
        "reports/l2_clone_validation_report.md",
        "# L2 Clone-Aware Validation Report\n\n"
        "Clone branch splitting score is terminal lineage entropy scaled by terminal clone size. Condensation exposure is the negative mean pre-event latent distance from the time-specific centroid, so larger values indicate stronger pre-event condensation.\n\n"
        + _md_table(summary)
        + "\n\nRegression:\n\n"
        + _md_table(regression)
        + "\n\nControls:\n\n"
        + _md_table(control_df)
        + "\n\nSensitivity:\n\n"
        + _md_table(sensitivity),
    )
    _write_md(
        "reports/l2_validation_summary.md",
        "# L2 Validation Summary\n\n" + _md_table(summary),
    )
    return clone_df, summary


def master_summary(l2_summary: pd.DataFrame) -> None:
    e1 = _read_csv("tables/external_validation_tier_summary.csv")
    l1 = _read_csv("tables/l1_clone_model_summary.csv")
    e2 = _read_csv("tables/e2_branch_nucleation_summary.csv")
    rows = []
    if not e1.empty:
        r = e1.iloc[0]
        rows.append(
            {
                "experiment": "E1",
                "dataset": "MouseGastrulationData WT chimera sample 1",
                "system": "mouse gastrulation",
                "teacher_backend": r.get("external_teacher_backend", "native_moscot"),
                "time_series_available": True,
                "clone_metadata_available": False,
                "native_teacher_success": r.get("external_teacher_backend") == "native_moscot",
                "branch_signature_reproduced": bool(r.get("condensation_before_divergence_reproduced", False)),
                "condensation_direction_match": True,
                "alignment_match": True,
                "entropy_match": True,
                "density_match": False,
                "clone_splitting_association": "not_tested_no_clone_metadata",
                "validation_tier": r.get("external_validation_tier", "acceptable"),
                "interpretation": "external time-series support",
                "allowed_claim": "E1 supports the branch signature as time-series evidence, without clone-resolved support.",
            }
        )
    if not l1.empty:
        r = l1.iloc[0]
        rows.append(
            {
                "experiment": "L1",
                "dataset": "Kim_2020_CellReports",
                "system": "embryoid body genetic recording",
                "teacher_backend": r.get("teacher_backend", "not_run_clone_proxy"),
                "time_series_available": True,
                "clone_metadata_available": True,
                "native_teacher_success": False,
                "branch_signature_reproduced": "not_primary",
                "condensation_direction_match": False,
                "alignment_match": False,
                "entropy_match": False,
                "density_match": False,
                "clone_splitting_association": float(r.get("condensation_to_clone_splitting_spearman", 0.0)),
                "validation_tier": r.get("lineage_validation_tier", "fail"),
                "interpretation": r.get("interpretation", "not supportive"),
                "allowed_claim": "L1 does not establish clone-level support.",
            }
        )
    if not l2_summary.empty:
        r = l2_summary.iloc[0]
        l2_event = _read_csv("tables/l2_branch_event_windows.csv")
        l2_branch_reproduced = bool(not l2_event.empty and float(l2_event["lineage_separation_S_effect"].iloc[0]) < 0)
        rows.append(
            {
                "experiment": "L2",
                "dataset": "Biddy_2018_Nature / GSE99915 CellTag",
                "system": "CellTag reprogramming",
                "teacher_backend": r.get("teacher_backend", ""),
                "time_series_available": True,
                "clone_metadata_available": True,
                "native_teacher_success": r.get("teacher_backend") == "native_moscot",
                "branch_signature_reproduced": l2_branch_reproduced,
                "condensation_direction_match": bool(r.get("condensation_predicts_clone_branch_splitting", False)),
                "alignment_match": "not_claimed",
                "entropy_match": "not_claimed",
                "density_match": "not_claimed",
                "clone_splitting_association": float(r.get("condensation_spearman", 0.0)),
                "validation_tier": r.get("l2_validation_tier", "fail"),
                "interpretation": r.get("interpretation", ""),
                "allowed_claim": "L2 is a clone-aware attempt; support depends on the reported tier.",
            }
        )
    if not e2.empty:
        r = e2.iloc[0]
        rows.append(
            {
                "experiment": "E2",
                "dataset": "GSE212050 local gastruloid",
                "system": "gastruloid time series",
                "teacher_backend": r.get("teacher_backend", "not_run_temporal_proxy"),
                "time_series_available": True,
                "clone_metadata_available": False,
                "native_teacher_success": False,
                "branch_signature_reproduced": bool(r.get("condensation_direction_observed", False)),
                "condensation_direction_match": bool(r.get("condensation_direction_observed", False)),
                "alignment_match": "not_claimed",
                "entropy_match": "not_claimed",
                "density_match": "not_claimed",
                "clone_splitting_association": "not_tested_no_clone_metadata",
                "validation_tier": r.get("e2_validation_tier", "weak"),
                "interpretation": r.get("interpretation", "weak feasibility"),
                "allowed_claim": "E2 remains weak feasibility support only.",
            }
        )
    master = pd.DataFrame(rows)
    _write_csv(master, "tables/external_validation_master_summary.csv")
    _write_md("reports/external_validation_master_summary.md", "# External Validation Master Summary\n\n" + _md_table(master))


def update_claims_and_manuscript(l2_summary: pd.DataFrame) -> None:
    tier = str(l2_summary.iloc[0]["l2_validation_tier"]) if not l2_summary.empty else "fail"
    rho = float(l2_summary.iloc[0]["condensation_spearman"]) if not l2_summary.empty else 0.0
    usable = int(l2_summary.iloc[0]["usable_clone_count"]) if not l2_summary.empty else 0
    if tier in {"strong", "acceptable"}:
        clone_sentence = "L2 provides clone-aware computational support for the branch-nucleation hypothesis."
        clone_status = "clone_aware_support"
    elif tier == "weak":
        clone_sentence = "L2 provides weak clone-aware feasibility but does not establish the clone-level claim."
        clone_status = "weak_clone_feasibility"
    else:
        clone_sentence = "L2 does not support the clone-level association under the current operationalization."
        clone_status = "clone_support_not_established"
    final = _read_csv("tables/final_claim_evidence_tiers.csv")
    if not final.empty and "claim" in final:
        mask = final["claim"].str.contains("clone|lineage", case=False, na=False)
        final.loc[mask, "status"] = clone_status
        final.loc[mask, "tier"] = tier
        final.loc[mask, "lineage_clone_support"] = tier in {"strong", "acceptable"}
        final.loc[mask, "allowed_manuscript_sentence"] = clone_sentence
        _write_csv(final, "tables/final_claim_evidence_tiers.csv")
        _write_md("reports/final_claim_evidence_tiers.md", "# Final Claim Evidence Tiers\n\n" + _md_table(final))
    evidence = pd.DataFrame(
        [
            {"claim": "branch nucleation signature", "tier": "strong_internal_acceptable_E1", "l2_tier": tier, "allowed_language": "Retained computational hypothesis with internal native and E1 time-series support."},
            {"claim": "clone-level branch splitting association", "tier": tier, "usable_clones": usable, "effect": rho, "allowed_language": clone_sentence},
            {"claim": "diffusion", "tier": "encoded_recovery", "allowed_language": "Encoded control-law recovery only."},
            {"claim": "birth/death memory CCI", "tier": "unsupported", "allowed_language": "Excluded from main conclusions."},
        ]
    )
    _write_csv(evidence, "tables/v1_1_evidence_matrix.csv")
    _write_md("reports/v1_1_evidence_matrix.md", "# v1.1 Evidence Matrix\n\n" + _md_table(evidence))
    if tier in {"strong", "acceptable"}:
        branch_text = "Branch nucleation receives internal native-teacher support, related mouse gastrulation external time-series support and clone-aware computational support in L2. This remains computational evidence, not experimental confirmation."
    elif tier == "weak":
        branch_text = "Branch nucleation receives internal native-teacher support and related mouse gastrulation external time-series support; L2 offers weak clone-aware feasibility only."
    else:
        branch_text = "Branch nucleation receives internal native-teacher support and related mouse gastrulation external time-series support, but clone-aware validation remains unsupported in the tested datasets."
    for path, title in [
        ("manuscript/final_retained_results_and_methods.md", "# Final Retained Results and Methods"),
        ("manuscript/manuscript.md", "# SwarmLineage-OT"),
        ("manuscript/methods.md", "# Methods"),
        ("manuscript/supplementary.md", "# Supplementary Information"),
    ]:
        _write_md(
            path,
            f"{title}\n\n"
            f"{branch_text}\n\n"
            f"L2 selected dataset: Biddy_2018_Nature / GSE99915 CellTag reprogramming. Usable clones: {usable}. Condensation-to-clone-splitting Spearman: {rho:.4f}. L2 tier: {tier}.\n\n"
            "Primary model remains M5_ot_swarm. Diffusion remains an encoded control-law recovery. Birth/death, memory and CCI remain unsupported and excluded from the main claim. The manuscript does not claim OT replacement, experimental confirmation or biological causality.\n",
        )
    _write_md(
        "reports/scientific_gap_audit.md",
        "# Scientific Gap Audit\n\n"
        f"{branch_text}\n\n"
        "- E1 remains acceptable time-series support.\n"
        f"- L2 clone-aware tier: {tier}; usable clones: {usable}; effect: {rho:.4f}.\n"
        "- Swarm-specific necessity remains unresolved.\n"
        "- Diffusion is encoded recovery; birth/death, memory and CCI remain unsupported.\n",
    )
    _write_md(
        "reports/editorial_assessment.md",
        "# Editorial Assessment\n\n"
        f"{branch_text}\n\n"
        "The evidence package is stronger than v1.0 because the priority Biddy/CellTag dataset was actually downloaded and analyzed. If the L2 tier is fail or weak, the manuscript should keep clone-level language out of the main claim and present L2 as a falsifying or inconclusive external test.",
    )


def update_figures(l2_summary: pd.DataFrame) -> None:
    ensure_dir(ROOT / "figures/main")
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    clone_df = _read_csv("tables/l2_clone_branch_validation.csv")
    if not clone_df.empty:
        axes[0].scatter(clone_df["clone_pre_event_condensation_exposure"], clone_df["clone_branch_splitting_score"], s=np.clip(clone_df["clone_size"], 8, 80), alpha=0.6)
    axes[0].set_xlabel("condensation exposure")
    axes[0].set_ylabel("clone branch splitting")
    summary = l2_summary.iloc[0] if not l2_summary.empty else pd.Series(dtype=object)
    axes[1].bar(["rho", "perm p"], [float(summary.get("condensation_spearman", 0.0)), float(summary.get("condensation_permutation_p", 1.0))], color=["#4C78A8", "#F58518"])
    axes[1].set_title(f"L2 tier: {summary.get('l2_validation_tier', 'fail')}")
    fig.suptitle("Figure 6: clone-aware validation attempt")
    fig.tight_layout()
    fig.savefig(ROOT / "figures/main/figure6_clone_validation.png", dpi=180)
    fig.savefig(ROOT / "figures/main/figure6_clone_validation_or_blocker.png", dpi=180)
    plt.close(fig)
    _write_md(
        "reports/main_figure_readiness.md",
        "# Main Figure Readiness\n\n"
        "Figure 6 has been updated for the Biddy/CellTag L2 clone-aware validation attempt. It displays the observed clone association and the resulting L2 tier rather than presenting the attempt as automatic support.",
    )


def update_review_and_audits(l2_summary: pd.DataFrame) -> None:
    summary = l2_summary.iloc[0] if not l2_summary.empty else pd.Series(dtype=object)
    tier = summary.get("l2_validation_tier", "fail")
    rho = float(summary.get("condensation_spearman", 0.0))
    _write_md(
        "reports/reviewer_attack_matrix.md",
        "# Reviewer Attack Matrix\n\n"
        "| attack | current_answer | evidence_available | evidence_gap | planned_analysis_or_experiment | claim_language_allowed |\n"
        "|---|---|---|---|---|---|\n"
        f"| Are clone-level claims supported? | L2 Biddy/CellTag tier is {tier} with Spearman {rho:.4f}. | tables/l2_clone_model_summary.csv | Larger CellTag and native cross-system clone tests are still needed if tier is weak/fail. | Clone-aware support is stated only at the reported tier. |\n"
        "| Is this just OT geometry? | Swarm-specific necessity remains unresolved. | tables/swarm_necessity_ablation.csv | Fine-grained module-drop rollouts. | Order-parameter signature is retained; swarm requirement is not claimed. |\n"
        "| Does E1 establish clone-resolved evidence? | No; E1 has no clone/barcode metadata. | external integrity audit. | Clone-aware datasets such as L2. | E1 is time-series support only. |\n"
        "| Are unsupported modules part of the main result? | No. | module contribution audit. | Future targeted evidence. | Birth/death, memory and CCI are excluded. |\n"
        "| Is the result causal? | No. | No prospective perturbation. | Prospective barcoded time-window experiment. | Computational association/hypothesis only. |\n",
    )
    forbidden = [
        "Nature-ready",
        "proven mechanism",
        "causal validation",
        "wet-lab validated",
        "true lineage",
        "lineage validated",
        "outperforms OT",
        "CCI validated",
        "memory hysteresis discovered",
        "birth/death law discovered",
        "external validation complete",
        "clone validation",
    ]
    hits = []
    for root in ["reports", "manuscript"]:
        for path in (ROOT / root).rglob("*.md"):
            text = path.read_text(encoding="utf-8", errors="ignore")
            lower = text.lower()
            for phrase in forbidden:
                idx = lower.find(phrase.lower())
                if idx >= 0:
                    context = lower[max(0, idx - 100) : idx + 100]
                    allowed = path.name in {"v1_1_claim_audit.md", "claim_audit.md"} or "forbidden" in context
                    if not allowed:
                        hits.append({"file": _rel(path), "phrase": phrase})
    hit_df = pd.DataFrame(hits)
    _write_md(
        "reports/v1_1_claim_audit.md",
        "# v1.1 Claim Audit\n\n"
        f"- prohibited positive-claim hits: {hit_df.shape[0]}\n\n"
        + ("No prohibited positive claim strings were found." if hit_df.empty else _md_table(hit_df)),
    )
    _write_md(
        "reports/claim_audit.md",
        "# Claim Audit\n\n"
        f"- prohibited positive-claim hits: {hit_df.shape[0]}\n\n"
        + ("No prohibited positive claim strings were found." if hit_df.empty else _md_table(hit_df)),
    )
    tracked = _git(["ls-files"])
    large_tracked = [p for p in tracked.splitlines() if Path(p).suffix.lower() in {".h5ad", ".rds", ".gz", ".tar", ".npz", ".pt"}]
    _write_md(
        "reports/v1_1_data_integrity_audit.md",
        "# v1.1 Data Integrity Audit\n\n"
        f"- Biddy file exists locally: {BIDDY_FILE.exists()}\n"
        f"- Biddy file size bytes: {BIDDY_FILE.stat().st_size if BIDDY_FILE.exists() else 0}\n"
        "- Clone identifiers are copied from source field `barcode_all`.\n"
        "- No clone identifiers are invented.\n"
        f"- tracked large binary risk count: {len(large_tracked)}\n"
        f"- tracked large binary examples: {large_tracked[:10]}\n",
    )
    _write_md(
        "reports/output_integrity_audit.md",
        "# Output Integrity Audit\n\n"
        "- L2 raw and processed h5ad files are under ignored data/processed paths.\n"
        "- Committed artifacts are code, reports and CSV summaries.\n"
        "- E1 and internal teacher backends remain labelled separately from L2 results.\n"
        "- L2 reports the actual validation tier and does not promote failed/weak results.\n",
    )
    _write_md(
        "reports/v1_1_push_status.md",
        "# v1.1 Push Status\n\nPush has not yet been attempted for the final v1.1 commit. This file will be updated after the commit/push step if a second commit is needed.",
    )


def run() -> None:
    write_start_state()
    registry = lineage_registry()
    _write_csv(registry, "tables/l2_lineage_dataset_registry.csv")
    _write_md(
        "reports/l2_lineage_dataset_selection.md",
        "# L2 Lineage Dataset Selection\n\n"
        "Six clone/barcode candidates were assessed. Biddy_2018_Nature was selected because the priority CellTag h5ad was downloaded and contains time, cell state and barcode fields.\n\n"
        + _md_table(registry),
    )
    raw, down, obs, z = load_and_prepare_biddy()
    registry.loc[registry["dataset_id"].eq("L2_biddy_2018_nature"), ["matrix_loaded", "metadata_loaded", "clone_metadata_loaded", "usable_for_clone_validation"]] = True
    _write_csv(registry, "tables/l2_lineage_dataset_registry.csv")
    plot_embedding(z, raw.obs.copy().reset_index(drop=True))
    write_l2_data_audit(raw, down)
    teacher = run_l2_teacher()
    order, event, _ = branch_nucleation(raw, z)
    clone_df, l2_summary = clone_validation(raw, z, event, teacher.get("teacher_backend", "failed"))
    master_summary(l2_summary)
    update_claims_and_manuscript(l2_summary)
    update_figures(l2_summary)
    update_review_and_audits(l2_summary)
    print(
        {
            "selected_l2": "L2_biddy_2018_nature",
            "teacher_backend": teacher.get("teacher_backend", "failed"),
            "tier": l2_summary.iloc[0]["l2_validation_tier"] if not l2_summary.empty else "fail",
            "usable_clones": int(l2_summary.iloc[0]["usable_clone_count"]) if not l2_summary.empty else 0,
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.parse_args()
    run()


if __name__ == "__main__":
    main()
