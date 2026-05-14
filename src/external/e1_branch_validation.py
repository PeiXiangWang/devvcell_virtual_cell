from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy import io, sparse
from sklearn.decomposition import TruncatedSVD

from src.discovery.branch_nucleation import run as run_branch_nucleation
from src.ot_teacher.build_teacher import build_teacher
from src.utils.config import ensure_dir, write_text


ROOT = Path(__file__).resolve().parents[2]
PRIMARY_COMPONENT = ROOT / "data/external/MouseGastrulationData/wt_chimera_sample1"
EXTERNAL_DATASET_ID = "E1_mouse_gastrulation_wt_chimera_sample1"
STAGE_MIN_CELLS = 80
MAX_CELLS_TOTAL = 1800


def _stage_to_float(value: object) -> float:
    text = str(value)
    match = pd.Series([text]).str.extract(r"([0-9]+\.?[0-9]*)", expand=False).iloc[0]
    try:
        return float(match)
    except Exception:
        return float("nan")


def _safe_read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def _write_yaml(path: Path, payload: dict) -> None:
    ensure_dir(path.parent)
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")


def _registry_row(
    *,
    dataset_id: str,
    dataset_name: str,
    source_type: str,
    accession: str,
    doi: str,
    url: str,
    publication: str,
    organism: str,
    system: str,
    time_or_stage_available: bool,
    cell_type_available: bool,
    lineage_or_barcode_available: bool,
    expression_matrix_available: bool,
    processed_h5ad_available: bool,
    raw_only: bool,
    download_attempted: bool,
    download_success: bool,
    matrix_loaded: bool,
    metadata_loaded: bool,
    usable_for_external_validation: bool,
    reason_if_not_usable: str,
    selected_as_external_E1: bool,
    notes: str,
) -> dict:
    return {
        "dataset_id": dataset_id,
        "dataset_name": dataset_name,
        "source_type": source_type,
        "accession": accession,
        "doi": doi,
        "url": url,
        "publication": publication,
        "organism": organism,
        "system": system,
        "time_or_stage_available": time_or_stage_available,
        "cell_type_available": cell_type_available,
        "lineage_or_barcode_available": lineage_or_barcode_available,
        "expression_matrix_available": expression_matrix_available,
        "processed_h5ad_available": processed_h5ad_available,
        "raw_only": raw_only,
        "download_attempted": download_attempted,
        "download_success": download_success,
        "matrix_loaded": matrix_loaded,
        "metadata_loaded": metadata_loaded,
        "usable_for_external_validation": usable_for_external_validation,
        "reason_if_not_usable": reason_if_not_usable,
        "selected_as_external_E1": selected_as_external_E1,
        "notes": notes,
    }


def build_registry(component_dir: Path = PRIMARY_COMPONENT) -> pd.DataFrame:
    obs = _safe_read_csv(component_dir / "obs.csv")
    var = _safe_read_csv(component_dir / "var.csv")
    matrix_exists = (component_dir / "matrix.mtx").exists()
    primary_usable = bool(matrix_exists and not obs.empty and not var.empty and {"stage.mapped", "celltype.mapped"}.issubset(obs.columns))
    rows = [
        _registry_row(
            dataset_id=EXTERNAL_DATASET_ID,
            dataset_name="Mouse gastrulation WT chimera sample 1",
            source_type="local_components_from_public_bioconductor_package",
            accession="MouseGastrulationData; raw WT chimera accessions E-MTAB-7324/E-MTAB-8812",
            doi="10.18129/B9.bioc.MouseGastrulationData; 10.1038/s41586-019-0933-9",
            url="https://bioconductor.org/packages/MouseGastrulationData/",
            publication="Pijuan-Sala et al., Nature 2019; MouseGastrulationData Bioconductor package",
            organism="Mus musculus",
            system="mouse gastrulation / early organogenesis chimera",
            time_or_stage_available="stage.mapped" in obs.columns,
            cell_type_available="celltype.mapped" in obs.columns,
            lineage_or_barcode_available=False,
            expression_matrix_available=matrix_exists,
            processed_h5ad_available=False,
            raw_only=False,
            download_attempted=False,
            download_success=primary_usable,
            matrix_loaded=False,
            metadata_loaded=not obs.empty,
            usable_for_external_validation=primary_usable,
            reason_if_not_usable="" if primary_usable else "Local matrix/metadata components missing or missing stage/cell-type columns.",
            selected_as_external_E1=primary_usable,
            notes="No network download was required in this run because matrix.mtx/obs.csv/var.csv components already exist locally.",
        ),
        _registry_row(
            dataset_id="E1_fallback_wot_reprogramming",
            dataset_name="Waddington-OT iPSC/fibroblast reprogramming time course",
            source_type="public_tutorial_remote_candidate",
            accession="WOT tutorial data links; Schiebinger et al. 2019",
            doi="10.1016/j.cell.2019.01.006",
            url="https://broadinstitute.github.io/wot/tutorial/",
            publication="Schiebinger et al., Cell 2019",
            organism="mouse",
            system="reprogramming time course",
            time_or_stage_available=True,
            cell_type_available=True,
            lineage_or_barcode_available=False,
            expression_matrix_available=True,
            processed_h5ad_available=False,
            raw_only=False,
            download_attempted=False,
            download_success=False,
            matrix_loaded=False,
            metadata_loaded=False,
            usable_for_external_validation=False,
            reason_if_not_usable="Fallback not downloaded because the primary local public MouseGastrulationData component was usable.",
            selected_as_external_E1=False,
            notes="Retained as fallback; no validation claim is made for this dataset.",
        ),
        _registry_row(
            dataset_id="E1_fallback_klein_esc_lif_withdrawal",
            dataset_name="mESC differentiation after LIF withdrawal",
            source_type="public_geo_remote_candidate",
            accession="GSE65525",
            doi="10.1016/j.cell.2015.04.044",
            url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65525",
            publication="Klein et al., Cell 2015",
            organism="mouse",
            system="mESC differentiation time course",
            time_or_stage_available=True,
            cell_type_available=False,
            lineage_or_barcode_available=False,
            expression_matrix_available=True,
            processed_h5ad_available=False,
            raw_only=False,
            download_attempted=False,
            download_success=False,
            matrix_loaded=False,
            metadata_loaded=False,
            usable_for_external_validation=False,
            reason_if_not_usable="Fallback not downloaded because the primary local public MouseGastrulationData component was usable; cell-type labels would need curation.",
            selected_as_external_E1=False,
            notes="Public fallback for unsupervised time-series support only.",
        ),
        _registry_row(
            dataset_id="E1_fallback_celltag_eb",
            dataset_name="Embryoid body differentiation with genetic recording",
            source_type="public_lineage_remote_candidate",
            accession="Publication-associated Cell Reports 2020 resources",
            doi="10.1016/j.celrep.2020.108222",
            url="https://www.sciencedirect.com/science/article/pii/S2211124720312110",
            publication="Kim et al., Cell Reports 2020",
            organism="mouse",
            system="embryoid body differentiation with genetic recording",
            time_or_stage_available=True,
            cell_type_available=True,
            lineage_or_barcode_available=True,
            expression_matrix_available=True,
            processed_h5ad_available=False,
            raw_only=False,
            download_attempted=False,
            download_success=False,
            matrix_loaded=False,
            metadata_loaded=False,
            usable_for_external_validation=False,
            reason_if_not_usable="Fallback not downloaded because the primary local public MouseGastrulationData component was usable; barcode processing would require a separate ingestion step.",
            selected_as_external_E1=False,
            notes="Potential lineage-aware follow-up; no lineage validation is claimed here.",
        ),
    ]
    out = pd.DataFrame(rows)
    ensure_dir("tables")
    out.to_csv(ROOT / "tables/external_dataset_registry.csv", index=False)
    return out


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
    return np.asarray(sorted(selected), dtype=int)


def _load_component_anndata(component_dir: Path) -> ad.AnnData:
    obs = pd.read_csv(component_dir / "obs.csv")
    var = pd.read_csv(component_dir / "var.csv")
    matrix = io.mmread(component_dir / "matrix.mtx").tocsr()
    if matrix.shape == (var.shape[0], obs.shape[0]):
        matrix = matrix.T.tocsr()
    elif matrix.shape != (obs.shape[0], var.shape[0]):
        raise ValueError(f"Matrix shape {matrix.shape} does not match obs={obs.shape[0]} and var={var.shape[0]}")
    var_names = var.get("SYMBOL", var.iloc[:, 0]).astype(str).replace({"nan": ""})
    var_names = np.where(var_names == "", var.get("ENSEMBL", var.iloc[:, 0]).astype(str), var_names)
    adata = ad.AnnData(matrix, obs=obs.copy(), var=var.copy())
    adata.var_names = pd.Index(var_names).astype(str)
    adata.obs_names = obs.get("cell", pd.Series([f"cell_{i}" for i in range(obs.shape[0])])).astype(str).to_numpy()
    adata.var_names_make_unique()
    return adata


def _preprocess_external(adata: ad.AnnData, seed: int = 17) -> ad.AnnData:
    obs = adata.obs.copy()
    obs["time_numeric"] = obs["stage.mapped"].map(_stage_to_float)
    obs["time_point"] = obs["stage.mapped"].astype(str)
    obs["lineage"] = obs["celltype.mapped"].astype(str).replace({"nan": "unknown", "": "unknown"})
    obs["cell_type"] = obs["lineage"]
    obs["external_dataset_id"] = EXTERNAL_DATASET_ID
    obs["external_source"] = "MouseGastrulationData/WTChimeraData sample1"
    obs["split_role"] = "train"
    obs["lineage_barcode"] = ""
    adata.obs = obs
    valid = np.isfinite(obs["time_numeric"].to_numpy(dtype=float))
    valid &= obs["lineage"].astype(str).to_numpy() != "Doublet"
    stage_counts = obs.loc[valid, "time_point"].value_counts()
    keep_stages = sorted(stage_counts[stage_counts >= STAGE_MIN_CELLS].index, key=_stage_to_float)
    valid &= obs["time_point"].isin(keep_stages).to_numpy()
    adata = adata[valid].copy()
    sample_idx = _stratified_sample(adata.obs, MAX_CELLS_TOTAL, seed)
    adata = adata[sample_idx].copy()
    x = adata.X.tocsr() if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
    totals = np.asarray(x.sum(axis=1)).ravel()
    n_genes = np.asarray((x > 0).sum(axis=1)).ravel()
    keep = (totals > 0) & (n_genes >= 50)
    adata = adata[keep].copy()
    x = adata.X.tocsr() if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
    totals = np.asarray(x.sum(axis=1)).ravel()
    adata.obs["n_counts"] = totals
    adata.obs["n_genes"] = np.asarray((x > 0).sum(axis=1)).ravel()
    scale = np.divide(10000.0, np.maximum(totals, 1), out=np.ones_like(totals, dtype=float), where=totals > 0)
    x = x.multiply(scale[:, None]).tocsr()
    x.data = np.log1p(x.data)
    adata.X = x
    mean = np.asarray(x.mean(axis=0)).ravel()
    mean_sq = np.asarray(x.power(2).mean(axis=0)).ravel()
    var = np.maximum(mean_sq - mean**2, 0.0)
    top = np.argsort(var)[::-1][: min(2000, adata.n_vars)]
    mask = np.zeros(adata.n_vars, dtype=bool)
    mask[top] = True
    adata.var["highly_variable"] = mask
    adata = adata[:, mask].copy()
    n_pcs = min(30, adata.n_obs - 1, adata.n_vars - 1)
    svd = TruncatedSVD(n_components=n_pcs, random_state=seed)
    adata.obsm["X_pca"] = svd.fit_transform(adata.X).astype(np.float32)
    adata.obsm["X_umap"] = adata.obsm["X_pca"][:, :2].copy()
    adata.obs["cell_cycle_score"] = 0.0
    adata.uns["external_e1"] = {
        "dataset_id": EXTERNAL_DATASET_ID,
        "source_component": str(PRIMARY_COMPONENT),
        "kept_stages": keep_stages,
        "stage_min_cells": STAGE_MIN_CELLS,
        "max_cells_total": MAX_CELLS_TOTAL,
        "lineage_barcode_available": False,
    }
    return adata


def _plot_external_overview(adata: ad.AnnData) -> None:
    fig_dir = ensure_dir(ROOT / "figures/external")
    emb = adata.obsm["X_umap"]
    time_codes = pd.Categorical(adata.obs["time_point"].astype(str), categories=sorted(adata.obs["time_point"].astype(str).unique(), key=_stage_to_float), ordered=True).codes
    fig, ax = plt.subplots(figsize=(6, 5), dpi=160)
    sc = ax.scatter(emb[:, 0], emb[:, 1], c=time_codes, s=6, cmap="viridis", linewidths=0, alpha=0.75)
    ax.set_title("E1 external MouseGastrulationData by stage")
    ax.set_xlabel("PCA 1")
    ax.set_ylabel("PCA 2")
    fig.colorbar(sc, ax=ax, fraction=0.035, label="stage code")
    fig.tight_layout()
    fig.savefig(fig_dir / "e1_external_umap_by_time.png")
    plt.close(fig)

    top_lineages = adata.obs["lineage"].value_counts().head(12).index.tolist()
    labels = pd.Categorical(np.where(adata.obs["lineage"].isin(top_lineages), adata.obs["lineage"], "other"))
    fig, ax = plt.subplots(figsize=(7, 5), dpi=160)
    sc = ax.scatter(emb[:, 0], emb[:, 1], c=labels.codes, s=6, cmap="tab20", linewidths=0, alpha=0.75)
    ax.set_title("E1 external MouseGastrulationData by lineage")
    ax.set_xlabel("PCA 1")
    ax.set_ylabel("PCA 2")
    fig.colorbar(sc, ax=ax, fraction=0.035, label="lineage code")
    fig.tight_layout()
    fig.savefig(fig_dir / "e1_external_umap_by_lineage.png")
    plt.close(fig)


def build_external_anndata(registry: pd.DataFrame) -> ad.AnnData | None:
    selected = registry[registry["selected_as_external_E1"].astype(bool)]
    if selected.empty:
        return None
    adata_raw = _load_component_anndata(PRIMARY_COMPONENT)
    adata = _preprocess_external(adata_raw)
    ensure_dir(ROOT / "data/external")
    ensure_dir(ROOT / "processed/external")
    adata.write_h5ad(ROOT / "data/external/e1_external_input.h5ad")
    adata.write_h5ad(ROOT / "processed/external/e1_swarmlineage_input.h5ad")
    registry = registry.copy()
    primary_mask = registry["dataset_id"] == EXTERNAL_DATASET_ID
    registry.loc[primary_mask, ["matrix_loaded", "metadata_loaded", "processed_h5ad_available", "download_success", "usable_for_external_validation"]] = True
    registry.loc[primary_mask, "notes"] = (
        "Local public MouseGastrulationData components were loaded and converted to "
        "data/external/e1_external_input.h5ad and processed/external/e1_swarmlineage_input.h5ad."
    )
    registry.to_csv(ROOT / "tables/external_dataset_registry.csv", index=False)
    counts = adata.obs.groupby(["time_point", "time_numeric", "lineage"], observed=False).size().reset_index(name="n")
    counts.to_csv(ROOT / "tables/external_time_lineage_counts.csv", index=False)
    _plot_external_overview(adata)
    write_text(
        ROOT / "reports/external_data_audit.md",
        "\n".join(
            [
                "# External Data Audit",
                "",
                f"- selected_dataset_id: {EXTERNAL_DATASET_ID}",
                "- source: local matrix.mtx/obs.csv/var.csv components derived from public MouseGastrulationData.",
                f"- cells_after_filtering: {adata.n_obs}",
                f"- genes_after_hvg: {adata.n_vars}",
                f"- time/stage field: `stage.mapped` -> `time_numeric` / `time_point`",
                "- lineage field: `celltype.mapped` -> `lineage` / `cell_type`",
                "- lineage barcode: unavailable; no lineage-validated claim is made.",
                "- embedding fit: external data only; no internal teacher or internal embedding was used.",
                "",
                "## Stage Counts",
                "",
                adata.obs["time_point"].value_counts().sort_index().to_markdown(),
                "",
            ]
        ),
    )
    write_text(
        ROOT / "reports/external_leakage_audit.md",
        "\n".join(
            [
                "# External Leakage Audit",
                "",
                "- External preprocessing was fit only on the selected external MouseGastrulationData component.",
                "- Internal SwarmLineage-OT data, internal moscot teacher and internal learned embeddings were not used to fit the external PCA.",
                "- `split_role` is set to train for all external cells because E1 tests branch-nucleation order parameters, not held-out prediction.",
                "- No lineage barcodes were present in this external component; lineage validation is not claimed.",
                "",
            ]
        ),
    )
    return adata


def write_configs() -> dict[str, Path]:
    ot_cfg = {
        "adata_path": "processed/external/e1_swarmlineage_input.h5ad",
        "teacher_path": "processed/external/e1_ot_teacher.h5ad",
        "couplings_dir": "processed/external/ot_couplings",
        "fate_probabilities_path": "processed/external/e1_ot_fate_probabilities.parquet",
        "time_key": "time_numeric",
        "time_label_key": "time_point",
        "cell_type_key": "lineage",
        "latent_key": "X_pca",
        "epsilon": 0.08,
        "max_cells_per_time": 120,
        "max_terminal_fates": 8,
        "native_moscot_timeout_seconds": 90,
        "use_native_moscot": True,
        "native_max_cells_per_time": 120,
        "native_max_iterations": 350,
        "native_jit": False,
        "native_device": "cpu",
        "random_seed": 17,
        "split_mode": "none",
        "holdout_time": 7.75,
        "teacher_backend": "native_moscot",
        "teacher_index_path": "processed/external/ot_couplings/teacher_coupling_index.csv",
        "figure_dir": "figures/external",
        "report_dir": "reports/external",
        "table_dir": "tables/external",
        "summary_path": "processed/external/e1_ot_teacher_summary.json",
    }
    model_cfg = {
        "teacher_path": "processed/external/e1_ot_teacher.h5ad",
        "model_dir": "results/external/e1/models",
        "simulation_dir": "results/external/e1/simulations",
        "metrics_path": "tables/external/e1_final_metrics.csv",
        "holdout_time": 7.75,
        "split_mode": "none",
        "time_key": "time_numeric",
        "cell_type_key": "lineage",
        "latent_key": "X_pca",
        "teacher_velocity_key": "X_ot_velocity",
        "seeds": [7, 17, 23, 42, 99],
        "epochs": 18,
        "quick_epochs": 8,
        "rollout_steps": 4,
        "dt": 0.25,
        "batch_size": 384,
        "learning_rate": 0.001,
        "hidden_dim": 96,
        "simulation_cells_per_seed": 300,
        "metric_sample_size": 300,
        "event_log_path": "tables/external/e1_birth_death_event_log.csv",
        "order_log_path": "tables/external/e1_rollout_order_parameters.csv",
        "baseline_execution_matrix_path": "reports/external/e1_baseline_execution_matrix.csv",
    }
    discovery_cfg = {
        "teacher_fidelity": {
            "strong": {"relative_sinkhorn_max": 1.10, "relative_mmd_max": 1.50, "composition_rmse_max": 0.005, "manifold_escape_rate_max": 0.08},
            "acceptable": {"relative_sinkhorn_max": 1.25, "relative_mmd_max": 2.50, "composition_rmse_max": 0.025, "manifold_escape_rate_max": 0.20},
            "weak": {"relative_sinkhorn_max": 1.50, "relative_mmd_max": 4.00, "composition_rmse_max": 0.05, "manifold_escape_rate_max": 0.35},
            "fail_above_weak": True,
        },
        "emergent_law": {
            "min_seed_count": 5,
            "quick_min_seed_count": 2,
            "min_effect_size_for_acceptable": 0.01,
            "min_effect_size_for_strong": 0.03,
            "max_permutation_q_for_acceptable": 0.10,
            "max_permutation_q_for_strong": 0.05,
            "require_negative_control_for_strong": True,
            "require_rollout_based_evidence_for_strong": True,
            "permutation_repeats": 100,
            "bootstrap_repeats": 500,
        },
        "output": {"tables_dir": "tables/external", "reports_dir": "reports/external", "figures_dir": "figures/external"},
    }
    train_cfg = {
        "model_config": "configs/external_e1_model.yaml",
        "discovery_config": "configs/external_e1_discovery.yaml",
        "ablation_metrics_path": "tables/external/e1_ablation_metrics.csv",
        "ablation_stats_path": "tables/external/e1_ablation_statistical_tests.csv",
        "negative_results_report": "reports/external/e1_negative_results.md",
        "ablation_report": "reports/external/e1_ablation_interpretation.md",
        "random_seed": 17,
        "baseline_execution_matrix_path": "reports/external/e1_baseline_execution_matrix.csv",
        "module_contribution_report": "reports/external/e1_module_contribution_audit.md",
        "scientific_gap_report": "reports/external/e1_scientific_gap_audit.md",
    }
    paths = {
        "ot": ROOT / "configs/external_e1_ot_teacher.yaml",
        "model": ROOT / "configs/external_e1_model.yaml",
        "discovery": ROOT / "configs/external_e1_discovery.yaml",
        "train": ROOT / "configs/external_e1_train.yaml",
    }
    _write_yaml(paths["ot"], ot_cfg)
    _write_yaml(paths["model"], model_cfg)
    _write_yaml(paths["discovery"], discovery_cfg)
    _write_yaml(paths["train"], train_cfg)
    return paths


def _native_python() -> Path:
    candidate = ROOT / ".venv_moscot_native/Scripts/python.exe"
    return candidate if candidate.exists() else Path(sys.executable)


def _run_command(cmd: list[str], timeout: int, native_env: bool = False) -> dict:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT)
    if native_env:
        env["JAX_PLATFORMS"] = "cpu"
        env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    start = time.perf_counter()
    try:
        out = subprocess.check_output(cmd, cwd=ROOT, text=True, stderr=subprocess.STDOUT, timeout=timeout, env=env)
        return {"success": True, "runtime_seconds": time.perf_counter() - start, "output": out[-4000:], "error": ""}
    except subprocess.TimeoutExpired as exc:
        return {"success": False, "runtime_seconds": time.perf_counter() - start, "output": (exc.output or "")[-4000:] if isinstance(exc.output, str) else "", "error": f"timeout after {timeout}s"}
    except subprocess.CalledProcessError as exc:
        return {"success": False, "runtime_seconds": time.perf_counter() - start, "output": (exc.output or "")[-4000:], "error": f"exit {exc.returncode}"}


def run_external_teacher(paths: dict[str, Path]) -> tuple[str, pd.DataFrame]:
    ensure_dir(ROOT / "processed/external/ot_couplings")
    cmd = [str(_native_python()), "-m", "src.ot_teacher.run_moscot", "--config", str(paths["ot"]), "--try-native"]
    result = _run_command(cmd, timeout=900, native_env=True)
    if not result["success"]:
        fallback_cmd = [sys.executable, "-m", "src.ot_teacher.run_moscot", "--config", str(paths["ot"])]
        fallback = _run_command(fallback_cmd, timeout=300, native_env=False)
        result = {**fallback, "native_failure": result}
    teacher_result = build_teacher(yaml.safe_load(paths["ot"].read_text(encoding="utf-8")))
    index = pd.read_csv(ROOT / "processed/external/ot_couplings/teacher_coupling_index.csv")
    backend = str(teacher_result.get("teacher_backend", index.get("teacher_backend", pd.Series(["failed"])).iloc[0] if not index.empty else "failed"))
    pairs = index.copy()
    pairs["external_teacher_backend"] = backend
    pairs["run_success"] = bool(result["success"])
    pairs["run_runtime_seconds"] = float(result["runtime_seconds"])
    pairs["run_error"] = result.get("error", "")
    ensure_dir(ROOT / "tables")
    pairs.to_csv(ROOT / "tables/external_native_teacher_pairs.csv", index=False)
    shutil.copyfile(ROOT / "reports/external/ot_teacher_report.md", ROOT / "reports/external_native_teacher_report.md")
    return backend, pairs


def run_external_training(paths: dict[str, Path]) -> dict:
    cmd = [sys.executable, "-m", "src.train.train_model", "--config", str(paths["model"])]
    return _run_command(cmd, timeout=1200, native_env=False)


def _copy_external_branch_outputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mapping = {
        "tables/external/swarm_order_parameters.csv": "tables/external_branch_nucleation_order_parameters.csv",
        "tables/external/branch_nucleation_event_windows.csv": "tables/external_branch_nucleation_event_windows.csv",
        "tables/external/branch_nucleation_negative_controls.csv": "tables/external_branch_nucleation_negative_controls.csv",
    }
    for src, dst in mapping.items():
        src_path = ROOT / src
        if src_path.exists():
            shutil.copyfile(src_path, ROOT / dst)
    return (
        _safe_read_csv(ROOT / "tables/external/branch_nucleation_model_comparison.csv"),
        _safe_read_csv(ROOT / "tables/external/branch_nucleation_event_windows.csv"),
        _safe_read_csv(ROOT / "tables/external/branch_nucleation_negative_controls.csv"),
    )


def _external_validation_tier(backend: str, model: pd.DataFrame, windows: pd.DataFrame, controls: pd.DataFrame) -> tuple[str, str, dict]:
    selected = model[model["variant"] == "M5_ot_swarm"] if not model.empty and "variant" in model else pd.DataFrame()
    if selected.empty or windows.empty:
        return "fail", "unsupported", {"reason": "No M5 branch event window detected."}
    row = selected.iloc[0].to_dict()
    effect = float(row.get("lineage_separation_effect", np.nan))
    alignment = float(row.get("local_velocity_alignment_A_effect", np.nan))
    entropy = float(row.get("fate_entropy_H_effect", np.nan))
    direction_match = bool(np.isfinite(effect) and effect < 0)
    branch_event = bool(np.isfinite(effect) and abs(effect) >= 0.01)
    shuffled = controls[controls["control"].astype(str).str.contains("shuffled", na=False)] if not controls.empty and "control" in controls else pd.DataFrame()
    control_pass = bool((shuffled.get("gate_tier", pd.Series(dtype=str)).astype(str) == "fail").any()) if not shuffled.empty else False
    seed_stable = bool(row.get("seed_stability_pass", False))
    if backend != "native_moscot":
        tier = "weak" if branch_event and direction_match else "fail"
    elif direction_match and branch_event and control_pass and seed_stable:
        tier = "acceptable"
    elif direction_match and branch_event:
        tier = "weak"
    else:
        tier = "fail"
    interpretation = "transient_condensation_before_divergence" if direction_match and branch_event else "unsupported"
    return tier, interpretation, {
        "effect": effect,
        "alignment_effect": alignment,
        "entropy_effect": entropy,
        "direction_match": direction_match,
        "branch_event_detected": branch_event,
        "negative_control_pass": control_pass,
        "seed_stability_pass": seed_stable,
    }


def _plot_branch_validation(model: pd.DataFrame, windows: pd.DataFrame, controls: pd.DataFrame, tier_summary: pd.DataFrame) -> None:
    fig_dir = ensure_dir(ROOT / "figures/external")
    order = _safe_read_csv(ROOT / "tables/external/e1_rollout_order_parameters.csv")
    fig, ax = plt.subplots(figsize=(7, 4), dpi=160)
    if not order.empty:
        subset = order[order["variant"] == "M5_ot_swarm"] if "variant" in order else order
        plot = subset.groupby("step", observed=False).mean(numeric_only=True).reset_index()
        for col in ["local_velocity_alignment_A", "lineage_separation_S", "fate_entropy_H", "branch_imbalance_B", "local_density_mean"]:
            vals = plot[col].to_numpy(dtype=float)
            vals = (vals - np.nanmin(vals)) / max(np.nanmax(vals) - np.nanmin(vals), 1e-8)
            ax.plot(plot["step"], vals, marker="o", label=col)
    ax.set_title("E1 branch-nucleation order parameters")
    ax.set_xlabel("rollout step")
    ax.set_ylabel("normalized value")
    ax.legend(fontsize=6, frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "e1_branch_nucleation_order_parameters.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4), dpi=160)
    if not windows.empty:
        effect_cols = [c for c in windows.columns if c.endswith("_effect")]
        means = windows[effect_cols].mean().sort_values()
        ax.barh(range(len(means)), means, color="#7da7c7")
        ax.set_yticks(range(len(means)), [c.replace("_effect", "") for c in means.index], fontsize=7)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_title("E1 branch event window effects")
    ax.set_xlabel("post - pre")
    fig.tight_layout()
    fig.savefig(fig_dir / "e1_branch_event_window.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4), dpi=160)
    if not controls.empty:
        vals = controls["effect_size"].astype(float)
        ax.barh(controls["control"].astype(str), vals, color="#c9a66b")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_title("E1 negative/control effects")
    ax.set_xlabel("effect size")
    fig.tight_layout()
    fig.savefig(fig_dir / "e1_branch_nucleation_controls.png")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(9, 7), dpi=160)
    counts = _safe_read_csv(ROOT / "tables/external_time_lineage_counts.csv")
    if not counts.empty:
        stage_counts = counts.groupby("time_point", observed=False)["n"].sum()
        stage_counts.plot(kind="bar", ax=axes[0, 0], color="#7da7c7", title="External stages")
    if not model.empty:
        axes[0, 1].barh(model["variant"].astype(str), model["lineage_separation_effect"].astype(float), color="#7da7c7")
        axes[0, 1].axvline(0, color="black", lw=0.8)
        axes[0, 1].set_title("External model effects")
    if not controls.empty:
        axes[1, 0].barh(controls["control"].astype(str), controls["effect_size"].astype(float), color="#c9a66b")
        axes[1, 0].axvline(0, color="black", lw=0.8)
        axes[1, 0].set_title("Controls")
    if not tier_summary.empty:
        axes[1, 1].axis("off")
        text = "\n".join([f"{k}: {tier_summary.iloc[0][k]}" for k in ["external_validation_tier", "external_teacher_backend", "condensation_before_divergence_reproduced"]])
        axes[1, 1].text(0.05, 0.7, text, va="top")
    fig.tight_layout()
    fig.savefig(ROOT / "figures/main/figure5_external_validation.png")
    plt.close(fig)


def summarize_external(backend: str, train_result: dict) -> dict:
    model, windows, controls = _copy_external_branch_outputs()
    tier, interpretation, evidence = _external_validation_tier(backend, model, windows, controls)
    selected = model[model["variant"] == "M5_ot_swarm"].iloc[0].to_dict() if not model.empty and (model["variant"] == "M5_ot_swarm").any() else {}
    summary = pd.DataFrame(
        [
            {
                "dataset_id": EXTERNAL_DATASET_ID,
                "external_validation_tier": tier,
                "external_teacher_backend": backend,
                "branch_event_detected": evidence.get("branch_event_detected", False),
                "condensation_before_divergence_reproduced": interpretation == "transient_condensation_before_divergence",
                "lineage_validated": False,
                "lineage_separation_effect": evidence.get("effect", np.nan),
                "alignment_effect": evidence.get("alignment_effect", np.nan),
                "entropy_effect": evidence.get("entropy_effect", np.nan),
                "negative_control_pass": evidence.get("negative_control_pass", False),
                "seed_stability_pass": evidence.get("seed_stability_pass", False),
                "interpretation": interpretation,
                "status": "analyzed" if tier != "fail" else "analyzed_but_not_supported",
                "blocker": "" if tier != "fail" else evidence.get("reason", "External branch signature was not supported."),
            }
        ]
    )
    summary.to_csv(ROOT / "tables/external_validation_tier_summary.csv", index=False)
    summary.to_csv(ROOT / "tables/external_branch_nucleation_summary.csv", index=False)

    internal = _safe_read_csv(ROOT / "tables/branch_nucleation_model_comparison.csv")
    internal_m5 = internal[internal["variant"] == "M5_ot_swarm"].iloc[0].to_dict() if not internal.empty and (internal["variant"] == "M5_ot_swarm").any() else {}
    comp = pd.DataFrame(
        [
            {
                "internal_effect_size": internal_m5.get("lineage_separation_effect", np.nan),
                "external_effect_size": selected.get("lineage_separation_effect", np.nan),
                "effect_direction_match": bool(np.sign(float(internal_m5.get("lineage_separation_effect", 0))) == np.sign(float(selected.get("lineage_separation_effect", 1)))) if selected else False,
                "alignment_effect_match": bool(np.sign(float(internal_m5.get("local_velocity_alignment_A_effect", 0))) == np.sign(float(selected.get("local_velocity_alignment_A_effect", 1)))) if selected else False,
                "separation_effect_match": bool(np.sign(float(internal_m5.get("lineage_separation_effect", 0))) == np.sign(float(selected.get("lineage_separation_effect", 1)))) if selected else False,
                "entropy_effect_match": bool(np.sign(float(internal_m5.get("fate_entropy_H_effect", 0))) == np.sign(float(selected.get("fate_entropy_H_effect", 1)))) if selected else False,
                "density_effect_match": bool(np.sign(float(internal_m5.get("local_density_mean_effect", 0))) == np.sign(float(selected.get("local_density_mean_effect", 1)))) if selected else False,
                "negative_control_status": "at_least_one_shuffled_control_failed_to_reproduce" if evidence.get("negative_control_pass", False) else "controls_insufficient",
                "teacher_backend_internal": "native_moscot",
                "teacher_backend_external": backend,
                "interpretation": interpretation if selected else "external_missing_or_failed",
            }
        ]
    )
    comp.to_csv(ROOT / "tables/internal_external_branch_nucleation_comparison.csv", index=False)

    fig, ax = plt.subplots(figsize=(6, 4), dpi=160)
    vals = [comp["internal_effect_size"].iloc[0], comp["external_effect_size"].iloc[0]]
    ax.bar(["internal", "external"], vals, color=["#7da7c7", "#c9a66b"])
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel("lineage separation event effect")
    ax.set_title("Internal vs external branch signature")
    fig.tight_layout()
    ensure_dir(ROOT / "figures/external")
    fig.savefig(ROOT / "figures/external/internal_external_branch_signature_comparison.png")
    plt.close(fig)

    _plot_branch_validation(model, windows, controls, summary)
    write_external_reports(summary, comp, model, windows, controls, train_result)
    return {"tier": tier, "interpretation": interpretation, "backend": backend}


def write_external_reports(summary: pd.DataFrame, comp: pd.DataFrame, model: pd.DataFrame, windows: pd.DataFrame, controls: pd.DataFrame, train_result: dict) -> None:
    tier = summary["external_validation_tier"].iloc[0]
    if tier in {"acceptable", "strong"}:
        conclusion = "External time-series support was observed for the branch-nucleation order-parameter signature. This remains computational evidence, not experimental validation."
    elif tier == "weak":
        conclusion = "External validation was initiated and provides partial feasibility support, but not decisive validation."
    else:
        conclusion = "External validation was attempted but remains unresolved or unsupported."
    write_text(
        ROOT / "reports/external_branch_nucleation_validation.md",
        "\n".join(
            [
                "# External Branch-Nucleation Validation",
                "",
                conclusion,
                "",
                "## Tier Summary",
                "",
                summary.to_markdown(index=False),
                "",
                "## Model Comparison",
                "",
                model.to_markdown(index=False) if not model.empty else "No external model comparison was produced.",
                "",
                "## Event Windows",
                "",
                windows.head(20).to_markdown(index=False) if not windows.empty else "No external event windows were detected.",
                "",
                "## Negative Controls",
                "",
                controls.to_markdown(index=False) if not controls.empty else "No external controls were produced.",
                "",
                f"Training command status: success={train_result.get('success')}, error={train_result.get('error')}",
                "",
                "No lineage barcode was available in this selected external component, so this is external time-series support only, not lineage validation.",
            ]
        ),
    )
    write_text(
        ROOT / "reports/external_validation_summary.md",
        "\n".join(
            [
                "# External Validation Summary",
                "",
                conclusion,
                "",
                summary.to_markdown(index=False),
                "",
                "Unsupported modules remain excluded from the main claim: birth/death, memory and CCI.",
                "",
            ]
        ),
    )
    write_text(
        ROOT / "reports/internal_external_branch_comparison.md",
        "\n".join(
            [
                "# Internal vs External Branch-Nucleation Comparison",
                "",
                comp.to_markdown(index=False),
                "",
                "Interpretation is restricted to computational order parameters. Causality and experimental mechanism are not established.",
                "",
            ]
        ),
    )
    write_text(
        ROOT / "reports/external_data_integrity_audit.md",
        "\n".join(
            [
                "# External Data Integrity Audit",
                "",
                "- Registry is not described as validation unless `external_validation_tier` is acceptable or strong.",
                "- External teacher backend is read from the actual generated teacher summary and coupling index.",
                "- Time-series support is not lineage validation because no lineage barcode was present in the selected component.",
                "- Computational branch signature is not described as a biological mechanism.",
                "- Failed or unattempted fallback downloads remain marked as not usable for this run.",
                "- Birth/death, memory and CCI remain unsupported and are not restored to the main claim.",
                "- Primary model remains evidence-selected M5, not M9 by default.",
                "",
                summary.to_markdown(index=False),
            ]
        ),
    )


def update_main_documents(result: dict) -> None:
    tier = result["tier"]
    if tier in {"acceptable", "strong"}:
        external_text = (
            "Branch nucleation, interpreted as transient condensation-before-divergence, is supported internally under native moscot teacher sensitivity "
            "and receives external time-series support in MouseGastrulationData WT chimera sample 1. This remains computational evidence, not experimental validation."
        )
    elif tier == "weak":
        external_text = "External validation was initiated. The external dataset provides partial feasibility support but not decisive validation."
    else:
        external_text = "External validation was attempted but remains unresolved. Internal branch nucleation remains a native-teacher computational hypothesis only."
    for path in [
        ROOT / "manuscript/final_retained_results_and_methods.md",
        ROOT / "manuscript/manuscript.md",
        ROOT / "manuscript/methods.md",
        ROOT / "reports/scientific_gap_audit.md",
        ROOT / "reports/editorial_assessment.md",
        ROOT / "reports/branch_nucleation_mechanism_summary.md",
    ]:
        existing = path.read_text(encoding="utf-8", errors="ignore") if path.exists() else f"# {path.stem}\n"
        marker = "\n\n## External Experiment E1\n\n"
        base = existing.split(marker)[0].rstrip()
        write_text(
            path,
            base
            + marker
            + external_text
            + "\n\n"
            + "- selected external dataset: MouseGastrulationData WT chimera sample 1.\n"
            + f"- external teacher backend: {result['backend']}.\n"
            + f"- external validation tier: {tier}.\n"
            + "- no experimental lineage tracing or experimental validation is claimed; causality and high-impact readiness are not established.\n"
            + "- diffusion remains encoded control-law recovery; birth/death, memory and CCI remain unsupported.\n",
        )


def claim_audit() -> None:
    forbidden = [
        "SwarmLineage-OT beats OT",
        "outperforms OT",
        "Nature-ready",
        "wet-lab validated",
        "causal validation",
        "true lineage",
        "CCI validated",
        "memory hysteresis discovered",
        "birth/death law discovered",
    ]
    rows = []
    for root in ["manuscript", "reports"]:
        for path in (ROOT / root).glob("*.md"):
            text = path.read_text(encoding="utf-8", errors="ignore")
            hits = [term for term in forbidden if term.lower() in text.lower()]
            if hits and path.name not in {"claim_audit.md", "external_data_integrity_audit.md"}:
                rows.append({"file": str(path.relative_to(ROOT)), "hits": "; ".join(hits), "status": "review_context_required"})
    audit = pd.DataFrame(rows)
    write_text(
        ROOT / "reports/claim_audit.md",
        "\n".join(
            [
                "# Claim Audit",
                "",
                "Prohibited positive claims were searched across manuscript and report markdown. Hits, if present, require negated limitation context.",
                "",
                audit.to_markdown(index=False) if not audit.empty else "No prohibited positive claim strings were found outside audit context.",
                "",
            ]
        ),
    )
    write_text(
        ROOT / "reports/output_integrity_audit.md",
        "\n".join(
            [
                "# Output Integrity Audit",
                "",
                "- Main internal outputs remain under top-level reports/tables.",
                "- External E1 outputs are under `processed/external`, `tables/external*`, `reports/external*`, and `figures/external`.",
                "- External generated h5ad/couplings are not committed because processed data are gitignored.",
                "- Registry-only fallback candidates are not described as validation results.",
                "- Native and fallback teacher backends are reported explicitly.",
                "",
            ]
        ),
    )


def run() -> dict:
    ensure_dir(ROOT / "reports")
    ensure_dir(ROOT / "tables")
    ensure_dir(ROOT / "figures/external")
    registry = build_registry()
    write_text(
        ROOT / "reports/external_dataset_selection.md",
        "\n".join(
            [
                "# External Dataset Selection",
                "",
                "E1-primary selected MouseGastrulationData WT chimera sample 1 because local matrix, metadata, stage and cell-type fields were present.",
                "Fallback candidates were registered but not downloaded because the primary public component was usable.",
                "",
                registry.to_markdown(index=False),
                "",
            ]
        ),
    )
    adata = build_external_anndata(registry)
    if adata is None:
        fail = {"tier": "fail", "interpretation": "unsupported", "backend": "failed"}
        update_main_documents(fail)
        claim_audit()
        return fail
    paths = write_configs()
    backend, _pairs = run_external_teacher(paths)
    train_result = run_external_training(paths)
    run_branch_nucleation(str(paths["train"]), quick_fixture=False)
    result = summarize_external(backend, train_result)
    update_main_documents(result)
    claim_audit()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.parse_args()
    print(json.dumps(run(), indent=2))


if __name__ == "__main__":
    main()
