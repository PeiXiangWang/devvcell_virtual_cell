from __future__ import annotations

import argparse
import logging
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD

from src.data.load import audit_data_files, discover_data_files, first_existing_column, load_h5ad
from src.data.fixtures import write_synthetic_fixture
from src.data.schema import CELL_CYCLE_MARKERS, gene_symbols
from src.utils.config import ensure_dir, load_config, write_json, write_text


def _setup_log(path: str) -> None:
    ensure_dir(Path(path).parent)
    logging.basicConfig(filename=path, level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", force=True)


def _stratified_sample(obs: pd.DataFrame, max_cells: int, seed: int, keys: list[str]) -> np.ndarray:
    if obs.shape[0] <= max_cells:
        return np.arange(obs.shape[0])
    rng = np.random.default_rng(seed)
    group_key = obs[keys].astype(str).agg("|".join, axis=1) if keys else pd.Series("all", index=obs.index)
    selected: list[int] = []
    per_group = max(15, max_cells // max(1, group_key.nunique()))
    for _, idx in group_key.groupby(group_key).groups.items():
        arr = np.asarray(list(idx))
        positions = obs.index.get_indexer(arr)
        take = min(len(positions), per_group)
        selected.extend(rng.choice(positions, size=take, replace=False).tolist())
    if len(selected) < max_cells:
        remaining = np.setdiff1d(np.arange(obs.shape[0]), np.asarray(selected), assume_unique=False)
        extra = rng.choice(remaining, size=min(max_cells - len(selected), len(remaining)), replace=False)
        selected.extend(extra.tolist())
    if len(selected) > max_cells:
        selected = rng.choice(np.asarray(selected), size=max_cells, replace=False).tolist()
    return np.asarray(sorted(selected), dtype=int)


def _standardize_obs(adata: ad.AnnData, cfg: dict) -> dict[str, str | None]:
    obs = adata.obs
    time_key = first_existing_column(obs, cfg.get("time_key_candidates", []))
    cell_type_key = first_existing_column(obs, cfg.get("cell_type_key_candidates", []))
    condition_key = first_existing_column(obs, cfg.get("condition_key_candidates", []))
    batch_key = first_existing_column(obs, cfg.get("batch_key_candidates", []))
    if time_key is None:
        raise ValueError("No time/stage/pseudotime column could be inferred for the main input.")
    if cell_type_key is None:
        cell_type_key = "unknown_lineage"
        obs[cell_type_key] = "unknown"
    numeric = pd.to_numeric(obs[time_key].astype(str).str.extract(r"([0-9]+\.?[0-9]*)", expand=False), errors="coerce")
    if numeric.isna().all():
        codes, uniques = pd.factorize(obs[time_key].astype(str), sort=True)
        numeric = pd.Series(codes.astype(float), index=obs.index)
    obs["time_numeric"] = numeric.astype(float)
    obs["time_point"] = obs[time_key].astype(str)
    obs["lineage"] = obs[cell_type_key].astype(str).replace({"nan": "unknown", "": "unknown"})
    obs["condition"] = obs[condition_key].astype(str) if condition_key else "observed"
    obs["batch"] = obs[batch_key].astype(str) if batch_key else "batch0"
    obs["is_perturbed"] = obs["condition"].str.lower().str.contains("perturb|ko|drug|stim|treat", regex=True)
    return {"time_key": time_key, "cell_type_key": cell_type_key, "condition_key": condition_key, "batch_key": batch_key}


def _split_roles(adata: ad.AnnData, cfg: dict) -> dict:
    split_mode = str(cfg.get("split_mode", "none"))
    holdout = cfg.get("holdout_time")
    role = np.full(adata.n_obs, "train", dtype=object)
    if split_mode in {"strict_time_holdout", "teacher_edge_holdout"} and holdout is not None:
        time = pd.to_numeric(adata.obs["time_numeric"], errors="coerce").to_numpy(dtype=float)
        role[np.isclose(time, float(holdout))] = "eval_holdout"
    adata.obs["split_role"] = role
    train_mask = role == "train"
    return {
        "split_mode": split_mode,
        "holdout_time": None if holdout is None else float(holdout),
        "n_train_cells": int(train_mask.sum()),
        "n_eval_holdout_cells": int((role == "eval_holdout").sum()),
    }


def _qc_and_embed(adata: ad.AnnData, cfg: dict) -> ad.AnnData:
    symbols = gene_symbols(adata)
    adata.var["swarm_gene_symbol"] = symbols
    adata.var_names_make_unique()
    x = adata.X
    totals = np.asarray(x.sum(axis=1)).ravel() if sparse.issparse(x) else np.asarray(x).sum(axis=1)
    n_genes = np.asarray((x > 0).sum(axis=1)).ravel() if sparse.issparse(x) else (np.asarray(x) > 0).sum(axis=1)
    adata.obs["n_counts"] = totals
    adata.obs["n_genes"] = n_genes
    mt_mask = pd.Series(symbols).str.lower().str.startswith(("mt-", "mt.")).to_numpy()
    if mt_mask.any():
        mt_counts = np.asarray(x[:, mt_mask].sum(axis=1)).ravel() if sparse.issparse(x) else np.asarray(x[:, mt_mask]).sum(axis=1)
        adata.obs["pct_mito"] = np.divide(mt_counts, np.maximum(totals, 1), out=np.zeros_like(totals, dtype=float), where=totals > 0)
    else:
        adata.obs["pct_mito"] = np.nan
    target_sum = float(cfg.get("normalize_total", 10000.0))
    scale = np.divide(target_sum, np.maximum(totals, 1), out=np.ones_like(totals, dtype=float), where=totals > 0)
    if sparse.issparse(adata.X):
        adata.X = adata.X.multiply(scale[:, None]).tocsr()
        adata.X.data = np.log1p(adata.X.data)
    else:
        adata.X = np.log1p(np.asarray(adata.X) * scale[:, None])
    marker_set = {g.lower() for values in CELL_CYCLE_MARKERS.values() for g in values}
    symbol_to_idx = {g.lower(): i for i, g in enumerate(symbols)}
    marker_idx = [symbol_to_idx[g] for g in marker_set if g in symbol_to_idx]
    if marker_idx:
        vals = adata.X[:, marker_idx]
        score = np.asarray(vals.mean(axis=1)).ravel() if sparse.issparse(vals) else np.asarray(vals).mean(axis=1)
        adata.obs["cell_cycle_score"] = score
    else:
        adata.obs["cell_cycle_score"] = 0.0
    train_mask = adata.obs.get("split_role", pd.Series("train", index=adata.obs_names)).astype(str).to_numpy() == "train"
    if not train_mask.any():
        train_mask = np.ones(adata.n_obs, dtype=bool)
    n_top = min(int(cfg.get("n_top_genes", 2000)), adata.n_vars)
    try:
        x_norm = adata.X
        if sparse.issparse(x_norm):
            x_fit = x_norm[train_mask]
            mean = np.asarray(x_fit.mean(axis=0)).ravel()
            mean_sq = np.asarray(x_fit.power(2).mean(axis=0)).ravel()
        else:
            arr = np.asarray(x_norm)[train_mask]
            mean = arr.mean(axis=0)
            mean_sq = (arr**2).mean(axis=0)
        var = np.maximum(mean_sq - mean**2, 0.0)
        top = np.argsort(var)[::-1][:n_top]
        keep = np.zeros(adata.n_vars, dtype=bool)
        keep[top] = True
        if marker_idx:
            keep[marker_idx] = True
        adata.var["highly_variable"] = keep
        adata = adata[:, keep].copy()
    except Exception as exc:
        logging.warning("HVG selection failed and all genes are retained: %s", exc)
    n_pcs = min(int(cfg.get("n_pcs", 30)), adata.n_vars - 1, adata.n_obs - 1)
    svd = TruncatedSVD(n_components=n_pcs, random_state=int(cfg.get("random_seed", 17)))
    svd.fit(adata.X[train_mask])
    adata.obsm["X_pca"] = svd.transform(adata.X).astype(np.float32)
    adata.uns["pca_variance_ratio"] = svd.explained_variance_ratio_.astype(float)
    if "X_umap" not in adata.obsm:
        adata.obsm["X_umap"] = adata.obsm["X_pca"][:, :2].copy()
    return adata


def _plot_overview(adata: ad.AnnData, figures_dir: str) -> None:
    ensure_dir(figures_dir)
    emb = adata.obsm["X_umap"] if "X_umap" in adata.obsm else adata.obsm["X_pca"][:, :2]
    fig, ax = plt.subplots(figsize=(7, 5), dpi=160)
    times = pd.Categorical(adata.obs["time_point"].astype(str))
    sca = ax.scatter(emb[:, 0], emb[:, 1], c=times.codes, s=4, cmap="viridis", linewidths=0, alpha=0.75)
    ax.set_title("Data overview by developmental stage")
    ax.set_xlabel("UMAP/PCA 1")
    ax.set_ylabel("UMAP/PCA 2")
    cbar = fig.colorbar(sca, ax=ax, fraction=0.035)
    cbar.set_label("stage code")
    fig.tight_layout()
    fig.savefig(Path(figures_dir) / "data_overview_umap.png")
    plt.close(fig)

    counts = adata.obs.groupby(["time_point", "lineage"], observed=False).size().reset_index(name="n")
    pivot = counts.pivot(index="time_point", columns="lineage", values="n").fillna(0)
    pivot = pivot.loc[pd.Series(pd.to_numeric(pivot.index.astype(str).str.extract(r"([0-9]+\.?[0-9]*)", expand=False), errors="coerce"), index=pivot.index).sort_values().index]
    fig, ax = plt.subplots(figsize=(8, 5), dpi=160)
    bottom = np.zeros(pivot.shape[0])
    for col in pivot.columns[:15]:
        vals = pivot[col].to_numpy()
        ax.bar(np.arange(pivot.shape[0]), vals, bottom=bottom, label=str(col))
        bottom += vals
    ax.set_xticks(np.arange(pivot.shape[0]), pivot.index.astype(str), rotation=45, ha="right")
    ax.set_ylabel("cells")
    ax.set_title("Stage × lineage cell counts")
    ax.legend(fontsize=6, ncols=2, frameon=False)
    fig.tight_layout()
    fig.savefig(Path(figures_dir) / "time_celltype_counts.png")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/data.yaml")
    parser.add_argument("--quick-fixture", action="store_true", help="Generate and preprocess a small synthetic AnnData fixture.")
    args = parser.parse_args()
    cfg = load_config(args.config)
    if args.quick_fixture:
        cfg = dict(cfg)
        cfg["quick_fixture"] = True
        cfg["input_path"] = "data/fixtures/swarmlineage_synthetic.h5ad"
        cfg["output_path"] = "processed/quick_fixture/swarmlineage_input.h5ad"
        cfg["schema_json"] = "reports/quick_fixture/data_schema.json"
        cfg["audit_report"] = "reports/quick_fixture/data_audit.md"
        cfg["loading_log"] = "logs/quick_fixture_data_loading.log"
        cfg["figures_dir"] = "figures/quick_fixture"
        cfg["max_cells"] = 450
        cfg["n_top_genes"] = 100
        cfg["n_pcs"] = 12
        cfg["holdout_time"] = 14.0
        cfg["split_mode"] = "strict_time_holdout"
        write_synthetic_fixture(cfg["input_path"], seed=int(cfg.get("random_seed", 17)))
    _setup_log(cfg.get("loading_log", "logs/data_loading.log"))
    ensure_dir("reports")
    ensure_dir("figures")
    ensure_dir("processed")
    files = discover_data_files()
    records, summary, md = audit_data_files(files)
    write_json(cfg.get("schema_json", "reports/data_schema.json"), {"summary": summary, "files": records})
    write_text(cfg.get("audit_report", "reports/data_audit.md"), md)
    logging.info("Scanned %d data files", len(files))

    adata = load_h5ad(cfg["input_path"])
    inferred = _standardize_obs(adata, cfg)
    sample_idx = _stratified_sample(adata.obs, int(cfg.get("max_cells", adata.n_obs)), int(cfg.get("random_seed", 17)), ["time_point", "lineage"])
    if sample_idx.size < adata.n_obs:
        adata = adata[sample_idx].copy()
    split_summary = _split_roles(adata, cfg)
    adata = _qc_and_embed(adata, cfg)
    adata.uns["swarmlineage_preprocess"] = {
        "input_path": cfg["input_path"],
        "inferred_columns": inferred,
        "n_cells_after_sampling": int(adata.n_obs),
        "n_genes_after_hvg": int(adata.n_vars),
        "stage_is_fallback_time": inferred["time_key"] not in {"time_numeric", "time", "day"},
        **split_summary,
    }
    _plot_overview(adata, cfg.get("figures_dir", "figures"))
    out = Path(cfg.get("output_path", "processed/swarmlineage_input.h5ad"))
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    leakage = [
        "# Leakage Audit",
        "",
        f"- split_mode: `{split_summary['split_mode']}`",
        f"- holdout_time: `{split_summary['holdout_time']}`",
        f"- preprocessing HVG/SVD fit cells: {split_summary['n_train_cells']} train cells only",
        f"- evaluation-only held-out cells present in output for scoring: {split_summary['n_eval_holdout_cells']}",
        "- held-out cells are transformed by the train-fitted SVD but do not contribute to HVG selection or SVD fitting.",
        "- teacher construction must additionally exclude `split_role == eval_holdout` cells.",
        "",
    ]
    leak_path = "reports/quick_fixture/leakage_audit.md" if cfg.get("quick_fixture") else "reports/leakage_audit.md"
    write_text(leak_path, "\n".join(leakage))
    logging.info("Wrote preprocessed AnnData to %s", out)
    print({"output": str(out), "n_obs": adata.n_obs, "n_vars": adata.n_vars, "time_key": inferred["time_key"]})


if __name__ == "__main__":
    main()
