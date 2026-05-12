from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

from src.model.cci import cci_context
from src.ot_teacher.diagnostics import coupling_lineage_edges, js_divergence
from src.utils.config import ensure_dir, load_config, write_json, write_text


def _load_npz(path: str | Path) -> dict:
    raw = np.load(path, allow_pickle=True)
    return {k: raw[k] for k in raw.files}


def _assign_unsampled(
    z: np.ndarray,
    obs: pd.DataFrame,
    sampled_idx: np.ndarray,
    values: np.ndarray,
    time_value: float,
    time_key: str,
) -> tuple[np.ndarray, np.ndarray]:
    idx_all = np.where(obs[time_key].to_numpy(dtype=float) == float(time_value))[0]
    out = np.full((idx_all.size, values.shape[1]), np.nan, dtype=float)
    sampled_positions = pd.Index(idx_all).get_indexer(sampled_idx)
    valid = sampled_positions >= 0
    out[sampled_positions[valid]] = values[valid]
    missing = np.isnan(out[:, 0])
    if missing.any() and valid.any():
        nn = NearestNeighbors(n_neighbors=1).fit(z[sampled_idx[valid]])
        nearest = nn.kneighbors(z[idx_all[missing]], return_distance=False).ravel()
        out[missing] = values[valid][nearest]
    return idx_all, out


def _plot_lineage_graph(edges: pd.DataFrame, out_path: Path) -> None:
    graph = nx.DiGraph()
    if not edges.empty:
        for row in edges.itertuples(index=False):
            src = f"{row.source_time:g}:{row.source_lineage}"
            tgt = f"{row.target_time:g}:{row.target_lineage}"
            graph.add_edge(src, tgt, weight=float(row.mass))
    fig, ax = plt.subplots(figsize=(9, 6), dpi=160)
    if graph.number_of_edges() == 0:
        ax.text(0.5, 0.5, "No lineage edges", ha="center", va="center")
    else:
        pos = nx.spring_layout(graph, seed=7, k=0.7)
        weights = np.array([graph[u][v]["weight"] for u, v in graph.edges()])
        widths = 0.5 + 5 * weights / max(weights.max(), 1e-12)
        nx.draw_networkx_nodes(graph, pos, node_size=420, node_color="#d8e6f3", edgecolors="#2b4c6f", ax=ax)
        nx.draw_networkx_edges(graph, pos, width=widths, arrows=True, alpha=0.6, edge_color="#4c6f91", ax=ax)
        nx.draw_networkx_labels(graph, pos, font_size=6, ax=ax)
        ax.set_title("OT teacher lineage graph")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def _plot_umap_scalar(adata: ad.AnnData, value: np.ndarray, title: str, out_path: Path, cmap: str = "viridis") -> None:
    emb = adata.obsm["X_umap"] if "X_umap" in adata.obsm else adata.obsm["X_pca"][:, :2]
    fig, ax = plt.subplots(figsize=(6.5, 5), dpi=160)
    sca = ax.scatter(emb[:, 0], emb[:, 1], c=value, s=4, cmap=cmap, linewidths=0, alpha=0.78)
    ax.set_title(title)
    ax.set_xlabel("UMAP/PCA 1")
    ax.set_ylabel("UMAP/PCA 2")
    fig.colorbar(sca, ax=ax, fraction=0.035)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def build_teacher(cfg: dict) -> dict:
    ensure_dir("processed")
    ensure_dir("figures")
    ensure_dir("reports")
    ensure_dir("tables")
    out_path = cfg.get("teacher_path", "processed/ot_teacher.h5ad")
    quick_outputs = "quick_fixture" in str(out_path).replace("\\", "/")
    figure_dir = ensure_dir("figures/quick_fixture") if quick_outputs else ensure_dir("figures")
    report_dir = ensure_dir("reports/quick_fixture") if quick_outputs else ensure_dir("reports")
    table_dir = ensure_dir("tables/quick_fixture") if quick_outputs else ensure_dir("tables")
    summary_path = Path(out_path).with_name("ot_teacher_summary.json") if quick_outputs else Path("processed/ot_teacher_summary.json")
    adata = ad.read_h5ad(cfg["adata_path"])
    z = np.asarray(adata.obsm[cfg.get("latent_key", "X_pca")], dtype=float)
    obs = adata.obs.copy()
    time_key = cfg.get("time_key", "time_numeric")
    cell_type_key = cfg.get("cell_type_key", "lineage")
    coupling_dir = Path(cfg.get("couplings_dir", "processed/ot_couplings"))
    index_path = Path(cfg.get("teacher_index_path", coupling_dir / "teacher_coupling_index.csv"))
    if not index_path.exists():
        legacy = coupling_dir / "moscot_coupling_index.csv"
        if legacy.exists():
            index_path = legacy
        else:
            raise FileNotFoundError(f"Missing coupling index: {index_path}. Run run_moscot first.")
    index = pd.read_csv(index_path)
    index = index[index["file"].astype(str).str.len() > 0].copy()
    if index.empty:
        raise ValueError("No teacher couplings available after held-out edge filtering.")
    terminal_time = float(pd.to_numeric(obs[time_key], errors="coerce").max())
    terminal_counts = obs.loc[pd.to_numeric(obs[time_key], errors="coerce") == terminal_time, cell_type_key].astype(str).value_counts()
    terminal_fates = terminal_counts.head(int(cfg.get("max_terminal_fates", 8))).index.tolist()
    if not terminal_fates:
        terminal_fates = sorted(obs[cell_type_key].astype(str).unique())[: int(cfg.get("max_terminal_fates", 8))]
    fate_cols = [f"fate_prob_{f}" for f in terminal_fates]
    n, d = adata.n_obs, z.shape[1]
    bary = np.full((n, d), np.nan, dtype=np.float32)
    velocity = np.zeros((n, d), dtype=np.float32)
    entropy = np.full(n, np.nan, dtype=np.float32)
    growth = np.full(n, np.nan, dtype=np.float32)
    target_time = np.full(n, np.nan, dtype=np.float32)
    fate_probs = np.zeros((n, len(terminal_fates)), dtype=np.float32)
    edge_frames = []
    pair_reports = []
    for row in index.itertuples(index=False):
        item = _load_npz(row.file)
        src = item["source_indices"].astype(int)
        tgt = item["target_indices"].astype(int)
        transition = item["transition"].astype(float)
        target_types = obs.iloc[tgt][cell_type_key].astype(str).to_numpy()
        pair_fate = np.zeros((src.size, len(terminal_fates)), dtype=float)
        for j, fate in enumerate(terminal_fates):
            pair_fate[:, j] = transition[:, target_types == fate].sum(axis=1)
        row_sum = np.maximum(pair_fate.sum(axis=1, keepdims=True), 1e-12)
        pair_fate = pair_fate / row_sum
        src_all, bary_all = _assign_unsampled(z, obs, src, item["barycentric"].astype(float), float(row.source_time), time_key)
        _, vel_all = _assign_unsampled(z, obs, src, item["barycentric"].astype(float) - z[src], float(row.source_time), time_key)
        _, fate_all = _assign_unsampled(z, obs, src, pair_fate, float(row.source_time), time_key)
        _, ent_all = _assign_unsampled(z, obs, src, item["entropy"].reshape(-1, 1), float(row.source_time), time_key)
        _, growth_all = _assign_unsampled(z, obs, src, item["growth"].reshape(-1, 1), float(row.source_time), time_key)
        bary[src_all] = bary_all
        velocity[src_all] = vel_all
        entropy[src_all] = ent_all[:, 0]
        growth[src_all] = growth_all[:, 0]
        target_time[src_all] = float(row.target_time)
        fate_probs[src_all] = fate_all
        edges = coupling_lineage_edges(item, obs, cell_type_key)
        if not edges.empty:
            edges["source_time"] = float(row.source_time)
            edges["target_time"] = float(row.target_time)
            edge_frames.append(edges)
        pair_reports.append({"source_time": row.source_time, "target_time": row.target_time, "mean_entropy": row.mean_entropy, "transport_cost": row.transport_cost})
    final_mask = pd.to_numeric(obs[time_key], errors="coerce").to_numpy(dtype=float) == terminal_time
    bary[final_mask] = z[final_mask]
    velocity[final_mask] = 0.0
    entropy[final_mask] = 0.0
    growth[final_mask] = 1.0
    target_time[final_mask] = terminal_time
    final_types = obs.loc[final_mask, cell_type_key].astype(str).to_numpy()
    for j, fate in enumerate(terminal_fates):
        fate_probs[final_mask, j] = (final_types == fate).astype(float)
    missing = ~np.isfinite(entropy)
    if missing.any():
        entropy[missing] = float(np.nanmedian(entropy)) if np.isfinite(entropy).any() else 0.5
        growth[missing] = 1.0
        bary[missing] = z[missing]
        velocity[missing] = 0.0
        fate_probs[missing] = 1.0 / len(terminal_fates)
    adata.obsm["X_ot_barycentric"] = bary
    adata.obsm["X_ot_velocity"] = velocity
    adata.obs["ot_transition_entropy"] = entropy
    adata.obs["ot_growth"] = growth
    adata.obs["ot_target_time"] = target_time
    for j, col in enumerate(fate_cols):
        adata.obs[col] = fate_probs[:, j]
    adata.obs["ot_fate_max"] = np.array(terminal_fates, dtype=object)[np.argmax(fate_probs, axis=1)]
    cci_signal, lr_pairs = cci_context(adata, obs[cell_type_key].astype(str).to_numpy())
    adata.obs["cci_signal"] = cci_signal
    adata.uns["swarmlineage_ot_teacher"] = {
        "backend": str(index.get("teacher_backend", pd.Series(["toy_sinkhorn_fallback"])).iloc[0]),
        "terminal_time": terminal_time,
        "terminal_fates": terminal_fates,
        "coupling_index": str(index_path),
        "lr_pairs_detected": [f"{a}-{b}" for a, b in lr_pairs],
        "warning": "Fate probabilities are OT-teacher pseudo-labels, not experimentally traced lineage. toy_sinkhorn_fallback is not a native moscot/WOT result.",
    }
    adata.write_h5ad(out_path)
    fate_frame = adata.obs[[time_key, cell_type_key, "ot_transition_entropy", "ot_growth", "ot_fate_max"] + fate_cols].copy()
    fate_frame.to_parquet(cfg.get("fate_probabilities_path", "processed/ot_fate_probabilities.parquet"), index=True)
    edges = pd.concat(edge_frames, ignore_index=True) if edge_frames else pd.DataFrame()
    if not edges.empty:
        edges.to_csv(table_dir / "ot_lineage_edges.csv", index=False)
    _plot_lineage_graph(edges, figure_dir / "ot_lineage_graph.png")
    _plot_umap_scalar(adata, entropy, "OT transition entropy / fate uncertainty", figure_dir / "ot_fate_umap.png", cmap="magma")
    _plot_umap_scalar(adata, growth, "OT mass expansion proxy", figure_dir / "ot_growth_map.png", cmap="viridis")

    reliability_score = float(max(0.0, min(1.0, 1.0 - np.nanmean(entropy) * 0.5 - np.nanstd(growth) * 0.1)))
    sensitivity = []
    wot_index_path = coupling_dir / "wot_coupling_index.csv"
    if wot_index_path.exists():
        wot_index = pd.read_csv(wot_index_path)
        for mos, wot in zip(index.itertuples(index=False), wot_index.itertuples(index=False)):
            m = _load_npz(mos.file)
            w = _load_npz(wot.file)
            m_mass = m["plan"].sum(axis=1)
            w_mass = w["plan"].sum(axis=1)
            size = min(m_mass.size, w_mass.size)
            sensitivity.append({"source_time": mos.source_time, "target_time": mos.target_time, "mass_js": js_divergence(m_mass[:size], w_mass[:size])})
    report = [
        "# OT Teacher Report",
        "",
        "The teacher was built from adjacent developmental-stage entropic OT couplings. If native moscot is unavailable, couplings are explicitly labelled toy_sinkhorn_fallback and cannot support high-level moscot claims.",
        "",
        f"- cells: {adata.n_obs}",
        f"- latent dimensions: {d}",
        f"- terminal pseudo-fates: {', '.join(terminal_fates)}",
        f"- mean transition entropy: {float(np.nanmean(entropy)):.4f}",
        f"- mean growth proxy: {float(np.nanmean(growth)):.4f}",
        f"- teacher backend: {adata.uns['swarmlineage_ot_teacher']['backend']}",
        f"- teacher reliability score, heuristic not a publication claim: {reliability_score:.4f}",
        "",
        "## Caveats",
        "",
        "- These couplings are OT-inferred pseudo-lineage, not true lineage tracing.",
        "- `stage_num`/Theiler stage is treated as real ordered developmental stage for this dataset; it is still coarser than dense experimental time.",
        "- Native moscot/WOT baselines should be rerun without quick fallback before any high-impact claim.",
        "",
    ]
    if sensitivity:
        sens = pd.DataFrame(sensitivity)
        sens.to_csv(table_dir / "ot_teacher_sensitivity.csv", index=False)
        report += ["## WOT-Style Sensitivity", "", sens.to_markdown(index=False), ""]
    write_text(report_dir / "ot_teacher_report.md", "\n".join(report))
    payload = {
        "teacher_path": out_path,
        "teacher_backend": adata.uns["swarmlineage_ot_teacher"]["backend"],
        "high_level_claims_allowed": adata.uns["swarmlineage_ot_teacher"]["backend"] in {"native_moscot", "native_wot"},
        "reliability_score": reliability_score,
        "terminal_fates": terminal_fates,
        "n_pairs": int(index.shape[0]),
    }
    write_json(summary_path, payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/ot_teacher.yaml")
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    cfg = load_config(args.config)
    if args.quick_fixture:
        cfg = dict(cfg)
        cfg["adata_path"] = "processed/quick_fixture/swarmlineage_input.h5ad"
        cfg["teacher_path"] = "processed/quick_fixture/ot_teacher.h5ad"
        cfg["couplings_dir"] = "processed/quick_fixture/ot_couplings"
        cfg["fate_probabilities_path"] = "processed/quick_fixture/ot_fate_probabilities.parquet"
        cfg["teacher_index_path"] = "processed/quick_fixture/ot_couplings/teacher_coupling_index.csv"
        cfg["holdout_time"] = 14.0
    result = build_teacher(cfg)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
