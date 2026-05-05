"""Run a TF/GRN stimulus-response head in the cell-level latent space.

This script is deliberately framed as a zero-shot perturbation head: it uses
existing TF priority, temporal sensitivity and GRN edge evidence to construct
TF knockdown directions, then projects those directions into the saved cell
state latent space. It does not claim wet-lab perturbation validation.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import anndata as ad
import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from devvcell.io import load_json, resolve_project_path, write_json  # noqa: E402
from devvcell.stages import canonical_stage, stage_number  # noqa: E402


plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial Unicode MS", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/cell_level_baseline.json")
    parser.add_argument(
        "--model-results-dir",
        default=None,
        help="Directory containing saved cell-level model artifacts. Defaults to config results_dir.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory where stimulus outputs are written. Defaults to config results_dir.",
    )
    return parser.parse_args()


def as_csr_float32(matrix) -> sparse.csr_matrix:
    if sparse.issparse(matrix):
        return matrix.tocsr().astype(np.float32)
    return sparse.csr_matrix(np.asarray(matrix, dtype=np.float32))


def minmax(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").astype(float)
    lo = values.min(skipna=True)
    hi = values.max(skipna=True)
    if not np.isfinite(lo) or not np.isfinite(hi) or math.isclose(lo, hi):
        return pd.Series(np.zeros(len(values)), index=series.index)
    return (values - lo) / (hi - lo)


def gene_names(adata: ad.AnnData) -> list[str]:
    if "gene_short_name" in adata.var.columns:
        return adata.var["gene_short_name"].astype(str).tolist()
    if "feature_name" in adata.var.columns:
        return adata.var["feature_name"].astype(str).tolist()
    return adata.var_names.astype(str).tolist()


def latent_delta_from_gene_delta(gene_delta: np.ndarray, svd, scaler) -> np.ndarray:
    scale = np.asarray(scaler.scale_, dtype=np.float32).copy()
    scale[scale == 0] = 1.0
    latent_delta = np.matmul(gene_delta.astype(np.float32), svd.components_.astype(np.float32).T) / scale
    norm = float(np.linalg.norm(latent_delta))
    if not np.isfinite(norm) or norm == 0:
        return np.zeros_like(latent_delta, dtype=np.float32)
    return (latent_delta / norm).astype(np.float32)


def build_tf_gene_delta(
    tf: str,
    system_name: str,
    config: dict,
    gene_to_idx: dict[str, int],
    n_genes: int,
    grn: pd.DataFrame,
    system_edges: pd.DataFrame,
) -> tuple[np.ndarray, dict[str, object]]:
    stim_cfg = config["stimulus"]
    top_targets = int(stim_cfg["top_target_genes"])
    gene_delta = np.zeros(n_genes, dtype=np.float32)

    tf_grn = grn.loc[grn["source_tf"] == tf].copy()
    tf_grn["abs_weight"] = pd.to_numeric(tf_grn["abs_weight"], errors="coerce").fillna(0.0)
    tf_grn["weight"] = pd.to_numeric(tf_grn["weight"], errors="coerce").fillna(0.0)
    tf_grn = tf_grn.sort_values("abs_weight", ascending=False).head(top_targets)

    matched_global = 0
    for row in tf_grn.itertuples(index=False):
        idx = gene_to_idx.get(str(row.target_gene))
        if idx is None:
            continue
        # TF knockdown: activator targets decrease and repressed targets increase.
        gene_delta[idx] += -float(row.weight) * float(stim_cfg["global_grn_weight"])
        matched_global += 1

    mapped_systems = set(stim_cfg["system_map"].get(system_name, [system_name]))
    tf_system = system_edges.loc[
        (system_edges["tf"] == tf) & (system_edges["system"].isin(mapped_systems))
    ].copy()
    tf_system["weight"] = pd.to_numeric(tf_system["weight"], errors="coerce").fillna(0.0)
    tf_system["abs_weight"] = tf_system["weight"].abs()
    tf_system = tf_system.sort_values("abs_weight", ascending=False).head(top_targets)

    matched_system = 0
    for row in tf_system.itertuples(index=False):
        idx = gene_to_idx.get(str(row.target))
        if idx is None:
            continue
        gene_delta[idx] += -float(row.weight) * float(stim_cfg["system_edge_weight"])
        matched_system += 1

    norm = float(np.linalg.norm(gene_delta))
    if norm > 0:
        gene_delta = gene_delta / norm

    info = {
        "matched_global_grn_targets": matched_global,
        "matched_system_targets": matched_system,
        "gene_delta_norm": norm,
        "mapped_systems": ";".join(sorted(mapped_systems)),
    }
    return gene_delta.astype(np.float32), info


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom == 0 or not np.isfinite(denom):
        return float("nan")
    return float(np.dot(a, b) / denom)


def read_inputs(results_dir: Path) -> dict[str, pd.DataFrame]:
    data_dir = PROJECT_ROOT / "data" / "rdeg_neural_cell_mvp"
    return {
        "priority": pd.read_csv(PROJECT_ROOT / "results" / "tables" / "perturbation_priority.csv"),
        "stage_vulnerability": pd.read_csv(PROJECT_ROOT / "results" / "tables" / "stage_vulnerability.csv"),
        "temporal": pd.read_csv(data_dir / "temporal_sensitivity.csv"),
        "grn": pd.read_csv(data_dir / "grn_learned_network.csv"),
        "system_edges": pd.read_csv(data_dir / "system_specific_edges.csv"),
    }


def prepare_temporal(temporal: pd.DataFrame) -> pd.DataFrame:
    temporal = temporal.copy()
    temporal["stage"] = temporal["stage"].map(canonical_stage)
    temporal["sensitivity"] = pd.to_numeric(temporal["sensitivity"], errors="coerce").fillna(0.0)
    temporal["temporal_sensitivity_norm"] = minmax(temporal["sensitivity"])
    return temporal


def make_heatmap(table: pd.DataFrame, figures_dir: Path) -> None:
    pivot = (
        table.groupby(["tf_name", "devvcell_system"])["stimulus_response_norm"]
        .mean()
        .reset_index()
        .pivot(index="tf_name", columns="devvcell_system", values="stimulus_response_norm")
        .fillna(0.0)
    )
    order = table.groupby("tf_name")["stimulus_response_norm"].mean().sort_values(ascending=False).index
    pivot = pivot.loc[order]

    fig, ax = plt.subplots(figsize=(7.4, max(4.2, 0.32 * len(pivot))))
    im = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="magma")
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_xlabel("DevVCell 系统")
    ax.set_ylabel("TF 敲低")
    ax.set_title("细胞级 TF/GRN 刺激响应 proxy")
    fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, label="响应强度")
    fig.tight_layout()
    fig.savefig(figures_dir / "cell_level_tf_grn_stimulus_heatmap.png", dpi=220)
    plt.close(fig)


def make_recovery_scatter(table: pd.DataFrame, figures_dir: Path) -> None:
    summary = (
        table.groupby("tf_name", as_index=False)
        .agg(
            mean_response=("stimulus_response_norm", "mean"),
            mean_recovery_probability=("cell_level_recovery_probability_proxy", "mean"),
            mean_fate_displacement=("fate_displacement_from_baseline", "mean"),
        )
        .sort_values("mean_response", ascending=False)
    )

    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    ax.scatter(
        summary["mean_response"],
        summary["mean_recovery_probability"],
        s=55 + 260 * minmax(summary["mean_fate_displacement"]),
        c=summary["mean_fate_displacement"],
        cmap="viridis",
        edgecolor="white",
        linewidth=0.6,
        alpha=0.9,
    )
    for row in summary.head(8).itertuples(index=False):
        ax.text(row.mean_response, row.mean_recovery_probability, row.tf_name, fontsize=8)
    ax.set_xlabel("平均刺激响应强度")
    ax.set_ylabel("细胞级恢复概率 proxy")
    ax.set_title("TF/GRN 刺激响应与恢复 proxy")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(figures_dir / "cell_level_tf_recovery_scatter.png", dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    config = load_json(args.config)
    seed = int(config["seed"])
    rng = np.random.default_rng(seed)

    subset_path = resolve_project_path(config["processed_subset_path"])
    model_results_dir = resolve_project_path(args.model_results_dir or config["results_dir"])
    results_dir = resolve_project_path(args.output_dir or config["results_dir"])
    model_dir = model_results_dir / "models"
    tables_dir = results_dir / "tables"
    figures_dir = results_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    svd = joblib.load(model_dir / "state_svd.joblib")
    scaler = joblib.load(model_dir / "state_scaler.joblib")
    ridge = joblib.load(model_dir / "transition_ridge.joblib")

    adata = ad.read_h5ad(subset_path)
    X = as_csr_float32(adata.X)
    obs = adata.obs.reset_index(drop=True).copy()
    latent = scaler.transform(svd.transform(X)).astype(np.float32)
    names = gene_names(adata)
    gene_to_idx = {name: idx for idx, name in enumerate(names)}

    inputs = read_inputs(results_dir)
    priority = inputs["priority"].head(int(config["stimulus"]["top_n_tfs"])).copy()
    priority_numeric = [
        "devvcell_priority_score",
        "response_amplitude_proxy",
        "fate_displacement_proxy",
        "feedback_cost",
        "recovery_probability_proxy",
    ]
    for col in priority_numeric:
        priority[col] = pd.to_numeric(priority[col], errors="coerce").fillna(0.0)

    temporal = prepare_temporal(inputs["temporal"])
    vulnerability = inputs["stage_vulnerability"].copy()
    vulnerability["stage"] = vulnerability["stage"].map(canonical_stage)
    vulnerability["vulnerability_score"] = pd.to_numeric(vulnerability["vulnerability_score"], errors="coerce").fillna(0.0)
    vulnerability_map = dict(zip(vulnerability["stage"], vulnerability["vulnerability_score"]))

    max_cells = int(config["stimulus"]["max_cells_per_stage_system"])
    rows: list[dict[str, object]] = []
    tf_delta_cache: dict[tuple[str, str], tuple[np.ndarray, dict[str, object]]] = {}

    for tf_row in priority.itertuples(index=False):
        tf = str(tf_row.tf_name)
        tf_temporal = temporal.loc[temporal["tf_name"] == tf, ["stage", "temporal_sensitivity_norm"]]
        temporal_map = dict(zip(tf_temporal["stage"], tf_temporal["temporal_sensitivity_norm"]))

        for system_name in sorted(obs["devvcell_system"].astype(str).unique()):
            if system_name == "unassigned":
                continue
            cache_key = (tf, system_name)
            if cache_key not in tf_delta_cache:
                gene_delta, target_info = build_tf_gene_delta(
                    tf,
                    system_name,
                    config,
                    gene_to_idx,
                    int(adata.n_vars),
                    inputs["grn"],
                    inputs["system_edges"],
                )
                latent_direction = latent_delta_from_gene_delta(gene_delta, svd, scaler)
                tf_delta_cache[cache_key] = (latent_direction, target_info)
            else:
                latent_direction, target_info = tf_delta_cache[cache_key]

            if np.linalg.norm(latent_direction) == 0:
                continue

            system_mask = obs["devvcell_system"].astype(str) == system_name
            for stage_num in sorted(obs.loc[system_mask, "stage_num"].astype(int).unique()):
                tgt_stage_num = int(stage_num) + 1
                if tgt_stage_num not in set(obs.loc[system_mask, "stage_num"].astype(int)):
                    continue
                stage_label = f"Theiler stage {int(stage_num)}"
                idx = obs.index[system_mask & (obs["stage_num"].astype(int) == int(stage_num))].to_numpy()
                if len(idx) == 0:
                    continue
                if len(idx) > max_cells:
                    idx = rng.choice(idx, size=max_cells, replace=False)

                x = latent[idx]
                baseline = ridge.predict(x).astype(np.float32)
                transition = baseline - x
                mean_transition = transition.mean(axis=0)
                transition_norm = float(np.linalg.norm(mean_transition))
                if transition_norm == 0 or not np.isfinite(transition_norm):
                    transition_norm = float(np.mean(np.linalg.norm(transition, axis=1)))

                temporal_modifier = float(temporal_map.get(stage_label, temporal["temporal_sensitivity_norm"].mean()))
                vulnerability_modifier = float(vulnerability_map.get(stage_label, vulnerability["vulnerability_score"].mean()))
                response_proxy = float(tf_row.response_amplitude_proxy)
                grn_evidence_scale = 0.5 + math.sqrt(max(0.0, float(target_info.get("gene_delta_norm", 0.0))))
                strength = (
                    float(config["stimulus"]["stimulus_strength"])
                    * transition_norm
                    * (0.4 + response_proxy)
                    * (0.6 + temporal_modifier)
                    * (0.7 + vulnerability_modifier)
                    * grn_evidence_scale
                )
                perturbation_delta = latent_direction * strength
                perturbed = baseline + perturbation_delta[None, :]

                response_norm = float(np.mean(np.linalg.norm(perturbed - baseline, axis=1)))
                fate_displacement = float(np.linalg.norm(perturbed.mean(axis=0) - baseline.mean(axis=0)))
                recovery_cost = response_norm * float(tf_row.feedback_cost)
                recovery_probability = float(np.exp(-min(8.0, 4.0 * recovery_cost)))
                alignment = cosine(latent_direction, mean_transition)

                rows.append(
                    {
                        "tf_name": tf,
                        "devvcell_system": system_name,
                        "stage": stage_label,
                        "stage_num": int(stage_num),
                        "target_stage": f"Theiler stage {tgt_stage_num}",
                        "n_cells": int(len(idx)),
                        "stimulus_mode": "tf_knockdown_zero_shot_grn_projection",
                        "stimulus_response_norm": response_norm,
                        "fate_displacement_from_baseline": fate_displacement,
                        "cell_level_feedback_cost_proxy": recovery_cost,
                        "cell_level_recovery_probability_proxy": recovery_probability,
                        "alignment_with_normal_transition": alignment,
                        "transition_norm": transition_norm,
                        "stimulus_strength": strength,
                        "grn_evidence_scale": grn_evidence_scale,
                        "temporal_sensitivity_norm": temporal_modifier,
                        "stage_vulnerability_score": vulnerability_modifier,
                        "devvcell_priority_score": float(tf_row.devvcell_priority_score),
                        "response_amplitude_proxy": response_proxy,
                        "feedback_cost": float(tf_row.feedback_cost),
                        **target_info,
                    }
                )

    response = pd.DataFrame(rows).sort_values("stimulus_response_norm", ascending=False)
    if response.empty:
        raise RuntimeError("No stimulus response rows were produced.")

    response["stimulus_response_rank"] = np.arange(1, len(response) + 1)
    response.to_csv(tables_dir / "cell_level_tf_grn_stimulus_response.csv", index=False)

    tf_summary = (
        response.groupby("tf_name", as_index=False)
        .agg(
            mean_stimulus_response_norm=("stimulus_response_norm", "mean"),
            max_stimulus_response_norm=("stimulus_response_norm", "max"),
            mean_fate_displacement=("fate_displacement_from_baseline", "mean"),
            mean_recovery_probability=("cell_level_recovery_probability_proxy", "mean"),
            mean_alignment_with_normal_transition=("alignment_with_normal_transition", "mean"),
            n_system_stage_tests=("stimulus_response_norm", "size"),
            matched_global_grn_targets=("matched_global_grn_targets", "max"),
            matched_system_targets=("matched_system_targets", "max"),
        )
        .sort_values("mean_stimulus_response_norm", ascending=False)
    )
    tf_summary.to_csv(tables_dir / "cell_level_tf_grn_stimulus_summary.csv", index=False)

    make_heatmap(response, figures_dir)
    make_recovery_scatter(response, figures_dir)

    def rel(path: Path) -> str:
        try:
            return str(path.relative_to(PROJECT_ROOT)).replace("\\", "/")
        except ValueError:
            return str(path)

    summary = {
        "analysis": "cell_level_tf_grn_stimulus_response_head",
        "note": "Zero-shot GRN projection in the saved cell-level latent space; no external perturbation labels were used.",
        "n_rows": int(len(response)),
        "n_tfs": int(response["tf_name"].nunique()),
        "n_systems": int(response["devvcell_system"].nunique()),
        "top_rows": response.head(10).to_dict(orient="records"),
        "top_tf_summary": tf_summary.head(10).to_dict(orient="records"),
        "outputs": {
            "response_table": rel(tables_dir / "cell_level_tf_grn_stimulus_response.csv"),
            "tf_summary": rel(tables_dir / "cell_level_tf_grn_stimulus_summary.csv"),
            "heatmap": rel(figures_dir / "cell_level_tf_grn_stimulus_heatmap.png"),
            "recovery_scatter": rel(figures_dir / "cell_level_tf_recovery_scatter.png"),
        },
    }
    write_json(results_dir / "stimulus_response_summary.json", summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
