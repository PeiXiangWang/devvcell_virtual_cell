"""Summarize developmental impact directly from existing OT+GRN outputs.

This is an interpretation layer over the original RDEG artifacts. It does not
train a new predictor.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data" / "rdeg_neural_cell_mvp"
RESULTS_DIR = PROJECT_ROOT / "results"
TABLES_DIR = RESULTS_DIR / "tables"
FIGURES_DIR = RESULTS_DIR / "figures"


SYSTEM_PROCESS = {
    "neural cell": "neural patterning / neuronal differentiation",
    "myoblast": "myogenic commitment and skeletal muscle differentiation",
    "myotube": "muscle fiber maturation",
    "muscle precursor cell": "mesoderm-to-muscle precursor specification",
    "cardiac muscle cell": "cardiac muscle differentiation",
    "erythroid progenitor cell": "definitive erythroid differentiation",
    "primitive erythroid progenitor": "primitive erythroid differentiation",
    "hematopoietic stem cell": "hematopoietic progenitor maintenance and lineage choice",
}

PROCESS_MARKERS = {
    "Wnt/FGF signaling": {"Lef1", "Tcf7", "Wnt1", "Wnt8b", "Fgf15"},
    "neural adhesion and axon guidance": {"Cntnap2", "Epha5", "Epha3", "Cdh8", "Cdh10", "Cadps2"},
    "neural regionalization": {"Rfx4", "Otx2", "Pax3", "Hoxd4", "Zbtb16", "Dach2", "Ebf1", "Ebf2", "Ebf3"},
    "myogenic / mesoderm morphogenesis": {"Tbx15", "Myl2", "Mecom", "Adam12", "Svep1", "Angpt1", "Smoc1", "Frem1"},
    "hematopoietic / erythroid regulation": {"Tal1", "Erg", "Foxo3", "Rgcc", "Ifi203"},
    "metabolic or stress response": {"Foxo3", "Acoxl", "Pde3a", "Rgcc"},
    "ECM and tissue morphogenesis": {"Frem1", "Smoc1", "Svep1", "Angpt1", "Adam12"},
}


def read_inputs() -> dict[str, pd.DataFrame]:
    return {
        "grn": pd.read_csv(DATA_DIR / "grn_learned_network.csv"),
        "system_edges": pd.read_csv(DATA_DIR / "system_specific_edges.csv"),
        "tf_ko": pd.read_csv(DATA_DIR / "tf_knockout_results.csv"),
        "temporal": pd.read_csv(DATA_DIR / "temporal_sensitivity.csv"),
        "fate_mi": pd.read_csv(DATA_DIR / "tf_fate_mutual_info.csv"),
        "edge_importance": pd.read_csv(DATA_DIR / "edge_importance_scores.csv"),
        "hub": pd.read_csv(DATA_DIR / "hub_tf_ranking.csv"),
    }


def minmax(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").astype(float)
    lo = values.min(skipna=True)
    hi = values.max(skipna=True)
    if not np.isfinite(lo) or not np.isfinite(hi) or np.isclose(lo, hi):
        return pd.Series(np.zeros(len(values)), index=values.index)
    return (values - lo) / (hi - lo)


def signed_target_list(frame: pd.DataFrame, top_n: int = 8) -> str:
    rows = []
    for _, row in frame.sort_values("abs_weight", ascending=False).head(top_n).iterrows():
        sign = "+" if row["weight"] >= 0 else "-"
        rows.append(f"{row['target']}({sign}{abs(float(row['weight'])):.3f})")
    return "; ".join(rows)


def target_process_tags(targets: list[str], system: str) -> str:
    tags = []
    target_set = set(targets)
    for label, markers in PROCESS_MARKERS.items():
        if target_set & markers:
            tags.append(label)

    system_process = SYSTEM_PROCESS.get(system)
    if system_process:
        tags.insert(0, system_process)

    if not tags:
        tags.append("system-specific developmental regulation")

    seen = []
    for tag in tags:
        if tag not in seen:
            seen.append(tag)
    return "; ".join(seen)


def summarize_grn(grn: pd.DataFrame) -> pd.DataFrame:
    grn = grn.copy()
    grn["abs_weight"] = pd.to_numeric(grn["abs_weight"], errors="coerce").fillna(0.0)
    rows = []
    for tf, sub in grn.groupby("source_tf"):
        top = sub.sort_values("abs_weight", ascending=False).head(8)
        rows.append(
            {
                "tf": tf,
                "grn_n_direct_targets": int(sub["target_gene"].nunique()),
                "grn_mean_abs_weight": float(sub["abs_weight"].mean()),
                "grn_max_abs_weight": float(sub["abs_weight"].max()),
                "top_grn_targets": "; ".join(
                    f"{r.target_gene}({r.regulation_type},{float(r.abs_weight):.3f})"
                    for r in top.itertuples(index=False)
                ),
            }
        )
    return pd.DataFrame(rows)


def summarize_temporal(temporal: pd.DataFrame) -> pd.DataFrame:
    temporal = temporal.copy()
    temporal["sensitivity"] = pd.to_numeric(temporal["sensitivity"], errors="coerce").fillna(0.0)
    idx = temporal.groupby("tf_name")["sensitivity"].idxmax()
    peak = temporal.loc[idx, ["tf_name", "stage", "sensitivity"]].rename(
        columns={"tf_name": "tf", "stage": "peak_sensitivity_stage", "sensitivity": "peak_sensitivity"}
    )
    mean = temporal.groupby("tf_name", as_index=False).agg(mean_temporal_sensitivity=("sensitivity", "mean"))
    mean = mean.rename(columns={"tf_name": "tf"})
    return mean.merge(peak, on="tf", how="left")


def build_tf_system_impact(inputs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    system_edges = inputs["system_edges"].copy()
    system_edges["weight"] = pd.to_numeric(system_edges["weight"], errors="coerce").fillna(0.0)
    system_edges["abs_weight"] = system_edges["weight"].abs()

    rows = []
    for (tf, system), sub in system_edges.groupby(["tf", "system"]):
        targets = sub.sort_values("abs_weight", ascending=False)["target"].astype(str).tolist()
        pos = int((sub["weight"] > 0).sum())
        neg = int((sub["weight"] < 0).sum())
        rows.append(
            {
                "tf": tf,
                "system_or_cell_type": system,
                "system_process": SYSTEM_PROCESS.get(system, "system-specific developmental regulation"),
                "n_system_specific_edges": int(len(sub)),
                "sum_abs_system_weight": float(sub["abs_weight"].sum()),
                "mean_abs_system_weight": float(sub["abs_weight"].mean()),
                "max_abs_system_weight": float(sub["abs_weight"].max()),
                "n_positive_edges": pos,
                "n_negative_edges": neg,
                "top_system_targets": signed_target_list(sub),
                "biological_process_proxy": target_process_tags(targets, system),
            }
        )

    impact = pd.DataFrame(rows)

    grn_summary = summarize_grn(inputs["grn"])
    temporal_summary = summarize_temporal(inputs["temporal"])

    tf_ko = inputs["tf_ko"].rename(columns={"tf_name": "tf"}).copy()
    for col in ["global_effect_score", "developmental_delay_score", "mass_shift_mean"]:
        tf_ko[col] = pd.to_numeric(tf_ko[col], errors="coerce")

    fate_mi = inputs["fate_mi"].rename(columns={"tf_name": "tf"}).copy()
    fate_mi["mutual_information"] = pd.to_numeric(fate_mi["mutual_information"], errors="coerce")

    hub = inputs["hub"].copy()
    hub["mean_out_degree"] = pd.to_numeric(hub["mean_out_degree"], errors="coerce")

    impact = impact.merge(grn_summary, on="tf", how="left")
    impact = impact.merge(tf_ko[["tf", "global_effect_score", "developmental_delay_score", "mass_shift_mean"]], on="tf", how="left")
    impact = impact.merge(temporal_summary, on="tf", how="left")
    impact = impact.merge(fate_mi[["tf", "mutual_information", "p_value", "significant"]], on="tf", how="left")
    impact = impact.merge(hub, on="tf", how="left")

    for col in [
        "sum_abs_system_weight",
        "global_effect_score",
        "developmental_delay_score",
        "mean_temporal_sensitivity",
        "mutual_information",
        "mean_out_degree",
        "grn_mean_abs_weight",
    ]:
        impact[f"{col}_norm"] = minmax(impact[col].fillna(0.0))

    # This is an evidence-integration rank over existing OT+GRN outputs, not a new predictor.
    impact["ot_grn_evidence_score"] = (
        0.35 * impact["sum_abs_system_weight_norm"]
        + 0.20 * impact["global_effect_score_norm"]
        + 0.15 * impact["mean_temporal_sensitivity_norm"]
        + 0.15 * impact["mutual_information_norm"]
        + 0.10 * impact["mean_out_degree_norm"]
        + 0.05 * impact["grn_mean_abs_weight_norm"]
    )

    ordered_cols = [
        "tf",
        "system_or_cell_type",
        "system_process",
        "ot_grn_evidence_score",
        "n_system_specific_edges",
        "sum_abs_system_weight",
        "mean_abs_system_weight",
        "max_abs_system_weight",
        "n_positive_edges",
        "n_negative_edges",
        "top_system_targets",
        "biological_process_proxy",
        "grn_n_direct_targets",
        "grn_mean_abs_weight",
        "top_grn_targets",
        "global_effect_score",
        "developmental_delay_score",
        "mean_temporal_sensitivity",
        "peak_sensitivity_stage",
        "peak_sensitivity",
        "mutual_information",
        "p_value",
        "significant",
        "mean_out_degree",
        "mass_shift_mean",
    ]
    return impact[ordered_cols].sort_values("ot_grn_evidence_score", ascending=False).reset_index(drop=True)


def build_tf_summary(impact: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for tf, sub in impact.groupby("tf"):
        top_systems = sub.sort_values("ot_grn_evidence_score", ascending=False).head(4)
        process_tags = []
        for text in top_systems["biological_process_proxy"].astype(str):
            for part in text.split("; "):
                if part and part not in process_tags:
                    process_tags.append(part)
        first = top_systems.iloc[0]
        rows.append(
            {
                "tf": tf,
                "n_implicated_systems": int(sub["system_or_cell_type"].nunique()),
                "total_system_edges": int(sub["n_system_specific_edges"].sum()),
                "total_abs_system_weight": float(sub["sum_abs_system_weight"].sum()),
                "best_system_or_cell_type": first["system_or_cell_type"],
                "best_system_process": first["system_process"],
                "best_system_evidence_score": float(first["ot_grn_evidence_score"]),
                "top_systems": "; ".join(top_systems["system_or_cell_type"].astype(str).tolist()),
                "process_summary": "; ".join(process_tags[:6]),
                "top_grn_targets": first.get("top_grn_targets", ""),
                "global_effect_score": first.get("global_effect_score", np.nan),
                "developmental_delay_score": first.get("developmental_delay_score", np.nan),
                "peak_sensitivity_stage": first.get("peak_sensitivity_stage", ""),
                "mutual_information": first.get("mutual_information", np.nan),
                "fate_mi_significant": first.get("significant", ""),
            }
        )
    return pd.DataFrame(rows).sort_values("best_system_evidence_score", ascending=False).reset_index(drop=True)


def build_system_summary(impact: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for system, sub in impact.groupby("system_or_cell_type"):
        top = sub.sort_values("ot_grn_evidence_score", ascending=False).head(8)
        rows.append(
            {
                "system_or_cell_type": system,
                "system_process": SYSTEM_PROCESS.get(system, "system-specific developmental regulation"),
                "n_tfs": int(sub["tf"].nunique()),
                "n_system_specific_edges": int(sub["n_system_specific_edges"].sum()),
                "sum_abs_system_weight": float(sub["sum_abs_system_weight"].sum()),
                "top_tfs": "; ".join(top["tf"].astype(str).tolist()),
                "top_tf_target_examples": " | ".join(
                    f"{r.tf}: {str(r.top_system_targets).split('; ')[0]}"
                    for r in top.itertuples(index=False)
                ),
            }
        )
    return pd.DataFrame(rows).sort_values("sum_abs_system_weight", ascending=False).reset_index(drop=True)


def plot_system_heatmap(impact: pd.DataFrame) -> None:
    top_tfs = impact.groupby("tf")["ot_grn_evidence_score"].max().sort_values(ascending=False).head(16).index
    pivot = (
        impact.loc[impact["tf"].isin(top_tfs)]
        .pivot_table(index="tf", columns="system_or_cell_type", values="sum_abs_system_weight", aggfunc="sum", fill_value=0.0)
        .loc[top_tfs]
    )
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    im = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="mako" if "mako" in plt.colormaps() else "viridis")
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.set_title("OT+GRN system-specific regulatory weight")
    ax.set_xlabel("System / cell type")
    ax.set_ylabel("TF")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label="sum abs system weight")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "ot_grn_tf_system_heatmap.png", dpi=220)
    plt.close(fig)


def plot_top_tf_bar(tf_summary: pd.DataFrame) -> None:
    top = tf_summary.head(14).iloc[::-1]
    fig, ax = plt.subplots(figsize=(9.0, 5.5))
    ax.barh(top["tf"], top["best_system_evidence_score"], color="#456990")
    ax.set_xlabel("OT+GRN evidence score")
    ax.set_title("Top TF developmental impact by original OT+GRN evidence")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "ot_grn_top_tf_developmental_impact.png", dpi=220)
    plt.close(fig)


def write_summary_json(tf_summary: pd.DataFrame, system_summary: pd.DataFrame, impact: pd.DataFrame) -> None:
    summary = {
        "analysis": "original_ot_grn_developmental_impact",
        "note": "Evidence integration over existing OT+GRN artifacts; no new predictor was trained.",
        "n_tf_system_rows": int(len(impact)),
        "n_tfs": int(tf_summary["tf"].nunique()),
        "n_systems": int(system_summary["system_or_cell_type"].nunique()),
        "top_tf_rows": tf_summary.head(8).to_dict(orient="records"),
        "top_system_rows": system_summary.head(8).to_dict(orient="records"),
        "outputs": {
            "tf_system_impact": "results/tables/ot_grn_tf_system_developmental_impact.csv",
            "tf_summary": "results/tables/ot_grn_tf_developmental_summary.csv",
            "system_summary": "results/tables/ot_grn_system_developmental_summary.csv",
        },
    }
    with (TABLES_DIR / "ot_grn_developmental_impact_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    inputs = read_inputs()
    impact = build_tf_system_impact(inputs)
    tf_summary = build_tf_summary(impact)
    system_summary = build_system_summary(impact)

    impact.to_csv(TABLES_DIR / "ot_grn_tf_system_developmental_impact.csv", index=False)
    tf_summary.to_csv(TABLES_DIR / "ot_grn_tf_developmental_summary.csv", index=False)
    system_summary.to_csv(TABLES_DIR / "ot_grn_system_developmental_summary.csv", index=False)
    write_summary_json(tf_summary, system_summary, impact)
    plot_system_heatmap(impact)
    plot_top_tf_bar(tf_summary)

    print("OT+GRN developmental impact analysis completed.")
    print(f"Top TF: {tf_summary.iloc[0]['tf']} -> {tf_summary.iloc[0]['best_system_or_cell_type']}")
    print(f"Top system: {system_summary.iloc[0]['system_or_cell_type']}")


if __name__ == "__main__":
    main()
