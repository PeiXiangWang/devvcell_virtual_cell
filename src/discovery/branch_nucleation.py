from __future__ import annotations

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

from src.discovery.common import cell_feature_frame, configure_paths, load_teacher, output_dirs, write_report


def _alignment(z: np.ndarray, v: np.ndarray, k: int = 12) -> np.ndarray:
    if z.shape[0] <= 2:
        return np.ones(z.shape[0])
    k = min(k, z.shape[0] - 1)
    neigh = NearestNeighbors(n_neighbors=k + 1).fit(z).kneighbors(z, return_distance=False)[:, 1:]
    vn = v / np.maximum(np.linalg.norm(v, axis=1, keepdims=True), 1e-8)
    return np.array([(vn[i] * vn[neigh[i]].mean(axis=0)).sum() for i in range(z.shape[0])])


def run(config: str = "configs/train.yaml", quick_fixture: bool = False):
    train_cfg, model_cfg = configure_paths(config, quick_fixture)
    table_dir, report_dir, fig_dir = output_dirs(train_cfg)
    adata = load_teacher(model_cfg)
    frame = cell_feature_frame(adata, model_cfg)
    z = np.asarray(adata.obsm[model_cfg.get("latent_key", "X_pca")], dtype=float)
    v = np.asarray(adata.obsm.get("X_ot_velocity", np.zeros_like(z)), dtype=float)
    frame["alignment_cell"] = _alignment(z, v)
    rows = []
    for time, group in frame.groupby("time_numeric", observed=False):
        labels = group["lineage"].astype(str)
        counts = labels.value_counts(normalize=True)
        imbalance = float(counts.max() - counts.min()) if len(counts) > 1 else 1.0
        rows.append(
            {
                "time_numeric": float(time),
                "local_velocity_alignment_A": float(group["alignment_cell"].mean()),
                "branch_cohesion_C": float(1.0 - group.groupby("lineage")["local_density"].std().fillna(0).mean()),
                "lineage_separation_S": float(group.groupby("lineage")["ot_displacement"].mean().std()),
                "fate_entropy_H": float(group["fate_entropy"].mean()),
                "branch_imbalance_B": imbalance,
                "n_cells": int(group.shape[0]),
            }
        )
    out = pd.DataFrame(rows).sort_values("time_numeric")
    out.to_csv(table_dir / "swarm_order_parameters.csv", index=False)
    fig, ax = plt.subplots(figsize=(7, 4), dpi=160)
    for col in ["local_velocity_alignment_A", "branch_cohesion_C", "lineage_separation_S", "fate_entropy_H", "branch_imbalance_B"]:
        vals = out[col].to_numpy(dtype=float)
        vals = (vals - np.nanmin(vals)) / max(np.nanmax(vals) - np.nanmin(vals), 1e-8)
        ax.plot(out["time_numeric"], vals, marker="o", label=col)
    ax.set_xlabel("developmental time")
    ax.set_ylabel("normalized order parameter")
    ax.set_title("Branch nucleation order parameters")
    ax.legend(fontsize=6, frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "branch_nucleation_order_parameters.png")
    plt.close(fig)
    gate = bool(out["local_velocity_alignment_A"].std() > 1e-3 and out["fate_entropy_H"].std() > 1e-3)
    write_report(
        report_dir / "discovery_branch_nucleation.md",
        "Discovery Branch Nucleation",
        [
            "Order parameters track velocity alignment, branch cohesion, lineage separation, fate entropy and branch imbalance over developmental time.",
            "",
            out.to_markdown(index=False),
            "",
            f"- branch_nucleation_gate: {gate}",
        ],
    )
    return {"law": "branch_nucleation", "gate": gate, "table": str(table_dir / "swarm_order_parameters.csv")}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/train.yaml")
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    print(run(args.config, args.quick_fixture))


if __name__ == "__main__":
    main()

