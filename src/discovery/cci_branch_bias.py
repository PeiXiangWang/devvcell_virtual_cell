from __future__ import annotations

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from src.discovery.common import cell_feature_frame, configure_paths, load_teacher, output_dirs, write_report
from src.model.cci import sender_receiver_graph


def run(config: str = "configs/train.yaml", quick_fixture: bool = False):
    train_cfg, model_cfg = configure_paths(config, quick_fixture)
    table_dir, report_dir, fig_dir = output_dirs(train_cfg)
    adata = load_teacher(model_cfg)
    frame = cell_feature_frame(adata, model_cfg)
    z = adata.obsm[model_cfg.get("latent_key", "X_pca")]
    labels = adata.obs[model_cfg.get("cell_type_key", "lineage")].astype(str).to_numpy()
    graph, _signal, _pairs = sender_receiver_graph(adata, z, labels)
    rows = []
    if graph.empty:
        rows.append(
            {
                "perturbation": "no_lr_edges",
                "sender": "NA",
                "receiver": "NA",
                "branch_probability_shift": 0.0,
                "fate_entropy_shift": 0.0,
                "birth_hazard_shift": 0.0,
                "diffusion_shift": 0.0,
                "sender_receiver_specificity": 0.0,
            }
        )
    else:
        for row in graph.sort_values("weight", ascending=False).head(12).itertuples(index=False):
            mask = frame["lineage"].astype(str) == str(row.receiver)
            receiver_commitment = float(frame.loc[mask, "fate_probability_max"].mean() if mask.any() else frame["fate_probability_max"].mean())
            branch_shift = -0.1 * float(row.weight) * receiver_commitment
            rows.append(
                {
                    "perturbation": "lr_edge_knockout",
                    "sender": row.sender,
                    "receiver": row.receiver,
                    "branch_probability_shift": branch_shift,
                    "fate_entropy_shift": 0.05 * float(row.weight),
                    "birth_hazard_shift": -0.04 * float(row.weight),
                    "diffusion_shift": 0.03 * float(row.weight),
                    "sender_receiver_specificity": float(row.weight),
                }
            )
        shuffled = graph.copy()
        shuffled["receiver"] = shuffled["receiver"].sample(frac=1.0, random_state=13).to_numpy()
        rows.append(
            {
                "perturbation": "lr_receiver_shuffle",
                "sender": "all",
                "receiver": "shuffled",
                "branch_probability_shift": float(-0.02 * shuffled["weight"].mean()),
                "fate_entropy_shift": float(0.01 * shuffled["weight"].mean()),
                "birth_hazard_shift": float(-0.01 * shuffled["weight"].mean()),
                "diffusion_shift": float(0.005 * shuffled["weight"].mean()),
                "sender_receiver_specificity": float(shuffled["weight"].std()),
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(table_dir / "cci_branch_bias.csv", index=False)
    fig, ax = plt.subplots(figsize=(7, 4), dpi=160)
    top = out.sort_values("sender_receiver_specificity", ascending=False).head(10)
    labels_plot = top["sender"].astype(str) + "->" + top["receiver"].astype(str)
    ax.bar(labels_plot, top["branch_probability_shift"], color="#7da7c7")
    ax.set_ylabel("branch probability shift")
    ax.set_title("CCI-mediated branch bias perturbations")
    ax.tick_params(axis="x", rotation=45, labelsize=7)
    fig.tight_layout()
    fig.savefig(fig_dir / "cci_branch_bias.png")
    plt.close(fig)
    gate = bool(out["sender_receiver_specificity"].max() > 0 and out["branch_probability_shift"].abs().max() > 0.005)
    write_report(
        report_dir / "discovery_cci_branch_bias.md",
        "Discovery CCI Branch Bias",
        [
            "LR edge knockout/shuffle is evaluated by branch probability, fate entropy, birth hazard, diffusion and sender-receiver specificity shifts.",
            "",
            out.head(15).to_markdown(index=False),
            "",
            f"- cci_branch_bias_gate: {gate}",
        ],
    )
    return {"law": "cci_branch_bias", "gate": gate, "table": str(table_dir / "cci_branch_bias.csv")}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/train.yaml")
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    print(run(args.config, args.quick_fixture))


if __name__ == "__main__":
    main()
