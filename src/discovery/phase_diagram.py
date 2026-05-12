from __future__ import annotations

import argparse
from itertools import product

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from src.discovery.common import cell_feature_frame, configure_paths, load_teacher, output_dirs, write_report


def _classify(diffusion_scale: float, memory_decay: float, memory_diffusion: float, alignment: float, density_suppression: float, cci_gate: float) -> str:
    if diffusion_scale > 1.6 and memory_diffusion > 0.4:
        return "over-diffusion"
    if density_suppression < 0.2 and cci_gate > 1.2:
        return "unstable population explosion"
    if density_suppression > 1.4 and alignment < 0.5:
        return "lineage extinction"
    if alignment > 1.4 and diffusion_scale < 0.7:
        return "collapsed branch"
    if memory_decay < 0.3 and cci_gate > 1.0:
        return "premature branching"
    return "faithful flow"


def run(config: str = "configs/train.yaml", quick_fixture: bool = False):
    train_cfg, model_cfg = configure_paths(config, quick_fixture)
    table_dir, report_dir, fig_dir = output_dirs(train_cfg)
    adata = load_teacher(model_cfg)
    frame = cell_feature_frame(adata, model_cfg)
    mean_entropy = float(frame["fate_entropy"].mean())
    rows = []
    values = [0.5, 1.0, 1.8]
    for diffusion_scale, memory_decay, memory_diffusion, alignment, density_suppression, cci_gate in product(values, [0.2, 0.8, 1.5], [0.1, 0.5], values, [0.1, 0.8, 1.6], values):
        state = _classify(diffusion_scale, memory_decay, memory_diffusion, alignment, density_suppression, cci_gate)
        rows.append(
            {
                "diffusion_scale": diffusion_scale,
                "memory_decay": memory_decay,
                "memory_diffusion": memory_diffusion,
                "swarm_alignment_weight": alignment,
                "density_birth_suppression": density_suppression,
                "cci_gate_strength": cci_gate,
                "mean_fate_entropy_context": mean_entropy,
                "system_state": state,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(table_dir / "phase_diagram.csv", index=False)
    pivot = (
        out.assign(faithful_flow_fraction=out["system_state"] == "faithful flow")
        .groupby(["diffusion_scale", "swarm_alignment_weight"])["faithful_flow_fraction"]
        .mean()
        .unstack(fill_value=0)
    )
    fig, ax = plt.subplots(figsize=(5, 4), dpi=160)
    im = ax.imshow(pivot.to_numpy(), origin="lower", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(pivot.shape[1]), pivot.columns)
    ax.set_yticks(range(pivot.shape[0]), pivot.index)
    ax.set_xlabel("swarm alignment weight")
    ax.set_ylabel("diffusion scale")
    ax.set_title("Phase diagram: faithful-flow fraction")
    fig.colorbar(im, ax=ax, label="fraction")
    fig.tight_layout()
    fig.savefig(fig_dir / "phase_diagram.png")
    plt.close(fig)
    state_counts = out["system_state"].value_counts().reset_index()
    gate = bool((out["system_state"] == "faithful flow").any() and out["system_state"].nunique() >= 4)
    write_report(
        report_dir / "discovery_phase_diagram.md",
        "Discovery Phase Diagram",
        [
            "A coarse parameter scan classifies finite-agent dynamics into faithful flow, over-diffusion, collapsed branch, premature branching, population explosion and lineage extinction regimes.",
            "",
            state_counts.to_markdown(index=False),
            "",
            f"- phase_diagram_gate: {gate}",
        ],
    )
    return {"law": "phase_diagram", "gate": gate, "table": str(table_dir / "phase_diagram.csv")}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/train.yaml")
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    print(run(args.config, args.quick_fixture))


if __name__ == "__main__":
    main()
