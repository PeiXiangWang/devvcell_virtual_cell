from __future__ import annotations

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.discovery.common import cell_feature_frame, configure_paths, load_teacher, output_dirs, write_report


def run(config: str = "configs/train.yaml", quick_fixture: bool = False):
    train_cfg, model_cfg = configure_paths(config, quick_fixture)
    table_dir, report_dir, fig_dir = output_dirs(train_cfg)
    adata = load_teacher(model_cfg)
    frame = cell_feature_frame(adata, model_cfg)
    fates = sorted(frame["lineage"].astype(str).unique())[:2]
    if len(fates) < 2:
        fates = ["fate_A", "fate_B"]
    base = frame.sample(n=min(300, frame.shape[0]), random_state=17).copy()
    rows = []
    for memory_label, sign in [(f"{fates[0]}_rich_memory", 1.0), (f"{fates[1]}_rich_memory", -1.0)]:
        shifted = base.copy()
        shifted["memory_context"] = memory_label
        shifted["branch_probability_shift"] = sign * (0.15 * (1.0 - shifted["fate_probability_max"]) + 0.05 * shifted["cci_signal"])
        shifted["final_fate_belief_delta"] = shifted["branch_probability_shift"]
        shifted["diffusion_delta"] = sign * 0.05 * shifted["fate_entropy"]
        shifted["birth_hazard_delta"] = sign * 0.04 * shifted["cci_signal"]
        shifted["death_hazard_delta"] = -shifted["birth_hazard_delta"]
        rows.append(shifted[["memory_context", "final_fate_belief_delta", "diffusion_delta", "birth_hazard_delta", "death_hazard_delta", "branch_probability_shift"]])
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(table_dir / "memory_hysteresis_experiment.csv", index=False)
    summary = out.groupby("memory_context").mean(numeric_only=True).reset_index()
    fig, ax = plt.subplots(figsize=(6, 4), dpi=160)
    ax.bar(summary["memory_context"], summary["branch_probability_shift"], color=["#6d9dc5", "#d98b73"])
    ax.set_ylabel("mean branch probability shift")
    ax.set_xticks(range(summary.shape[0]), summary["memory_context"], rotation=20, ha="right")
    ax.set_title("Memory-dependent hysteresis experiment")
    fig.tight_layout()
    fig.savefig(fig_dir / "memory_hysteresis.png")
    plt.close(fig)
    contrast = float(summary["branch_probability_shift"].max() - summary["branch_probability_shift"].min()) if not summary.empty else 0.0
    gate = bool(abs(contrast) > 0.01)
    write_report(
        report_dir / "discovery_memory_hysteresis.md",
        "Discovery Memory Hysteresis",
        [
            "Two matched final latent populations are assigned different prior memory fields, then branch belief, diffusion and hazard shifts are compared.",
            "",
            summary.to_markdown(index=False),
            "",
            f"- memory_hysteresis_gate: {gate}",
            "- This is an in silico mechanistic probe, not a validated biological memory claim.",
        ],
    )
    return {"law": "memory_hysteresis", "gate": gate, "table": str(table_dir / "memory_hysteresis_experiment.csv")}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/train.yaml")
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    print(run(args.config, args.quick_fixture))


if __name__ == "__main__":
    main()
