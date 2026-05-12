from __future__ import annotations

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from src.discovery.common import cell_feature_frame, configure_paths, linear_effects, load_teacher, output_dirs, write_report


def run(config: str = "configs/train.yaml", quick_fixture: bool = False):
    train_cfg, model_cfg = configure_paths(config, quick_fixture)
    table_dir, report_dir, fig_dir = output_dirs(train_cfg)
    adata = load_teacher(model_cfg)
    frame = cell_feature_frame(adata, model_cfg)
    predictors = ["ot_transition_entropy", "local_density", "fate_probability_max", "cell_cycle_score", "cci_signal"]
    stats = linear_effects(frame, "learned_sigma", predictors)
    stats.to_csv(table_dir / "diffusion_law_regression.csv", index=False)
    fig, ax = plt.subplots(figsize=(6, 4), dpi=160)
    sc = ax.scatter(frame["ot_transition_entropy"], frame["learned_sigma"], c=frame["local_density"], s=7, cmap="viridis", alpha=0.7)
    ax.set_xlabel("OT transition entropy")
    ax.set_ylabel("learned sigma")
    ax.set_title("Diffusion law: uncertainty and density")
    fig.colorbar(sc, ax=ax, label="local density")
    fig.tight_layout()
    fig.savefig(fig_dir / "diffusion_entropy_density.png")
    plt.close(fig)
    top = stats.sort_values("abs_coef", ascending=False).head(3)
    gate = bool((stats["abs_coef"] > 1e-4).sum() >= 2 and stats["r2"].max() > 0.01) if not stats.empty else False
    write_report(
        report_dir / "discovery_diffusion_law.md",
        "Discovery Diffusion Law",
        [
            "The analysis regresses learned diffusion scale against OT transition entropy, local density, fate commitment, cell-cycle score and CCI signal.",
            "",
            "## Top Effects",
            "",
            top.to_markdown(index=False) if not top.empty else "No stable regression fit.",
            "",
            f"- diffusion_law_gate: {gate}",
            "- Interpretation is mechanistic only if the signal is stable across seeds and not driven by teacher fallback artifacts.",
        ],
    )
    return {"law": "diffusion", "gate": gate, "table": str(table_dir / "diffusion_law_regression.csv")}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/train.yaml")
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    print(run(args.config, args.quick_fixture))


if __name__ == "__main__":
    main()

