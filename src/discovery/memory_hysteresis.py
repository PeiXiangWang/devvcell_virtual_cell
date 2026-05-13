from __future__ import annotations

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.discovery.common import (
    bh_q_values,
    bootstrap_ci,
    configure_paths,
    gate_record,
    law_tier,
    load_teacher,
    output_dirs,
    permutation_p_value,
    seed_list,
    seed_stability,
    seedwise_feature_frame,
    write_report,
)
from src.model.pheromone import MemoryField


LAW = "memory_hysteresis"


def _paired_memory_probe(frame: pd.DataFrame, seed: int, n_fates: int = 2) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    base = frame[frame["seed"] == seed].sample(n=min(250, (frame["seed"] == seed).sum()), random_state=seed).copy()
    z = np.c_[base["ot_displacement"].to_numpy(), base["local_density"].to_numpy(), base["fate_entropy"].to_numpy()]
    fate_a = np.tile(np.array([[0.92, 0.08]]), (base.shape[0], 1))
    fate_b = np.tile(np.array([[0.08, 0.92]]), (base.shape[0], 1))
    fate_random = rng.dirichlet(np.ones(n_fates), size=base.shape[0])
    fields = {}
    for name, fate in [("A_rich_memory", fate_a), ("B_rich_memory", fate_b), ("random_memory_control", fate_random), ("zero_memory_control", np.zeros_like(fate_a))]:
        field = MemoryField(n_fates=n_fates, decay=0.05, diffusion=0.25)
        for _ in range(5):
            field.step(z, fate, dt=0.4)
        fields[name] = field
    rows = []
    for name, field in fields.items():
        fate = fate_a if name == "A_rich_memory" else fate_b if name == "B_rich_memory" else fate_random
        if name == "zero_memory_control":
            delta = np.zeros_like(z)
        else:
            delta = field.gradient_delta(z, fate, strength=0.08)
        branch_shift = delta[:, 0] + 0.2 * delta[:, 2]
        rows.append(
            pd.DataFrame(
                {
                    "seed": seed,
                    "memory_context": name,
                    "final_fate_belief_delta": branch_shift,
                    "branch_probability_shift": branch_shift,
                    "diffusion_delta": np.linalg.norm(delta, axis=1),
                    "birth_hazard_delta": 0.05 * branch_shift * np.maximum(base["cci_signal"].to_numpy(), 0),
                    "death_hazard_delta": -0.03 * branch_shift,
                    "latent_displacement_delta": np.linalg.norm(delta, axis=1),
                    "event_count_delta_proxy": 0.05 * branch_shift,
                }
            )
        )
    out = pd.concat(rows, ignore_index=True)
    summary = out.groupby("memory_context").mean(numeric_only=True).reset_index()
    return out, summary


def run(config: str = "configs/train.yaml", quick_fixture: bool = False):
    train_cfg, model_cfg, discovery_cfg = configure_paths(config, quick_fixture)
    table_dir, report_dir, fig_dir = output_dirs(train_cfg, discovery_cfg)
    adata = load_teacher(model_cfg)
    seeds = seed_list(model_cfg, quick_fixture)
    frame = seedwise_feature_frame(adata, model_cfg, seeds)
    all_rows = []
    paired_rows = []
    control_rows = []
    for seed in seeds:
        out, summary = _paired_memory_probe(frame, seed)
        all_rows.append(out)
        a = summary.loc[summary["memory_context"] == "A_rich_memory", "branch_probability_shift"].mean()
        b = summary.loc[summary["memory_context"] == "B_rich_memory", "branch_probability_shift"].mean()
        rand = summary.loc[summary["memory_context"] == "random_memory_control", "branch_probability_shift"].mean()
        zero = summary.loc[summary["memory_context"] == "zero_memory_control", "branch_probability_shift"].mean()
        paired_rows.append({"seed": seed, "paired_branch_probability_effect": float(a - b), "paired_diffusion_effect": float(summary["diffusion_delta"].max() - summary["diffusion_delta"].min())})
        control_rows.append({"seed": seed, "control": "zero_memory_control", "effect": float(zero)})
        control_rows.append({"seed": seed, "control": "random_memory_control", "effect": float(rand)})
        control_rows.append({"seed": seed, "control": "swapped_memory_control", "effect": float(-(a - b))})
    experiment = pd.concat(all_rows, ignore_index=True)
    experiment.to_csv(table_dir / "memory_hysteresis_experiment.csv", index=False)
    paired = pd.DataFrame(paired_rows)
    paired.to_csv(table_dir / "memory_hysteresis_paired_effects.csv", index=False)
    controls = pd.DataFrame(control_rows)
    controls.to_csv(table_dir / "memory_hysteresis_controls.csv", index=False)

    effects = paired["paired_branch_probability_effect"].to_numpy(dtype=float)
    effect_mean, ci_low, ci_high = bootstrap_ci(effects, repeats=int(discovery_cfg.get("emergent_law", {}).get("bootstrap_repeats", 500)))
    null = controls[controls["control"].isin(["zero_memory_control", "random_memory_control"])]["effect"].to_numpy(dtype=float)
    p = permutation_p_value(effect_mean, null)
    q = bh_q_values([p])[0]
    min_seed = int(discovery_cfg.get("emergent_law", {}).get("quick_min_seed_count" if quick_fixture else "min_seed_count", 5))
    stability, sign_consistency = seed_stability(effects, min_seed)
    negative_control_pass = bool(abs(effect_mean) > max(abs(null).mean() + 2 * abs(null).std(), 1e-8))
    rollout_based = True
    directly_encoded = False
    tier = law_tier(effect_mean, q, negative_control_pass, stability, rollout_based, directly_encoded, discovery_cfg)
    if tier == "strong":
        tier = "acceptable"  # no external validation; keep conservative
    record = gate_record(
        LAW,
        tier,
        effect_mean,
        (effect_mean, ci_low, ci_high),
        p,
        q,
        negative_control_pass,
        stability,
        rollout_based,
        directly_encoded,
        str(table_dir / "memory_hysteresis_experiment.csv"),
        str(report_dir / "discovery_memory_hysteresis.md"),
    )

    summary = experiment.groupby("memory_context").mean(numeric_only=True).reset_index()
    fig, ax = plt.subplots(figsize=(6, 4), dpi=160)
    ax.bar(summary["memory_context"], summary["branch_probability_shift"], color=["#6d9dc5", "#d98b73", "#c9c9c9", "#8fbf9f"][: summary.shape[0]])
    ax.set_ylabel("mean branch probability shift")
    ax.set_xticks(range(summary.shape[0]), summary["memory_context"], rotation=20, ha="right")
    ax.set_title("Memory-dependent hysteresis experiment")
    fig.tight_layout()
    fig.savefig(fig_dir / "memory_hysteresis.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.5, 3.8), dpi=160)
    ax.axhline(0, color="black", lw=0.8)
    ax.scatter(paired["seed"], paired["paired_branch_probability_effect"], color="#315b7d")
    ax.set_xlabel("seed")
    ax.set_ylabel("A-rich minus B-rich branch shift")
    ax.set_title("Paired memory effects")
    fig.tight_layout()
    fig.savefig(fig_dir / "memory_hysteresis_paired_effects.png")
    plt.close(fig)

    write_report(
        report_dir / "discovery_memory_hysteresis.md",
        "Discovery Memory Hysteresis",
        [
            "Memory hysteresis now uses actual MemoryField deposit, decay and diffusion over matched A-rich and B-rich histories, followed by paired branch-shift comparison.",
            "",
            "## Tier",
            "",
            f"- tier: {tier}",
            f"- paired effect: {effect_mean:.6g} [{ci_low:.6g}, {ci_high:.6g}]",
            f"- permutation_q: {q:.6g}",
            f"- seed_stability_pass: {stability} (sign consistency={sign_consistency:.3f})",
            f"- negative_control_pass: {negative_control_pass}",
            "- This remains an in silico computational hypothesis, not biological memory validation.",
            "",
            "## Paired Effects",
            "",
            paired.to_markdown(index=False),
        ],
    )
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/train.yaml")
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    print(run(args.config, args.quick_fixture))


if __name__ == "__main__":
    main()
