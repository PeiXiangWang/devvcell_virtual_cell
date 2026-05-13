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
    linear_effects,
    load_teacher,
    output_dirs,
    permutation_p_value,
    seed_list,
    seed_stability,
    seedwise_feature_frame,
    standardized_coef,
    write_report,
)


LAW = "diffusion"


def _permutation_null(frame: pd.DataFrame, repeats: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    base = frame.copy()
    for rep in range(repeats):
        shuffled = base.copy()
        shuffled["ot_transition_entropy"] = rng.permutation(shuffled["ot_transition_entropy"].to_numpy())
        coef = standardized_coef(
            shuffled,
            "learned_sigma",
            "ot_transition_entropy",
            ["local_density", "fate_probability_max", "fate_entropy", "cell_cycle_score", "cci_signal"],
        )
        rows.append({"permutation": rep, "entropy_coef_null": coef})
    return pd.DataFrame(rows)


def run(config: str = "configs/train.yaml", quick_fixture: bool = False):
    train_cfg, model_cfg, discovery_cfg = configure_paths(config, quick_fixture)
    table_dir, report_dir, fig_dir = output_dirs(train_cfg, discovery_cfg)
    adata = load_teacher(model_cfg)
    seeds = seed_list(model_cfg, quick_fixture)
    frame = seedwise_feature_frame(adata, model_cfg, seeds)
    predictors = ["ot_transition_entropy", "local_density", "fate_probability_max", "fate_entropy", "cell_cycle_score", "cci_signal"]
    standard = linear_effects(frame, "learned_sigma", predictors)
    standard["analysis_mode"] = "standard"
    no_entropy = linear_effects(frame, "learned_sigma", [p for p in predictors if p != "ot_transition_entropy"])
    no_entropy["analysis_mode"] = "no_entropy_input_proxy"
    empirical = linear_effects(frame, "learned_displacement", predictors)
    empirical["analysis_mode"] = "empirical_displacement_proxy"
    stats = pd.concat([standard, no_entropy, empirical], ignore_index=True)
    stats.to_csv(table_dir / "diffusion_law_regression.csv", index=False)

    seed_rows = []
    for seed, group in frame.groupby("seed", observed=False):
        effect = standardized_coef(group, "learned_sigma", "ot_transition_entropy", [p for p in predictors if p != "ot_transition_entropy"])
        no_entropy_r2 = linear_effects(group, "learned_sigma", [p for p in predictors if p != "ot_transition_entropy"])["r2"].max()
        empirical_effect = standardized_coef(group, "learned_displacement", "ot_transition_entropy", [p for p in predictors if p != "ot_transition_entropy"])
        seed_rows.append({"seed": int(seed), "entropy_effect": effect, "no_entropy_proxy_r2": float(no_entropy_r2), "empirical_displacement_entropy_effect": empirical_effect})
    seedwise = pd.DataFrame(seed_rows)
    seedwise.to_csv(table_dir / "diffusion_law_seedwise.csv", index=False)

    repeats = int(discovery_cfg.get("emergent_law", {}).get("permutation_repeats", 100))
    null = _permutation_null(frame, repeats, 101)
    null.to_csv(table_dir / "diffusion_law_permutation_null.csv", index=False)

    effects = seedwise["entropy_effect"].to_numpy(dtype=float)
    effect_mean, ci_low, ci_high = bootstrap_ci(effects, repeats=int(discovery_cfg.get("emergent_law", {}).get("bootstrap_repeats", 500)))
    p = permutation_p_value(effect_mean, null["entropy_coef_null"].to_numpy(dtype=float))
    q = bh_q_values([p])[0]
    min_seed = int(discovery_cfg.get("emergent_law", {}).get("quick_min_seed_count" if quick_fixture else "min_seed_count", 5))
    stability, sign_consistency = seed_stability(effects, min_seed)
    negative_control_pass = bool(q <= float(discovery_cfg.get("emergent_law", {}).get("max_permutation_q_for_acceptable", 0.10)))
    directly_encoded = True
    rollout_based = False
    tier = law_tier(effect_mean, q, negative_control_pass, stability, rollout_based, directly_encoded, discovery_cfg)
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
        str(table_dir / "diffusion_law_regression.csv"),
        str(report_dir / "discovery_diffusion_law.md"),
    )

    plot_frame = frame[frame["seed"] == seeds[0]]
    fig, ax = plt.subplots(figsize=(6, 4), dpi=160)
    sc = ax.scatter(plot_frame["ot_transition_entropy"], plot_frame["learned_sigma"], c=plot_frame["local_density"], s=7, cmap="viridis", alpha=0.7)
    ax.set_xlabel("OT transition entropy")
    ax.set_ylabel("learned sigma")
    ax.set_title("Diffusion law: uncertainty and density")
    fig.colorbar(sc, ax=ax, label="local density")
    fig.tight_layout()
    fig.savefig(fig_dir / "diffusion_entropy_density.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.5, 3.8), dpi=160)
    ax.axhline(0, color="black", lw=0.8)
    ax.errorbar([0], [effect_mean], yerr=[[effect_mean - ci_low], [ci_high - effect_mean]], fmt="o", color="#315b7d", capsize=4)
    ax.scatter(np.full(seedwise.shape[0], 0.08), seedwise["entropy_effect"], s=24, color="#6d9dc5", edgecolor="#1b2935", linewidth=0.3)
    ax.set_xlim(-0.5, 0.6)
    ax.set_xticks([0], ["entropy effect"])
    ax.set_ylabel("standardized coefficient")
    ax.set_title("Seed-wise diffusion effect")
    fig.tight_layout()
    fig.savefig(fig_dir / "diffusion_effect_seed_ci.png")
    plt.close(fig)

    top = stats.sort_values("abs_coef", ascending=False).head(8)
    write_report(
        report_dir / "discovery_diffusion_law.md",
        "Discovery Diffusion Law",
        [
            "Diffusion is hardened with standard regression, entropy-shuffle negative control, no-entropy proxy analysis and seed-wise coefficient stability.",
            "",
            "The current training objective directly includes a sigma-to-entropy calibration term. Therefore entropy-associated learned sigma is classified as encoded control-law recovery unless independently supported by rollout displacement.",
            "",
            "## Tier",
            "",
            f"- tier: {tier}",
            f"- effect_size: {effect_mean:.6g}",
            f"- 95% CI: [{ci_low:.6g}, {ci_high:.6g}]",
            f"- permutation_p: {p:.6g}",
            f"- permutation_q: {q:.6g}",
            f"- seed_stability_pass: {stability} (sign consistency={sign_consistency:.3f})",
            f"- negative_control_pass: {negative_control_pass}",
            f"- directly_supervised_or_encoded: {directly_encoded}",
            "",
            "## Top Regression Terms",
            "",
            top.to_markdown(index=False) if not top.empty else "No stable regression fit.",
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
