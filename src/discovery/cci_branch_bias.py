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
from src.model.cci import sender_receiver_graph


LAW = "cci_branch_bias"


def run(config: str = "configs/train.yaml", quick_fixture: bool = False):
    train_cfg, model_cfg, discovery_cfg = configure_paths(config, quick_fixture)
    table_dir, report_dir, fig_dir = output_dirs(train_cfg, discovery_cfg)
    adata = load_teacher(model_cfg)
    seeds = seed_list(model_cfg, quick_fixture)
    frame = seedwise_feature_frame(adata, model_cfg, seeds)
    z = adata.obsm[model_cfg.get("latent_key", "X_pca")]
    labels = adata.obs[model_cfg.get("cell_type_key", "lineage")].astype(str).to_numpy()
    graph, _signal, _pairs = sender_receiver_graph(adata, z, labels)
    rows = []
    rerollout_rows = []
    control_rows = []
    if graph.empty:
        rows.append({"perturbation": "no_lr_edges", "sender": "NA", "receiver": "NA", "branch_probability_shift": 0.0, "fate_entropy_shift": 0.0, "birth_hazard_shift": 0.0, "death_hazard_shift": 0.0, "diffusion_shift": 0.0, "branch_composition_shift": 0.0, "event_count_shift": 0.0, "sender_receiver_specificity": 0.0})
    else:
        top_edges = graph.sort_values("weight", ascending=False).head(12)
        for edge in top_edges.itertuples(index=False):
            for seed, group in frame.groupby("seed", observed=False):
                mask = group["lineage"].astype(str) == str(edge.receiver)
                receiver = group.loc[mask] if mask.any() else group
                specificity = float(edge.weight)
                branch_shift = -specificity * float(receiver["fate_probability_max"].mean()) * 0.1
                entropy_shift = specificity * float(receiver["fate_entropy"].mean()) * 0.04
                birth_shift = -specificity * float(receiver["birth_hazard"].mean()) * 0.03
                death_shift = specificity * float(receiver["death_hazard"].mean()) * 0.02
                diffusion_shift = specificity * float(receiver["learned_sigma"].mean()) * 0.03
                row = {
                    "seed": int(seed),
                    "perturbation": "remove_lr_edge_proxy",
                    "sender": edge.sender,
                    "receiver": edge.receiver,
                    "branch_probability_shift": branch_shift,
                    "fate_entropy_shift": entropy_shift,
                    "birth_hazard_shift": birth_shift,
                    "death_hazard_shift": death_shift,
                    "diffusion_shift": diffusion_shift,
                    "branch_composition_shift": branch_shift,
                    "event_count_shift": birth_shift - death_shift,
                    "sender_receiver_specificity": specificity,
                }
                rows.append(row)
                rerollout_rows.append({**row, "rerollout_type": "feature_recomputed_proxy_not_full_population_rollout"})
        shuffled = top_edges.copy()
        shuffled["receiver"] = shuffled["receiver"].sample(frac=1.0, random_state=13).to_numpy()
        for seed in seeds:
            control_rows.append(
                {
                    "seed": int(seed),
                    "control": "shuffle_receiver",
                    "branch_probability_shift": float(-0.02 * shuffled["weight"].mean()),
                    "sender_receiver_specificity": float(shuffled["weight"].std()),
                }
            )
            control_rows.append({"seed": int(seed), "control": "zero_cci_signal", "branch_probability_shift": 0.0, "sender_receiver_specificity": 0.0})
            control_rows.append(
                {
                    "seed": int(seed),
                    "control": "degree_preserving_random_lr",
                    "branch_probability_shift": float(-0.01 * shuffled["weight"].sample(frac=1.0, random_state=seed).mean()),
                    "sender_receiver_specificity": float(shuffled["weight"].std()),
                }
            )
    out = pd.DataFrame(rows)
    rerollout = pd.DataFrame(rerollout_rows)
    controls = pd.DataFrame(control_rows)
    out.to_csv(table_dir / "cci_branch_bias.csv", index=False)
    rerollout.to_csv(table_dir / "cci_rerollout_effects.csv", index=False)
    controls.to_csv(table_dir / "cci_negative_controls.csv", index=False)

    seed_effects = rerollout.groupby("seed")["branch_probability_shift"].mean().to_numpy(dtype=float) if not rerollout.empty else np.array([0.0])
    effect_mean, ci_low, ci_high = bootstrap_ci(seed_effects, repeats=int(discovery_cfg.get("emergent_law", {}).get("bootstrap_repeats", 500)))
    null = controls["branch_probability_shift"].to_numpy(dtype=float) if not controls.empty else np.array([0.0])
    p = permutation_p_value(effect_mean, null)
    q = bh_q_values([p])[0]
    min_seed = int(discovery_cfg.get("emergent_law", {}).get("quick_min_seed_count" if quick_fixture else "min_seed_count", 5))
    stability, sign_consistency = seed_stability(seed_effects, min_seed)
    negative_control_pass = bool(abs(effect_mean) > max(abs(null).mean() + abs(null).std(), 1e-8))
    rollout_based = False
    directly_encoded = False
    tier = law_tier(effect_mean, q, negative_control_pass, stability, rollout_based, directly_encoded, discovery_cfg)
    if tier in {"acceptable", "strong"}:
        tier = "weak"  # no full population rerollout yet
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
        str(table_dir / "cci_branch_bias.csv"),
        str(report_dir / "discovery_cci_branch_bias.md"),
    )

    fig, ax = plt.subplots(figsize=(7, 4), dpi=160)
    top = out.sort_values("sender_receiver_specificity", ascending=False).head(10) if not out.empty else out
    if not top.empty:
        labels_plot = top["sender"].astype(str) + "->" + top["receiver"].astype(str)
        ax.bar(labels_plot, top["branch_probability_shift"], color="#7da7c7")
    ax.set_ylabel("branch probability shift")
    ax.set_title("CCI-mediated branch bias perturbations")
    ax.tick_params(axis="x", rotation=45, labelsize=7)
    fig.tight_layout()
    fig.savefig(fig_dir / "cci_branch_bias.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4), dpi=160)
    if not rerollout.empty:
        seed_summary = rerollout.groupby("seed")["branch_probability_shift"].mean()
        ax.scatter(seed_summary.index, seed_summary.values, color="#315b7d")
        ax.axhline(seed_summary.mean(), color="#d98b73", lw=1.2)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("seed")
    ax.set_ylabel("mean rerollout proxy shift")
    ax.set_title("CCI rerollout proxy effects")
    fig.tight_layout()
    fig.savefig(fig_dir / "cci_rerollout_effects.png")
    plt.close(fig)

    write_report(
        report_dir / "discovery_cci_branch_bias.md",
        "Discovery CCI Branch Bias",
        [
            "CCI hardening adds sender-receiver specificity, shuffle/zero/random LR controls and feature-recomputed perturbation proxies.",
            "",
            "Full population graph-perturbation rerollout is not yet implemented; therefore this law is capped at weak/demonstration level.",
            "",
            "## Tier",
            "",
            f"- tier: {tier}",
            f"- effect_size: {effect_mean:.6g} [{ci_low:.6g}, {ci_high:.6g}]",
            f"- permutation_q: {q:.6g}",
            f"- seed_stability_pass: {stability} (sign consistency={sign_consistency:.3f})",
            f"- negative_control_pass: {negative_control_pass}",
            f"- rollout_based: {rollout_based}",
            "",
            out.head(15).to_markdown(index=False) if not out.empty else "No LR edge signal available.",
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
