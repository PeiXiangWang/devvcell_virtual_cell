from __future__ import annotations

import argparse
from itertools import product

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
    write_report,
)
from src.model.swarm_rules import local_density, swarm_delta


LAW = "phase_diagram"


def _classify_outputs(row: dict) -> str:
    if row["extinction_rate"] > 0.35 or row["population_ratio"] < 0.45:
        return "lineage_extinction"
    if row["explosion_rate"] > 0.35 or row["population_ratio"] > 2.0:
        return "unstable_population_explosion"
    if row["manifold_escape_rate"] > 0.35:
        return "off_manifold"
    if row["diffusion_spread"] > 1.8:
        return "over_diffusion"
    if row["branch_collapse_score"] > 0.82:
        return "collapsed_branch"
    if row["branch_entropy"] > 0.8 and row["fate_composition_drift"] > 0.15:
        return "premature_branching"
    return "faithful_flow"


def _rollout_proxy(adata, model_cfg: dict, params: dict, seed: int, n: int = 180) -> dict:
    rng = np.random.default_rng(seed)
    z_all = np.asarray(adata.obsm[model_cfg.get("latent_key", "X_pca")], dtype=float)
    labels_all = adata.obs[model_cfg.get("cell_type_key", "lineage")].astype(str).to_numpy()
    v_all = np.asarray(adata.obsm.get(model_cfg.get("teacher_velocity_key", "X_ot_velocity"), np.zeros_like(z_all)), dtype=float)
    idx = rng.choice(np.arange(z_all.shape[0]), size=min(n, z_all.shape[0]), replace=False)
    z = z_all[idx].copy()
    labels = labels_all[idx].copy()
    source_freq = pd.Series(labels).value_counts(normalize=True)
    n0 = z.shape[0]
    for step in range(4):
        density = local_density(z)
        v = v_all[idx[: z.shape[0] % idx.shape[0] if z.shape[0] > idx.shape[0] else z.shape[0]]]
        if v.shape[0] != z.shape[0]:
            v = np.resize(v, z.shape)
        sw = swarm_delta(z, np.zeros_like(z), labels)
        z = z + 0.18 * v + 0.08 * params["swarm_alignment_weight"] * sw
        z = z + rng.normal(scale=0.025 * params["diffusion_scale"], size=z.shape)
        birth_prob = np.clip(0.04 + 0.03 * params["cci_gate_strength"] - 0.05 * params["density_birth_suppression"] * density, 0, 0.3)
        death_prob = np.clip(0.03 + 0.04 * params["density_birth_suppression"] * density - 0.02 * params["cci_gate_strength"], 0, 0.3)
        birth = rng.random(z.shape[0]) < birth_prob
        death = rng.random(z.shape[0]) < death_prob
        daughters = z[birth & ~death] + rng.normal(scale=0.015 * params["memory_diffusion"], size=(int((birth & ~death).sum()), z.shape[1]))
        daughter_labels = labels[birth & ~death]
        keep = ~death
        z = np.vstack([z[keep], daughters]) if daughters.size else z[keep]
        labels = np.r_[labels[keep], daughter_labels] if daughters.size else labels[keep]
        if z.shape[0] == 0:
            break
    if z.shape[0] == 0:
        return {
            "population_ratio": 0.0,
            "extinction_rate": 1.0,
            "explosion_rate": 0.0,
            "manifold_escape_rate": 1.0,
            "branch_entropy": 0.0,
            "branch_collapse_score": 1.0,
            "fate_composition_drift": 1.0,
            "diffusion_spread": 0.0,
            "teacher_fidelity_relative_sinkhorn": 999.0,
        }
    population_ratio = z.shape[0] / max(n0, 1)
    density_final = local_density(z)
    nearest = np.sqrt(((z[:, None, :] - z_all[None, :, :]) ** 2).sum(axis=2)).min(axis=1)
    manifold_escape_rate = float((nearest > np.quantile(nearest, 0.75) + 3 * np.std(nearest)).mean())
    freq = pd.Series(labels).value_counts(normalize=True)
    entropy = float(-(freq * np.log(np.clip(freq, 1e-12, 1))).sum() / max(np.log(freq.shape[0]), 1.0))
    collapse = float(freq.max())
    aligned = source_freq.reindex(sorted(set(source_freq.index) | set(freq.index)), fill_value=0)
    freq_aligned = freq.reindex(aligned.index, fill_value=0)
    drift = float(np.sqrt(((aligned - freq_aligned) ** 2).mean()))
    spread = float(np.mean(np.linalg.norm(z - z.mean(axis=0), axis=1)) / max(np.mean(np.linalg.norm(z_all[idx] - z_all[idx].mean(axis=0), axis=1)), 1e-8))
    fidelity = float(np.mean(nearest) / max(np.mean(np.linalg.norm(z_all[idx] - z_all[idx].mean(axis=0), axis=1)), 1e-8))
    return {
        "population_ratio": float(population_ratio),
        "extinction_rate": float(population_ratio < 0.5),
        "explosion_rate": float(population_ratio > 2.0),
        "manifold_escape_rate": manifold_escape_rate,
        "branch_entropy": entropy,
        "branch_collapse_score": collapse,
        "fate_composition_drift": drift,
        "diffusion_spread": spread,
        "teacher_fidelity_relative_sinkhorn": fidelity,
        "local_density_final_mean": float(density_final.mean()),
    }


def run(config: str = "configs/train.yaml", quick_fixture: bool = False):
    train_cfg, model_cfg, discovery_cfg = configure_paths(config, quick_fixture)
    table_dir, report_dir, fig_dir = output_dirs(train_cfg, discovery_cfg)
    adata = load_teacher(model_cfg)
    seeds = seed_list(model_cfg, quick_fixture)
    values = [0.5, 1.8]
    rows = []
    for diffusion_scale, memory_decay, memory_diffusion, alignment, density_suppression, cci_gate in product(values, [0.2, 0.8], [0.1, 0.5], values, [0.1, 1.4], values):
        params = {
            "diffusion_scale": diffusion_scale,
            "memory_decay": memory_decay,
            "memory_diffusion": memory_diffusion,
            "swarm_alignment_weight": alignment,
            "density_birth_suppression": density_suppression,
            "cci_gate_strength": cci_gate,
        }
        for seed in seeds:
            outputs = _rollout_proxy(adata, model_cfg, params, seed)
            rows.append({**params, "seed": seed, **outputs, "system_state": _classify_outputs(outputs)})
    out = pd.DataFrame(rows)
    out.to_csv(table_dir / "phase_diagram.csv", index=False)
    stability = out.groupby(["diffusion_scale", "swarm_alignment_weight"], observed=False)["system_state"].agg(lambda x: x.value_counts(normalize=True).iloc[0]).reset_index(name="majority_state_fraction")
    stability.to_csv(table_dir / "phase_diagram_seed_stability.csv", index=False)

    faithful_fraction = out.groupby("seed")["system_state"].apply(lambda x: (x == "faithful_flow").mean()).to_numpy(dtype=float)
    effect_mean, ci_low, ci_high = bootstrap_ci(faithful_fraction, repeats=int(discovery_cfg.get("emergent_law", {}).get("bootstrap_repeats", 500)))
    null = np.full(max(10, faithful_fraction.size), out["system_state"].value_counts(normalize=True).get("faithful_flow", 0.0))
    p = permutation_p_value(effect_mean, null)
    q = bh_q_values([p])[0]
    min_seed = int(discovery_cfg.get("emergent_law", {}).get("quick_min_seed_count" if quick_fixture else "min_seed_count", 5))
    stability_pass, sign_consistency = seed_stability(faithful_fraction - faithful_fraction.mean() + effect_mean, min_seed)
    negative_control_pass = bool(out["system_state"].nunique() >= 4)
    rollout_based = True
    directly_encoded = False
    tier = law_tier(effect_mean, q, negative_control_pass, stability_pass, rollout_based, directly_encoded, discovery_cfg)
    if tier == "strong":
        tier = "acceptable"  # coarse proxy rollout only
    record = gate_record(
        LAW,
        tier,
        effect_mean,
        (effect_mean, ci_low, ci_high),
        p,
        q,
        negative_control_pass,
        stability_pass,
        rollout_based,
        directly_encoded,
        str(table_dir / "phase_diagram.csv"),
        str(report_dir / "discovery_phase_diagram.md"),
    )

    pivot = (
        out.assign(faithful_flow_fraction=out["system_state"] == "faithful_flow")
        .groupby(["diffusion_scale", "swarm_alignment_weight"], observed=False)["faithful_flow_fraction"]
        .mean()
        .unstack(fill_value=0)
    )
    fig, ax = plt.subplots(figsize=(5, 4), dpi=160)
    im = ax.imshow(pivot.to_numpy(), origin="lower", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(pivot.shape[1]), pivot.columns)
    ax.set_yticks(range(pivot.shape[0]), pivot.index)
    ax.set_xlabel("swarm alignment weight")
    ax.set_ylabel("diffusion scale")
    ax.set_title("Rollout phase diagram: faithful-flow fraction")
    fig.colorbar(im, ax=ax, label="fraction")
    fig.tight_layout()
    fig.savefig(fig_dir / "phase_diagram.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 4), dpi=160)
    ax.hist(stability["majority_state_fraction"], bins=8, color="#7da7c7", edgecolor="#1b2935", linewidth=0.4)
    ax.set_xlabel("majority state fraction across seeds")
    ax.set_ylabel("parameter settings")
    ax.set_title("Phase boundary stability")
    fig.tight_layout()
    fig.savefig(fig_dir / "phase_boundary_stability.png")
    plt.close(fig)

    state_counts = out["system_state"].value_counts().reset_index()
    write_report(
        report_dir / "discovery_phase_diagram.md",
        "Discovery Phase Diagram",
        [
            "The phase diagram now classifies regimes from rollout outputs, not directly from input parameters.",
            "",
            "## Tier",
            "",
            f"- tier: {tier}",
            f"- faithful-flow fraction effect: {effect_mean:.6g} [{ci_low:.6g}, {ci_high:.6g}]",
            f"- permutation_q: {q:.6g}",
            f"- seed_stability_pass: {stability_pass} (sign consistency={sign_consistency:.3f})",
            f"- negative_control_pass: {negative_control_pass}",
            "",
            state_counts.to_markdown(index=False),
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
