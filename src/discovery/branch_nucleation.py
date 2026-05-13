from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

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


LAW = "branch_nucleation"


def _alignment(z: np.ndarray, v: np.ndarray, k: int = 12) -> np.ndarray:
    if z.shape[0] <= 2:
        return np.ones(z.shape[0])
    k = min(k, z.shape[0] - 1)
    neigh = NearestNeighbors(n_neighbors=k + 1).fit(z).kneighbors(z, return_distance=False)[:, 1:]
    vn = v / np.maximum(np.linalg.norm(v, axis=1, keepdims=True), 1e-8)
    return np.array([(vn[i] * vn[neigh[i]].mean(axis=0)).sum() for i in range(z.shape[0])])


def _static_order_parameters(adata, model_cfg: dict) -> pd.DataFrame:
    frame = seedwise_feature_frame(adata, model_cfg, [7]).drop(columns=["seed"])
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
                "seed": 7,
                "variant": "static_teacher",
                "step": int(time),
                "time": float(time),
                "local_velocity_alignment_A": float(group["alignment_cell"].mean()),
                "branch_cohesion_C": float(1.0 - group.groupby("lineage")["local_density"].std().fillna(0).mean()),
                "lineage_separation_S": float(group.groupby("lineage")["ot_displacement"].mean().std()),
                "fate_entropy_H": float(group["fate_entropy"].mean()),
                "branch_imbalance_B": imbalance,
                "local_density_mean": float(group["local_density"].mean()),
                "local_density_var": float(group["local_density"].var()),
                "n_agents": int(group.shape[0]),
                "per_lineage_counts": "{}",
            }
        )
    return pd.DataFrame(rows)


def _event_windows(order: pd.DataFrame) -> pd.DataFrame:
    rows = []
    metrics = ["local_velocity_alignment_A", "branch_cohesion_C", "lineage_separation_S", "fate_entropy_H", "branch_imbalance_B"]
    for (seed, variant), group in order.groupby(["seed", "variant"], observed=False):
        group = group.sort_values("step").reset_index(drop=True)
        if group.shape[0] < 3:
            continue
        sep_delta = group["lineage_separation_S"].diff().fillna(0)
        entropy_delta = -group["fate_entropy_H"].diff().fillna(0)
        imbalance_delta = group["branch_imbalance_B"].diff().fillna(0)
        score = sep_delta.rank(pct=True) + entropy_delta.rank(pct=True) + imbalance_delta.rank(pct=True)
        event_pos = int(score.idxmax())
        pre = group.iloc[max(0, event_pos - 2) : event_pos]
        post = group.iloc[event_pos + 1 : min(group.shape[0], event_pos + 3)]
        if pre.empty or post.empty:
            continue
        row = {"seed": int(seed), "variant": variant, "branch_event_step": int(group.loc[event_pos, "step"]), "event_definition": "max_ranked_separation_entropy_imbalance_change"}
        for metric in metrics:
            row[f"{metric}_pre_mean"] = float(pre[metric].mean())
            row[f"{metric}_post_mean"] = float(post[metric].mean())
            row[f"{metric}_effect"] = float(post[metric].mean() - pre[metric].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def run(config: str = "configs/train.yaml", quick_fixture: bool = False):
    train_cfg, model_cfg, discovery_cfg = configure_paths(config, quick_fixture)
    table_dir, report_dir, fig_dir = output_dirs(train_cfg, discovery_cfg)
    adata = load_teacher(model_cfg)
    order_path = Path(model_cfg.get("order_log_path", table_dir / "rollout_order_parameters.csv"))
    if order_path.exists():
        order = pd.read_csv(order_path)
        rollout_based = bool("variant" in order and not order["variant"].astype(str).eq("static_teacher").all())
    else:
        order = _static_order_parameters(adata, model_cfg)
        rollout_based = False
    if "variant" in order:
        preferred = order[order["variant"].isin(["M9_full_memory", "M8_ot_swarm_birth_death_diffusion_cci", "M7_ot_swarm_birth_death_diffusion"])]
        if not preferred.empty:
            order = preferred.copy()
    order.to_csv(table_dir / "rollout_order_parameters.csv", index=False)
    stage_summary = order.groupby(["seed", "step"], observed=False).mean(numeric_only=True).reset_index()
    stage_summary.to_csv(table_dir / "swarm_order_parameters.csv", index=False)

    windows = _event_windows(order)
    windows.to_csv(table_dir / "branch_nucleation_event_windows.csv", index=False)
    effect_col = "lineage_separation_S_effect"
    effects = windows[effect_col].to_numpy(dtype=float) if effect_col in windows else np.array([])
    effect_mean, ci_low, ci_high = bootstrap_ci(effects, repeats=int(discovery_cfg.get("emergent_law", {}).get("bootstrap_repeats", 500)))

    rng = np.random.default_rng(23)
    null_rows = []
    for rep in range(int(discovery_cfg.get("emergent_law", {}).get("permutation_repeats", 100))):
        shuffled = order.copy()
        shuffled["step"] = shuffled.groupby(["seed", "variant"], observed=False)["step"].transform(lambda x: rng.permutation(x.to_numpy()))
        sh_win = _event_windows(shuffled)
        null_effect = float(sh_win[effect_col].mean()) if effect_col in sh_win and not sh_win.empty else 0.0
        null_rows.append({"control": "shuffled_temporal_order", "permutation": rep, "effect_null": null_effect})
    controls = pd.DataFrame(null_rows)
    controls.to_csv(table_dir / "branch_nucleation_negative_controls.csv", index=False)
    p = permutation_p_value(effect_mean, controls["effect_null"].to_numpy(dtype=float))
    q = bh_q_values([p])[0]
    min_seed = int(discovery_cfg.get("emergent_law", {}).get("quick_min_seed_count" if quick_fixture else "min_seed_count", 5))
    stability, sign_consistency = seed_stability(effects, min_seed)
    negative_control_pass = bool(q <= float(discovery_cfg.get("emergent_law", {}).get("max_permutation_q_for_acceptable", 0.10)))
    directly_encoded = False
    tier = law_tier(effect_mean, q, negative_control_pass, stability, rollout_based, directly_encoded, discovery_cfg)
    if not rollout_based and tier in {"acceptable", "strong"}:
        tier = "weak"
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
        str(table_dir / "swarm_order_parameters.csv"),
        str(report_dir / "discovery_branch_nucleation.md"),
    )

    fig, ax = plt.subplots(figsize=(7, 4), dpi=160)
    plot = stage_summary.groupby("step", observed=False).mean(numeric_only=True).reset_index()
    for col in ["local_velocity_alignment_A", "branch_cohesion_C", "lineage_separation_S", "fate_entropy_H", "branch_imbalance_B"]:
        if col not in plot:
            continue
        vals = plot[col].to_numpy(dtype=float)
        vals = (vals - np.nanmin(vals)) / max(np.nanmax(vals) - np.nanmin(vals), 1e-8)
        ax.plot(plot["step"], vals, marker="o", label=col)
    ax.set_xlabel("rollout step")
    ax.set_ylabel("normalized order parameter")
    ax.set_title("Branch nucleation order parameters")
    ax.legend(fontsize=6, frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "branch_nucleation_order_parameters.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4), dpi=160)
    effect_cols = [c for c in windows.columns if c.endswith("_effect")]
    if effect_cols and not windows.empty:
        means = windows[effect_cols].mean().sort_values()
        ax.barh(range(len(means)), means, color="#7da7c7")
        ax.set_yticks(range(len(means)), [c.replace("_effect", "") for c in means.index], fontsize=7)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("post - pre effect")
    ax.set_title("Branch event window effects")
    fig.tight_layout()
    fig.savefig(fig_dir / "branch_event_window_effects.png")
    plt.close(fig)

    write_report(
        report_dir / "discovery_branch_nucleation.md",
        "Discovery Branch Nucleation",
        [
            "Branch nucleation is hardened by rollout order-parameter traces, event-window detection and shuffled temporal-order controls.",
            "",
            "## Tier",
            "",
            f"- tier: {tier}",
            f"- rollout_based: {rollout_based}",
            f"- lineage_separation_event_effect: {effect_mean:.6g} [{ci_low:.6g}, {ci_high:.6g}]",
            f"- permutation_q: {q:.6g}",
            f"- seed_stability_pass: {stability} (sign consistency={sign_consistency:.3f})",
            f"- negative_control_pass: {negative_control_pass}",
            "",
            "## Event Windows",
            "",
            windows.head(20).to_markdown(index=False) if not windows.empty else "No event window could be detected.",
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
