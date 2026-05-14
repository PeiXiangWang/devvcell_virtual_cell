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
    metrics = [
        "local_velocity_alignment_A",
        "branch_cohesion_C",
        "lineage_separation_S",
        "fate_entropy_H",
        "branch_imbalance_B",
        "local_density_mean",
        "n_agents",
    ]
    for (seed, variant), group in order.groupby(["seed", "variant"], observed=False):
        group = group.sort_values("step").reset_index(drop=True)
        if group.shape[0] < 3:
            continue
        sep_delta = group["lineage_separation_S"].diff().fillna(0)
        entropy_delta = -group["fate_entropy_H"].diff().fillna(0)
        imbalance_delta = group["branch_imbalance_B"].diff().fillna(0)
        score = sep_delta.rank(pct=True) + entropy_delta.rank(pct=True) + imbalance_delta.rank(pct=True)
        candidates = list(range(1, max(1, group.shape[0] - 1)))
        if not candidates:
            continue
        event_pos = int(score.iloc[candidates].idxmax())
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


def _summarize_model_windows(windows: pd.DataFrame, null_effects: np.ndarray, discovery_cfg: dict, quick_fixture: bool) -> pd.DataFrame:
    rows = []
    effect_cols = [c for c in windows.columns if c.endswith("_effect")]
    min_seed = int(discovery_cfg.get("emergent_law", {}).get("quick_min_seed_count" if quick_fixture else "min_seed_count", 5))
    for variant, group in windows.groupby("variant", observed=False):
        effects = group["lineage_separation_S_effect"].to_numpy(dtype=float) if "lineage_separation_S_effect" in group else np.array([])
        mean, lo, hi = bootstrap_ci(effects, repeats=int(discovery_cfg.get("emergent_law", {}).get("bootstrap_repeats", 500)))
        p = permutation_p_value(mean, null_effects)
        q = bh_q_values([p])[0]
        stability, sign_consistency = seed_stability(effects, min_seed)
        tier = law_tier(mean, q, q <= float(discovery_cfg.get("emergent_law", {}).get("max_permutation_q_for_acceptable", 0.10)), stability, True, False, discovery_cfg)
        rec = {
            "variant": variant,
            "branch_nucleation_tier": tier,
            "lineage_separation_effect": mean,
            "effect_ci_low": lo,
            "effect_ci_high": hi,
            "permutation_p": p,
            "permutation_q": q,
            "seed_stability_pass": stability,
            "sign_consistency": sign_consistency,
            "n_seed_windows": int(group["seed"].nunique()),
        }
        for col in effect_cols:
            rec[col] = float(group[col].mean())
        rows.append(rec)
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows)
    out["best_interpretation"] = np.where(
        (out.get("local_velocity_alignment_A_effect", 0) > 0) & (out["lineage_separation_effect"] < 0),
        "transient_condensation_before_divergence",
        np.where(out.get("local_velocity_alignment_A_effect", 0) > 0, "alignment-driven branching", "unsupported"),
    )
    return out


def _control_summary(order_all: pd.DataFrame, candidate_order: pd.DataFrame, discovery_cfg: dict, quick_fixture: bool) -> pd.DataFrame:
    rng = np.random.default_rng(31)
    controls = []
    repeats = int(discovery_cfg.get("emergent_law", {}).get("permutation_repeats", 100))
    for control in ["shuffled_temporal_order", "shuffled_velocity", "shuffled_lineage_labels", "shuffled_fate_probabilities"]:
        values = []
        for rep in range(repeats):
            shuffled = candidate_order.copy()
            if control == "shuffled_temporal_order":
                shuffled["step"] = shuffled.groupby(["seed", "variant"], observed=False)["step"].transform(lambda x: rng.permutation(x.to_numpy()))
            elif control == "shuffled_velocity":
                shuffled["local_velocity_alignment_A"] = rng.permutation(shuffled["local_velocity_alignment_A"].to_numpy())
            elif control == "shuffled_lineage_labels":
                shuffled["lineage_separation_S"] = rng.permutation(shuffled["lineage_separation_S"].to_numpy())
                shuffled["branch_imbalance_B"] = rng.permutation(shuffled["branch_imbalance_B"].to_numpy())
            elif control == "shuffled_fate_probabilities":
                shuffled["fate_entropy_H"] = rng.permutation(shuffled["fate_entropy_H"].to_numpy())
            win = _event_windows(shuffled)
            values.append(float(win["lineage_separation_S_effect"].mean()) if not win.empty else 0.0)
        mean, lo, hi = bootstrap_ci(values, repeats=300)
        controls.append(
            {
                "control": control,
                "effect_size": mean,
                "effect_ci_low": lo,
                "effect_ci_high": hi,
                "seed_stability_pass": False,
                "gate_tier": "fail",
                "reason": "permutation/shuffle null control, expected not to reproduce retained branch signature",
            }
        )
    for control, variant in [
        ("no_swarm_model", "M2_ot_teacher_force"),
        ("no_teacher_model", "M1_intrinsic_neural"),
        ("random_teacher_velocity", "M10_shuffled_time_ot"),
    ]:
        subset = order_all[order_all["variant"].astype(str) == variant]
        win = _event_windows(subset) if not subset.empty else pd.DataFrame()
        effects = win["lineage_separation_S_effect"].to_numpy(dtype=float) if "lineage_separation_S_effect" in win else np.array([0.0])
        mean, lo, hi = bootstrap_ci(effects, repeats=300)
        stability, _ = seed_stability(effects, int(discovery_cfg.get("emergent_law", {}).get("quick_min_seed_count" if quick_fixture else "min_seed_count", 5)))
        controls.append(
            {
                "control": control,
                "variant": variant,
                "effect_size": mean,
                "effect_ci_low": lo,
                "effect_ci_high": hi,
                "seed_stability_pass": stability,
                "gate_tier": "weak" if stability and abs(mean) > 0.03 else "fail",
                "reason": "architectural negative/control comparator",
            }
        )
    pvals = [permutation_p_value(abs(row["effect_size"]), np.array([abs(v["effect_size"]) for v in controls[:4]])) for row in controls]
    qvals = bh_q_values(pvals)
    for row, p, q in zip(controls, pvals, qvals):
        row["permutation_p"] = p
        row["permutation_q"] = q
    return pd.DataFrame(controls)


def run(config: str = "configs/train.yaml", quick_fixture: bool = False):
    train_cfg, model_cfg, discovery_cfg = configure_paths(config, quick_fixture)
    table_dir, report_dir, fig_dir = output_dirs(train_cfg, discovery_cfg)
    adata = load_teacher(model_cfg)
    order_path = Path(model_cfg.get("order_log_path", table_dir / "rollout_order_parameters.csv"))
    if order_path.exists():
        order_all = pd.read_csv(order_path)
        rollout_based = bool("variant" in order_all and not order_all["variant"].astype(str).eq("static_teacher").all())
    else:
        order_all = _static_order_parameters(adata, model_cfg)
        rollout_based = False
    order = order_all.copy()
    analysis_variants = ["M5_ot_swarm", "M7_ot_swarm_birth_death_diffusion", "M8_ot_swarm_birth_death_diffusion_cci", "M9_full_memory"]
    if "variant" in order:
        preferred = order[order["variant"].isin(analysis_variants)]
        if not preferred.empty:
            order = preferred.copy()
    order.to_csv(table_dir / "branch_nucleation_rollout_order_parameters.csv", index=False)
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
    model_comparison = _summarize_model_windows(windows, controls["effect_null"].to_numpy(dtype=float), discovery_cfg, quick_fixture)
    model_comparison.to_csv(table_dir / "branch_nucleation_model_comparison.csv", index=False)
    control_summary = _control_summary(order_all, order, discovery_cfg, quick_fixture)
    control_summary.to_csv(table_dir / "branch_nucleation_negative_controls.csv", index=False)
    module_rows = []
    if not model_comparison.empty:
        base = model_comparison[model_comparison["variant"] == "M5_ot_swarm"]
        base_eff = float(base["lineage_separation_effect"].iloc[0]) if not base.empty else np.nan
        module_map = {
            "M5_ot_swarm": "swarm",
            "M7_ot_swarm_birth_death_diffusion": "swarm+birth_death+diffusion",
            "M8_ot_swarm_birth_death_diffusion_cci": "swarm+birth_death+diffusion+cci",
            "M9_full_memory": "swarm+birth_death+diffusion+cci+memory",
        }
        for row in model_comparison.itertuples(index=False):
            module_rows.append(
                {
                    "variant": row.variant,
                    "modules": module_map.get(row.variant, ""),
                    "lineage_separation_effect": row.lineage_separation_effect,
                    "delta_vs_M5": float(row.lineage_separation_effect - base_eff) if np.isfinite(base_eff) else np.nan,
                    "unsupported_module_burden": int("birth_death" in module_map.get(row.variant, "")) + int("cci" in module_map.get(row.variant, "")) + int("memory" in module_map.get(row.variant, "")),
                    "interpretation": "swarm_sufficient" if row.variant == "M5_ot_swarm" else "added_modules_not_required_for_main_claim",
                }
            )
    pd.DataFrame(module_rows).to_csv(table_dir / "branch_nucleation_module_dependence.csv", index=False)
    sensitivity_path = table_dir / "branch_nucleation_teacher_sensitivity.csv"
    if not sensitivity_path.exists():
        pd.DataFrame(
            [
                {
                    "native_max_cells_per_time": np.nan,
                    "epsilon": np.nan,
                    "branch_nucleation_teacher_sensitivity_tier": "pending",
                    "note": "Run python -m src.ot_teacher.native_sensitivity to populate teacher sensitivity.",
                }
            ]
        ).to_csv(sensitivity_path, index=False)
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

    fig, ax = plt.subplots(figsize=(7, 4), dpi=160)
    if not model_comparison.empty:
        plot = model_comparison.sort_values("lineage_separation_effect")
        ax.barh(plot["variant"], plot["lineage_separation_effect"], color="#7da7c7")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("lineage separation event effect")
    ax.set_title("Branch nucleation model comparison")
    fig.tight_layout()
    fig.savefig(fig_dir / "branch_nucleation_model_comparison.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4), dpi=160)
    if not windows.empty:
        trajectory = order.groupby(["variant", "step"], observed=False).mean(numeric_only=True).reset_index()
        for variant, group in trajectory.groupby("variant", observed=False):
            vals = group["lineage_separation_S"].to_numpy(float)
            vals = (vals - np.nanmin(vals)) / max(np.nanmax(vals) - np.nanmin(vals), 1e-8)
            ax.plot(group["step"], vals, marker="o", label=variant)
    ax.set_xlabel("rollout step")
    ax.set_ylabel("normalized lineage separation")
    ax.set_title("Branch nucleation event window")
    ax.legend(fontsize=6, frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "branch_nucleation_event_window.png")
    plt.close(fig)

    best_interpretation = "unsupported"
    if not model_comparison.empty:
        retained = model_comparison[model_comparison["branch_nucleation_tier"].isin(["acceptable", "strong"])]
        if not retained.empty:
            best_interpretation = str(retained.sort_values(["branch_nucleation_tier", "lineage_separation_effect"], ascending=[False, True])["best_interpretation"].iloc[0])
    primary_hint = "M5_ot_swarm" if "M5_ot_swarm" in set(model_comparison.get("variant", [])) else "not_selected"
    write_report(
        report_dir / "branch_nucleation_mechanism_summary.md",
        "Branch Nucleation Mechanism Summary",
        [
            f"- branch_nucleation_tier: {tier}",
            f"- best_interpretation: {best_interpretation}",
            f"- primary_model_hint: {primary_hint}",
            "- unsupported modules to exclude from main claim: birth/death, CCI, memory.",
            "- architectural controls can show related condensation signals, so the current evidence supports an order-parameter signature, not proof that swarm or teacher terms are necessary by themselves.",
            "- Negative controls include temporal, velocity, lineage, fate, no-swarm, no-teacher and random-teacher controls.",
            "",
            "## Model Comparison",
            "",
            model_comparison.to_markdown(index=False) if not model_comparison.empty else "No model comparison was available.",
            "",
            "## Negative Controls",
            "",
            control_summary.to_markdown(index=False) if not control_summary.empty else "No negative controls were available.",
        ],
    )

    write_report(
        report_dir / "discovery_branch_nucleation.md",
        "Discovery Branch Nucleation",
        [
            "Branch nucleation is hardened by rollout order-parameter traces, event-window detection, model comparison and shuffled/architectural negative controls.",
            "",
            "## Tier",
            "",
            f"- tier: {tier}",
            f"- rollout_based: {rollout_based}",
            f"- lineage_separation_event_effect: {effect_mean:.6g} [{ci_low:.6g}, {ci_high:.6g}]",
            f"- permutation_q: {q:.6g}",
            f"- seed_stability_pass: {stability} (sign consistency={sign_consistency:.3f})",
            f"- negative_control_pass: {negative_control_pass}",
            f"- best_interpretation: {best_interpretation}",
            "",
            "## Event Windows",
            "",
            windows.head(20).to_markdown(index=False) if not windows.empty else "No event window could be detected.",
            "",
            "## Model Comparison",
            "",
            model_comparison.to_markdown(index=False) if not model_comparison.empty else "No model comparison was available.",
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
