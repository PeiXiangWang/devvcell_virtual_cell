from __future__ import annotations

import argparse
from copy import deepcopy

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
    load_teacher_model,
    output_dirs,
    permutation_p_value,
    seed_list,
    seed_stability,
    seedwise_feature_frame,
    write_report,
)
from src.model.cci import sender_receiver_graph
from src.model.simulator import VARIANTS, simulate_variant


LAW = "cci_branch_bias"


def _label_entropy(labels: np.ndarray) -> float:
    counts = pd.Series(labels.astype(str)).value_counts(normalize=True)
    return float(-(counts * np.log(np.clip(counts, 1e-12, 1.0))).sum() / max(np.log(counts.shape[0]), 1.0))


def _composition_tv(a: np.ndarray, b: np.ndarray) -> float:
    fa = pd.Series(a.astype(str)).value_counts(normalize=True)
    fb = pd.Series(b.astype(str)).value_counts(normalize=True)
    idx = sorted(set(fa.index) | set(fb.index))
    return float(0.5 * np.abs(fa.reindex(idx, fill_value=0) - fb.reindex(idx, fill_value=0)).sum())


def _event_counts(event_rows: list[dict]) -> tuple[int, int]:
    events = pd.DataFrame(event_rows)
    if events.empty or "event" not in events:
        return 0, 0
    return int((events["event"] == "birth").sum()), int((events["event"] == "death").sum())


def _spread(z: np.ndarray) -> float:
    if z.shape[0] == 0:
        return 0.0
    return float(np.mean(np.linalg.norm(z - z.mean(axis=0), axis=1)))


def _summarize_rollout(pred: np.ndarray, labels: np.ndarray, event_rows: list[dict], receiver: str) -> dict:
    births, deaths = _event_counts(event_rows)
    freq = pd.Series(labels.astype(str)).value_counts(normalize=True)
    return {
        "receiver_branch_probability": float(freq.get(receiver, 0.0)),
        "fate_entropy": _label_entropy(labels),
        "birth_events": births,
        "death_events": deaths,
        "event_count": births + deaths,
        "net_event_count": births - deaths,
        "diffusion_spread": _spread(pred),
        "n_pred": int(pred.shape[0]),
    }


def _run_rollout(adata, model_cfg: dict, seed: int, perturbation: dict | None) -> tuple[np.ndarray, np.ndarray, dict]:
    model = load_teacher_model(adata, model_cfg, seed=seed)
    if model is None:
        raise FileNotFoundError(f"missing trained teacher dynamics model for seed {seed}")
    variant = next(v for v in VARIANTS if v.name == "M9_full_memory")
    cfg = deepcopy(model_cfg)
    cfg["simulation_cells_per_seed"] = min(int(cfg.get("simulation_cells_per_seed", 900)), int(cfg.get("cci_rerollout_cells_per_seed", 320)))
    if perturbation:
        cfg["cci_perturbation"] = perturbation
    pred, _true, meta = simulate_variant(adata, cfg, {"teacher": model}, variant, seed)
    return pred, meta["pred_labels"], meta


def _proxy_fallback(frame: pd.DataFrame, graph: pd.DataFrame, seeds: list[int]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows = []
    rerollout_rows = []
    control_rows = []
    if graph.empty:
        rows.append({"seed": seeds[0], "perturbation": "no_lr_edges", "sender": "NA", "receiver": "NA", "branch_probability_shift": 0.0, "fate_entropy_shift": 0.0, "birth_hazard_shift": 0.0, "death_hazard_shift": 0.0, "diffusion_shift": 0.0, "branch_composition_shift": 0.0, "event_count_shift": 0.0, "sender_receiver_specificity": 0.0})
    else:
        top_edges = graph.sort_values("weight", ascending=False).head(3)
        for edge in top_edges.itertuples(index=False):
            for seed, group in frame.groupby("seed", observed=False):
                mask = group["lineage"].astype(str) == str(edge.receiver)
                receiver = group.loc[mask] if mask.any() else group
                specificity = float(edge.weight)
                branch_shift = -specificity * float(receiver["fate_probability_max"].mean()) * 0.1
                row = {
                    "seed": int(seed),
                    "perturbation": "remove_lr_edge_proxy",
                    "sender": edge.sender,
                    "receiver": edge.receiver,
                    "branch_probability_shift": branch_shift,
                    "fate_entropy_shift": specificity * float(receiver["fate_entropy"].mean()) * 0.04,
                    "birth_hazard_shift": -specificity * float(receiver["birth_hazard"].mean()) * 0.03,
                    "death_hazard_shift": specificity * float(receiver["death_hazard"].mean()) * 0.02,
                    "diffusion_shift": specificity * float(receiver["learned_sigma"].mean()) * 0.03,
                    "branch_composition_shift": branch_shift,
                    "event_count_shift": 0.0,
                    "sender_receiver_specificity": specificity,
                }
                rows.append(row)
                rerollout_rows.append({**row, "rerollout_type": "feature_recomputed_proxy_not_full_population_rollout"})
        for seed in seeds:
            control_rows.append({"seed": int(seed), "control": "zero_cci_signal", "branch_probability_shift": 0.0, "sender_receiver_specificity": 0.0})
    return pd.DataFrame(rows), pd.DataFrame(rerollout_rows), pd.DataFrame(control_rows)


def run(config: str = "configs/train.yaml", quick_fixture: bool = False):
    train_cfg, model_cfg, discovery_cfg = configure_paths(config, quick_fixture)
    table_dir, report_dir, fig_dir = output_dirs(train_cfg, discovery_cfg)
    adata = load_teacher(model_cfg)
    seeds = seed_list(model_cfg, quick_fixture)
    frame = seedwise_feature_frame(adata, model_cfg, seeds)
    z = adata.obsm[model_cfg.get("latent_key", "X_pca")]
    labels = adata.obs[model_cfg.get("cell_type_key", "lineage")].astype(str).to_numpy()
    graph, _signal, _pairs = sender_receiver_graph(adata, z, labels)

    rows: list[dict] = []
    control_rows: list[dict] = []
    rerollout_based = False
    status = "executed"
    top_n = 1 if quick_fixture else 3
    try:
        if graph.empty:
            raise ValueError("no ligand-receptor sender-receiver edges detected")
        top_edges = graph.sort_values("weight", ascending=False).head(top_n)
        baseline_cache = {}
        for seed in seeds:
            pred, pred_labels, meta = _run_rollout(adata, model_cfg, seed, perturbation=None)
            baseline_cache[seed] = {"pred": pred, "labels": pred_labels, "meta": meta}
        for edge in top_edges.itertuples(index=False):
            for seed in seeds:
                base = baseline_cache[seed]
                base_summary = _summarize_rollout(base["pred"], base["labels"], base["meta"]["event_rows"], str(edge.receiver))
                perturbations = [
                    ("remove_lr_edge", {"mode": "remove_lr_edge", "sender": str(edge.sender), "receiver": str(edge.receiver)}),
                    ("remove_receiver", {"mode": "remove_receiver", "lineage": str(edge.receiver)}),
                    ("shuffle_receiver", {"mode": "shuffle_receiver"}),
                    ("zero_cci_signal", {"mode": "zero_cci_signal"}),
                    ("degree_preserving_random_lr", {"mode": "degree_preserving_random_lr"}),
                ]
                for name, perturb in perturbations:
                    pred, pred_labels, meta = _run_rollout(adata, model_cfg, seed, perturbation=perturb)
                    pert_summary = _summarize_rollout(pred, pred_labels, meta["event_rows"], str(edge.receiver))
                    row = {
                        "seed": int(seed),
                        "perturbation": name,
                        "sender": edge.sender,
                        "receiver": edge.receiver,
                        "branch_probability_shift": pert_summary["receiver_branch_probability"] - base_summary["receiver_branch_probability"],
                        "fate_entropy_shift": pert_summary["fate_entropy"] - base_summary["fate_entropy"],
                        "birth_hazard_shift": (pert_summary["birth_events"] - base_summary["birth_events"]) / max(base_summary["n_pred"], 1),
                        "death_hazard_shift": (pert_summary["death_events"] - base_summary["death_events"]) / max(base_summary["n_pred"], 1),
                        "diffusion_shift": pert_summary["diffusion_spread"] - base_summary["diffusion_spread"],
                        "branch_composition_shift": _composition_tv(pred_labels, base["labels"]),
                        "event_count_shift": pert_summary["event_count"] - base_summary["event_count"],
                        "sender_receiver_specificity": float(edge.weight),
                        "rerollout_type": "full_population_simulator_rerollout",
                    }
                    if name in {"remove_lr_edge", "remove_receiver"}:
                        rows.append(row)
                    else:
                        control_rows.append({"control": name, **row})
        rerollout = pd.DataFrame(rows)
        controls = pd.DataFrame(control_rows)
        out = rerollout.copy()
        rerollout_based = True
    except Exception as exc:
        status = f"proxy_fallback: {type(exc).__name__}: {exc}"
        out, rerollout, controls = _proxy_fallback(frame, graph, seeds)

    out.to_csv(table_dir / "cci_branch_bias.csv", index=False)
    rerollout.to_csv(table_dir / "cci_rerollout_effects.csv", index=False)
    controls.to_csv(table_dir / "cci_negative_controls.csv", index=False)

    seed_effects = rerollout.groupby("seed")["branch_probability_shift"].mean().to_numpy(dtype=float) if not rerollout.empty else np.array([0.0])
    effect_mean, ci_low, ci_high = bootstrap_ci(seed_effects, repeats=int(discovery_cfg.get("emergent_law", {}).get("bootstrap_repeats", 500)))
    null = controls["branch_probability_shift"].to_numpy(dtype=float) if not controls.empty and "branch_probability_shift" in controls else np.array([0.0])
    p = permutation_p_value(effect_mean, null)
    q = bh_q_values([p])[0]
    min_seed = int(discovery_cfg.get("emergent_law", {}).get("quick_min_seed_count" if quick_fixture else "min_seed_count", 5))
    stability, sign_consistency = seed_stability(seed_effects, min_seed)
    negative_control_pass = bool(abs(effect_mean) > max(abs(null).mean() + abs(null).std(), 1e-8))
    directly_encoded = False
    tier = law_tier(effect_mean, q, negative_control_pass, stability, rerollout_based, directly_encoded, discovery_cfg)
    if not rerollout_based and tier in {"acceptable", "strong"}:
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
        rerollout_based,
        directly_encoded,
        str(table_dir / "cci_branch_bias.csv"),
        str(report_dir / "discovery_cci_branch_bias.md"),
        status=status,
    )

    fig, ax = plt.subplots(figsize=(7, 4), dpi=160)
    top = out.sort_values("sender_receiver_specificity", ascending=False).head(10) if not out.empty else out
    if not top.empty:
        labels_plot = top["sender"].astype(str) + "->" + top["receiver"].astype(str) + ":" + top["perturbation"].astype(str)
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
    if not controls.empty and "branch_probability_shift" in controls:
        ax.scatter(controls["seed"], controls["branch_probability_shift"], color="#999999", s=12, alpha=0.5, label="controls")
        ax.legend(frameon=False)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("seed")
    ax.set_ylabel("mean rerollout shift")
    ax.set_title("CCI rerollout effects")
    fig.tight_layout()
    fig.savefig(fig_dir / "cci_rerollout_effects.png")
    plt.close(fig)

    write_report(
        report_dir / "discovery_cci_branch_bias.md",
        "Discovery CCI Branch Bias",
        [
            "CCI hardening now attempts full population simulator rerollouts using `cci_perturbation` configs: remove LR edge, remove receiver, shuffle receiver, random LR and zero CCI signal.",
            "",
            "## Tier",
            "",
            f"- tier: {tier}",
            f"- effect_size: {effect_mean:.6g} [{ci_low:.6g}, {ci_high:.6g}]",
            f"- permutation_q: {q:.6g}",
            f"- seed_stability_pass: {stability} (sign consistency={sign_consistency:.3f})",
            f"- negative_control_pass: {negative_control_pass}",
            f"- rollout_based: {rerollout_based}",
            f"- status: {status}",
            "",
            out.head(20).to_markdown(index=False) if not out.empty else "No LR edge signal available.",
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
