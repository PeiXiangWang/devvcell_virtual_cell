from __future__ import annotations

import argparse
from pathlib import Path

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


LAW = "birth_death"


def _density_bins(values: pd.Series, q: int = 8) -> pd.Series:
    ranked = values.rank(method="first")
    return pd.qcut(ranked, q=min(q, max(2, values.shape[0] // 20)), duplicates="drop")


def _event_rates(events: pd.DataFrame, frame: pd.DataFrame) -> pd.DataFrame:
    if events.empty:
        return pd.DataFrame()
    ev = events.copy()
    if "local_density" not in ev.columns:
        lineage_density = frame.groupby("lineage")["local_density"].mean()
        ev["local_density"] = ev["lineage"].map(lineage_density).fillna(float(frame["local_density"].mean()))
    if "n_agents_pre" not in ev.columns:
        ev["n_agents_pre"] = np.nan
    ev["density_bin"] = _density_bins(ev["local_density"], q=8).astype(str)
    rows = []
    for (seed, variant, density_bin), group in ev.groupby(["seed", "variant", "density_bin"], observed=False):
        births = int((group["event"] == "birth").sum())
        deaths = int((group["event"] == "death").sum())
        denom = float(group["n_agents_pre"].replace(0, np.nan).mean())
        if not np.isfinite(denom):
            denom = max(group.shape[0], 1)
        rows.append(
            {
                "seed": int(seed),
                "variant": variant,
                "density_bin": density_bin,
                "local_density": float(group["local_density"].mean()),
                "empirical_birth_rate_by_density_bin": births / denom,
                "empirical_death_rate_by_density_bin": deaths / denom,
                "empirical_net_growth_rate_by_density_bin": (births - deaths) / denom,
                "event_count": int(group.shape[0]),
            }
        )
    return pd.DataFrame(rows)


def _calibration(pred: pd.DataFrame, rates: pd.DataFrame) -> tuple[float, float, float, bool]:
    if rates.empty:
        return float("nan"), float("nan"), float("nan"), False
    pred_curve = pred.groupby("density_bin", observed=False).agg(predicted_net_growth_hazard=("predicted_net_growth_hazard", "mean")).reset_index()
    empirical = rates.groupby("density_bin", observed=False).agg(empirical_net_growth_rate_by_density_bin=("empirical_net_growth_rate_by_density_bin", "mean")).reset_index()
    merged = pred_curve.merge(empirical, on="density_bin", how="inner")
    if merged.shape[0] < 3:
        return float("nan"), float("nan"), float("nan"), False
    x = merged["predicted_net_growth_hazard"].to_numpy(dtype=float)
    y = merged["empirical_net_growth_rate_by_density_bin"].to_numpy(dtype=float)
    corr = float(np.corrcoef(x, y)[0, 1]) if np.std(x) > 0 and np.std(y) > 0 else 0.0
    slope = float(np.linalg.lstsq(np.c_[np.ones(x.size), x], y, rcond=None)[0][1])
    rmse = float(np.sqrt(np.mean((x - y) ** 2)))
    threshold = bool((merged["predicted_net_growth_hazard"].min() < 0 < merged["predicted_net_growth_hazard"].max()) or (merged["empirical_net_growth_rate_by_density_bin"].min() < 0 < merged["empirical_net_growth_rate_by_density_bin"].max()))
    return corr, slope, rmse, threshold


def run(config: str = "configs/train.yaml", quick_fixture: bool = False):
    train_cfg, model_cfg, discovery_cfg = configure_paths(config, quick_fixture)
    table_dir, report_dir, fig_dir = output_dirs(train_cfg, discovery_cfg)
    adata = load_teacher(model_cfg)
    seeds = seed_list(model_cfg, quick_fixture)
    frame = seedwise_feature_frame(adata, model_cfg, seeds)
    predictors = ["local_density", "ot_growth", "cell_cycle_score", "fate_entropy", "cci_signal"]
    stats = linear_effects(frame, "net_growth_hazard", predictors)

    pred = frame.copy()
    pred["density_bin"] = _density_bins(pred["local_density"], q=8).astype(str)
    curve = pred.groupby("density_bin", observed=False).agg(
        local_density=("local_density", "mean"),
        predicted_birth_hazard=("birth_hazard", "mean"),
        predicted_death_hazard=("death_hazard", "mean"),
        predicted_net_growth_hazard=("net_growth_hazard", "mean"),
    ).reset_index()
    curve.to_csv(table_dir / "birth_death_law.csv", index=False)

    events_path = Path(model_cfg.get("event_log_path", "tables/birth_death_event_log.csv"))
    events = pd.read_csv(events_path) if events_path.exists() else pd.DataFrame()
    rates = _event_rates(events, frame)
    rates.to_csv(table_dir / "birth_death_event_rates.csv", index=False)

    density_effects = []
    for seed, group in frame.groupby("seed", observed=False):
        density_effects.append(standardized_coef(group, "net_growth_hazard", "local_density", ["ot_growth", "cell_cycle_score", "fate_entropy", "cci_signal"]))
    stability_df = pd.DataFrame({"seed": seeds, "density_effect": density_effects[: len(seeds)]})
    if not rates.empty:
        count_seed = events.groupby(["seed", "event"]).size().unstack(fill_value=0)
        for col in ("birth", "death"):
            if col not in count_seed:
                count_seed[col] = 0
        count_seed["birth_death_ratio"] = count_seed["birth"] / np.maximum(count_seed["death"], 1)
        stability_df = stability_df.merge(count_seed.reset_index(), on="seed", how="left")
    stability_df.to_csv(table_dir / "birth_death_seed_stability.csv", index=False)

    corr, slope, rmse, threshold_like = _calibration(curve, rates)
    repeats = int(discovery_cfg.get("emergent_law", {}).get("permutation_repeats", 100))
    rng = np.random.default_rng(19)
    null_rows = []
    for rep in range(repeats):
        shuffled = frame.copy()
        shuffled["local_density"] = rng.permutation(shuffled["local_density"].to_numpy())
        coef = standardized_coef(shuffled, "net_growth_hazard", "local_density", ["ot_growth", "cell_cycle_score", "fate_entropy", "cci_signal"])
        null_rows.append({"permutation": rep, "density_effect_null": coef})
    null = pd.DataFrame(null_rows)
    null.to_csv(table_dir / "birth_death_density_permutation.csv", index=False)

    effect_mean, ci_low, ci_high = bootstrap_ci(density_effects, repeats=int(discovery_cfg.get("emergent_law", {}).get("bootstrap_repeats", 500)))
    p = permutation_p_value(effect_mean, null["density_effect_null"].to_numpy(dtype=float))
    q = bh_q_values([p])[0]
    min_seed = int(discovery_cfg.get("emergent_law", {}).get("quick_min_seed_count" if quick_fixture else "min_seed_count", 5))
    stability, sign_consistency = seed_stability(density_effects, min_seed)
    negative_control_pass = bool(q <= float(discovery_cfg.get("emergent_law", {}).get("max_permutation_q_for_acceptable", 0.10)))
    empirical_support = bool(np.isfinite(corr) and abs(corr) > 0.2)
    rollout_based = bool(not events.empty)
    directly_encoded = False
    tier = law_tier(effect_mean, q, negative_control_pass, stability, rollout_based, directly_encoded, discovery_cfg)
    if tier == "strong" and not empirical_support:
        tier = "acceptable"
    elif tier == "acceptable" and not empirical_support:
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
        str(table_dir / "birth_death_law.csv"),
        str(report_dir / "discovery_birth_death_law.md"),
    )

    fig, ax = plt.subplots(figsize=(6, 4), dpi=160)
    ax.plot(curve["local_density"], curve["predicted_birth_hazard"], marker="o", label="predicted birth")
    ax.plot(curve["local_density"], curve["predicted_death_hazard"], marker="o", label="predicted death")
    ax.plot(curve["local_density"], curve["predicted_net_growth_hazard"], marker="o", label="predicted net")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("local density")
    ax.set_ylabel("hazard")
    ax.set_title("Birth/death density phase curve")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "birth_density_phase_curve.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4), dpi=160)
    if not rates.empty:
        empirical = rates.groupby("density_bin", observed=False).agg(local_density=("local_density", "mean"), empirical_net=("empirical_net_growth_rate_by_density_bin", "mean")).sort_values("local_density")
        ax.plot(empirical["local_density"], empirical["empirical_net"], marker="o", color="#7da7c7", label="empirical net")
    ax.plot(curve["local_density"], curve["predicted_net_growth_hazard"], marker="o", color="#d98b73", label="predicted net")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("local density")
    ax.set_ylabel("net rate / hazard")
    ax.set_title("Event-rate calibration by density")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "birth_death_event_rate_bins.png")
    plt.close(fig)

    event_by_lineage = events.groupby(["lineage", "event"]).size().unstack(fill_value=0).reset_index() if not events.empty else pd.DataFrame()
    event_by_time = events.groupby(["time", "event"]).size().unstack(fill_value=0).reset_index() if not events.empty else pd.DataFrame()
    write_report(
        report_dir / "discovery_birth_death_law.md",
        "Discovery Birth Death Law",
        [
            "Birth/death hardening uses empirical stochastic event logs, predicted hazards, density permutation controls and seed-wise event stability.",
            "",
            "## Tier",
            "",
            f"- tier: {tier}",
            f"- density effect: {effect_mean:.6g} [{ci_low:.6g}, {ci_high:.6g}]",
            f"- permutation_q: {q:.6g}",
            f"- seed_stability_pass: {stability} (sign consistency={sign_consistency:.3f})",
            f"- empirical_hazard_correlation: {corr:.6g}",
            f"- calibration_slope: {slope:.6g}",
            f"- calibration_rmse: {rmse:.6g}",
            f"- carrying_capacity_like_threshold: {threshold_like}",
            f"- rollout_based_event_log: {rollout_based}",
            "",
            "## Hazard Regression",
            "",
            stats.to_markdown(index=False) if not stats.empty else "No stable regression fit.",
            "",
            "## Event Counts By Lineage",
            "",
            event_by_lineage.to_markdown(index=False) if not event_by_lineage.empty else "No event log available.",
            "",
            "## Event Counts By Time",
            "",
            event_by_time.head(20).to_markdown(index=False) if not event_by_time.empty else "No event log available.",
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
