"""Run the formal DevVCell Nature-level evidence pipeline.

The pipeline is configuration-driven and performs:

- input existence and schema validation,
- deterministic statistical tests and bootstrap confidence intervals,
- provenance capture with file hashes and runtime metadata,
- claim gate evaluation with explicit failure policy,
- publication-oriented tables, figures, manifest, and Methods text.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from devvcell.io import load_json, resolve_project_path, write_json  # noqa: E402


REQUIRED_COLUMNS: dict[str, dict[str, list[str]]] = {
    "stage_vulnerability": {
        "columns": ["stage", "stage_num", "window", "vulnerability_score"],
        "numeric": ["stage_num", "vulnerability_score"],
    },
    "stimulus_response": {
        "columns": [
            "tf_name",
            "devvcell_system",
            "stage",
            "stage_num",
            "target_stage",
            "stimulus_response_norm",
            "fate_displacement_from_baseline",
            "cell_level_feedback_cost_proxy",
            "cell_level_recovery_probability_proxy",
            "alignment_with_normal_transition",
            "stage_vulnerability_score",
        ],
        "numeric": [
            "stage_num",
            "stimulus_response_norm",
            "fate_displacement_from_baseline",
            "cell_level_feedback_cost_proxy",
            "cell_level_recovery_probability_proxy",
            "alignment_with_normal_transition",
            "stage_vulnerability_score",
        ],
    },
    "transition_metrics": {
        "columns": [
            "model",
            "split",
            "system",
            "src_stage",
            "tgt_stage",
            "n_pairs",
            "pair_latent_mse",
            "centroid_latent_mse",
            "rbf_mmd",
        ],
        "numeric": ["src_stage", "tgt_stage", "n_pairs", "pair_latent_mse", "centroid_latent_mse", "rbf_mmd"],
    },
    "transition_ablation": {
        "columns": [
            "pairing_method",
            "model",
            "mean_pair_latent_mse",
            "mean_centroid_latent_mse",
            "mean_rbf_mmd",
            "n_evaluations",
        ],
        "numeric": ["mean_pair_latent_mse", "mean_centroid_latent_mse", "mean_rbf_mmd", "n_evaluations"],
    },
    "stimulus_ablation": {
        "columns": ["variant", "mean_response", "max_response", "mean_recovery_probability", "n_tfs"],
        "numeric": ["mean_response", "max_response", "mean_recovery_probability", "n_tfs"],
    },
    "external_condition_metrics": {
        "columns": [
            "model",
            "condition",
            "gene",
            "context",
            "perturbation_label",
            "n_pairs",
            "pair_latent_mse",
            "centroid_latent_mse",
            "effect_cosine",
            "rbf_mmd",
        ],
        "numeric": ["n_pairs", "pair_latent_mse", "centroid_latent_mse", "effect_cosine", "rbf_mmd"],
    },
    "perturbation_priority": {
        "columns": [
            "tf_name",
            "devvcell_priority_score",
            "response_amplitude_proxy",
            "fate_displacement_proxy",
            "feedback_cost",
            "best_rescuer_tf",
            "best_rescue_fraction",
        ],
        "numeric": [
            "devvcell_priority_score",
            "response_amplitude_proxy",
            "fate_displacement_proxy",
            "feedback_cost",
            "best_rescue_fraction",
        ],
    },
}


@dataclass(frozen=True)
class PipelinePaths:
    root: Path
    tables: Path
    figures: Path
    methods: Path
    manifest: Path


class EvidencePipelineError(RuntimeError):
    """Raised when a configured hard failure condition is reached."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/nature_virtual_cell_evidence.json")
    parser.add_argument("--output-dir", default=None, help="Override output_dir in config.")
    parser.add_argument("--bootstrap-reps", type=int, default=None, help="Override bootstrap_reps in config.")
    parser.add_argument(
        "--fail-on-claim-threshold-failure",
        action="store_true",
        help="Treat claim gate failures as hard pipeline failures.",
    )
    return parser.parse_args()


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def make_paths(config: dict[str, Any], output_override: str | None) -> PipelinePaths:
    output_dir = resolve_project_path(output_override or config["output_dir"])
    paths = PipelinePaths(
        root=output_dir,
        tables=output_dir / "tables",
        figures=output_dir / "figures",
        methods=output_dir / "paper_methods_evidence_pipeline.md",
        manifest=output_dir / "evidence_manifest.json",
    )
    for path in [paths.root, paths.tables, paths.figures]:
        path.mkdir(parents=True, exist_ok=True)
    return paths


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT_ROOT)).replace("\\", "/")
    except ValueError:
        return str(path)


def input_file_inventory(config: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for name, path_like in config["input_tables"].items():
        path = resolve_project_path(path_like)
        rows.append(
            {
                "input_name": name,
                "path": rel(path),
                "exists": path.exists(),
                "size_bytes": int(path.stat().st_size) if path.exists() else None,
                "modified_utc": datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat() if path.exists() else None,
                "sha256": sha256_file(path) if path.exists() and path.is_file() else None,
            }
        )
    return rows


def read_csv_inputs(config: dict[str, Any]) -> dict[str, pd.DataFrame]:
    tables = {}
    for name in REQUIRED_COLUMNS:
        path = resolve_project_path(config["input_tables"][name])
        tables[name] = pd.read_csv(path)
    return tables


def validate_inputs(config: dict[str, Any], paths: PipelinePaths) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    failure_policy = config.get("failure_policy", {})
    inventory = input_file_inventory(config)
    inventory_df = pd.DataFrame(inventory)
    inventory_df.to_csv(paths.tables / "input_file_inventory.csv", index=False)

    missing = inventory_df.loc[~inventory_df["exists"], "path"].tolist()
    if missing and failure_policy.get("fail_on_missing_inputs", True):
        raise EvidencePipelineError(f"Missing required input files: {missing}")

    tables = read_csv_inputs(config)
    validation_rows = []
    min_rows = config.get("schema", {}).get("minimum_rows", {})
    nullable_numeric = config.get("schema", {}).get("nullable_numeric_columns", {})
    for name, frame in tables.items():
        spec = REQUIRED_COLUMNS[name]
        missing_cols = [col for col in spec["columns"] if col not in frame.columns]
        too_few_rows = len(frame) < int(min_rows.get(name, 1))
        invalid_numeric_cols = []
        missing_numeric_cols = []
        nullable_cols = set(nullable_numeric.get(name, []))
        for col in spec["numeric"]:
            if col not in frame.columns:
                continue
            values = pd.to_numeric(frame[col], errors="coerce")
            invalid_tokens = frame[col].notna() & values.isna()
            finite_values = values.dropna()
            has_infinite = bool(np.isinf(finite_values).any()) if not finite_values.empty else False
            if col in nullable_cols:
                if finite_values.empty or invalid_tokens.any() or has_infinite:
                    invalid_numeric_cols.append(col)
                continue
            if values.isna().any():
                missing_numeric_cols.append(col)
            if invalid_tokens.any() or has_infinite:
                invalid_numeric_cols.append(col)

        validation_rows.append(
            {
                "input_name": name,
                "n_rows": int(len(frame)),
                "n_columns": int(frame.shape[1]),
                "minimum_rows": int(min_rows.get(name, 1)),
                "missing_columns": ";".join(missing_cols),
                "nullable_numeric_columns": ";".join(sorted(nullable_cols)),
                "missing_required_numeric_columns": ";".join(missing_numeric_cols),
                "invalid_or_nonfinite_required_numeric_columns": ";".join(invalid_numeric_cols),
                "too_few_rows": bool(too_few_rows),
                "schema_pass": bool(
                    not missing_cols
                    and not too_few_rows
                    and not missing_numeric_cols
                    and not invalid_numeric_cols
                ),
            }
        )

    validation_df = pd.DataFrame(validation_rows)
    validation_df.to_csv(paths.tables / "input_schema_validation.csv", index=False)
    failed = validation_df[~validation_df["schema_pass"]]
    if not failed.empty and failure_policy.get("fail_on_schema_error", True):
        raise EvidencePipelineError(f"Input schema validation failed: {failed.to_dict(orient='records')}")

    return validation_df, tables


def validate_required_models(config: dict[str, Any], tables: dict[str, pd.DataFrame], paths: PipelinePaths) -> pd.DataFrame:
    checks = []
    transition_models = set(tables["transition_metrics"]["model"].astype(str))
    external_models = set(tables["external_condition_metrics"]["model"].astype(str))
    transition_ablation_models = set(tables["transition_ablation"]["model"].astype(str))

    for model in [config["primary_models"]["transition"], *config["baseline_models"]["transition"]]:
        checks.append(
            {
                "model_group": "transition",
                "model": model,
                "present_in_main_metrics": model in transition_models,
                "present_in_ablation": model in transition_ablation_models,
            }
        )
    for model in [config["primary_models"]["external_perturbation"], *config["baseline_models"]["external_perturbation"]]:
        checks.append(
            {
                "model_group": "external_perturbation",
                "model": model,
                "present_in_main_metrics": model in external_models,
                "present_in_ablation": None,
            }
        )

    checks_df = pd.DataFrame(checks)
    checks_df.to_csv(paths.tables / "required_model_validation.csv", index=False)
    missing = checks_df[checks_df["present_in_main_metrics"] == False]  # noqa: E712
    if not missing.empty and config.get("failure_policy", {}).get("fail_on_missing_required_models", True):
        raise EvidencePipelineError(f"Required model validation failed: {missing.to_dict(orient='records')}")
    return checks_df


def bootstrap_ci(values: np.ndarray, rng: np.random.Generator, reps: int) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan"), float("nan")
    if values.size == 1:
        return float(values[0]), float(values[0])
    samples = rng.choice(values, size=(reps, values.size), replace=True)
    means = samples.mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def paired_wilcoxon(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 2 or np.allclose(values, 0):
        return float("nan")
    try:
        return float(stats.wilcoxon(values).pvalue)
    except ValueError:
        return float("nan")


def mann_whitney_p(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return float("nan")
    return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)


def p_adjust_bh(p_values: list[float]) -> list[float]:
    p = np.asarray([np.nan if value is None else value for value in p_values], dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return q.tolist()
    valid_idx = np.where(valid)[0]
    order = valid_idx[np.argsort(p[valid])]
    ranked = p[order]
    m = len(ranked)
    adjusted = ranked * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    q[order] = adjusted
    return q.tolist()


def safe_relative_improvement(baseline: float, model: float) -> float:
    if not np.isfinite(baseline) or baseline == 0 or not np.isfinite(model):
        return float("nan")
    return float((baseline - model) / baseline)


def transition_analysis(config: dict[str, Any], tables: dict[str, pd.DataFrame], paths: PipelinePaths, rng: np.random.Generator) -> dict[str, Any]:
    reps = int(config["bootstrap_reps"])
    primary = config["primary_models"]["transition"]
    transition = tables["transition_metrics"].copy()

    model_summary = (
        transition.groupby("model", as_index=False)
        .agg(
            n_strata=("pair_latent_mse", "size"),
            n_pairs_total=("n_pairs", "sum"),
            mean_pair_latent_mse=("pair_latent_mse", "mean"),
            mean_centroid_latent_mse=("centroid_latent_mse", "mean"),
            mean_rbf_mmd=("rbf_mmd", "mean"),
        )
        .sort_values("mean_pair_latent_mse")
    )
    model_summary.to_csv(paths.tables / "transition_model_performance.csv", index=False)

    wide = transition.pivot_table(
        index=["split", "system", "src_stage", "tgt_stage"],
        columns="model",
        values="pair_latent_mse",
        aggfunc="mean",
    )
    paired_rows = []
    primary_mean = float(model_summary.loc[model_summary["model"] == primary, "mean_pair_latent_mse"].iloc[0])
    for competitor in config["baseline_models"]["transition"]:
        if competitor not in wide.columns or primary not in wide.columns:
            continue
        diff = (wide[competitor] - wide[primary]).dropna().to_numpy(dtype=float)
        ci_low, ci_high = bootstrap_ci(diff, rng, reps)
        competitor_mean = float(model_summary.loc[model_summary["model"] == competitor, "mean_pair_latent_mse"].iloc[0])
        paired_rows.append(
            {
                "claim": "stage_conditioned_transition",
                "primary_model": primary,
                "competitor_model": competitor,
                "n_paired_strata": int(len(diff)),
                "primary_mean_mse": primary_mean,
                "competitor_mean_mse": competitor_mean,
                "mean_competitor_minus_primary_mse": float(np.mean(diff)) if diff.size else float("nan"),
                "ci95_low": ci_low,
                "ci95_high": ci_high,
                "relative_improvement_vs_competitor": safe_relative_improvement(competitor_mean, primary_mean),
                "primary_better_strata": int(np.sum(diff > 0)) if diff.size else 0,
                "wilcoxon_p": paired_wilcoxon(diff),
            }
        )

    paired_df = pd.DataFrame(paired_rows)
    if not paired_df.empty:
        paired_df["wilcoxon_q_bh"] = p_adjust_bh(paired_df["wilcoxon_p"].tolist())
    paired_df.to_csv(paths.tables / "transition_paired_statistical_tests.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    plot_df = model_summary.sort_values("mean_pair_latent_mse")
    ax.bar(plot_df["model"], plot_df["mean_pair_latent_mse"], color="#4f7697")
    ax.set_ylabel("Mean heldout latent MSE")
    ax.set_xlabel("Model")
    ax.set_title("Formal transition benchmark")
    ax.grid(axis="y", alpha=0.25)
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(paths.figures / "figure_transition_formal_benchmark.png", dpi=240)
    plt.close(fig)

    return {
        "primary_model": primary,
        "primary_mean_mse": primary_mean,
        "model_summary_rows": int(len(model_summary)),
        "paired_tests_rows": int(len(paired_df)),
    }


def competence_analysis(config: dict[str, Any], tables: dict[str, pd.DataFrame], paths: PipelinePaths, rng: np.random.Generator) -> dict[str, Any]:
    reps = int(config["bootstrap_reps"])
    stimulus = tables["stimulus_response"].copy()
    stage = tables["stage_vulnerability"][["stage", "stage_num", "window", "vulnerability_score"]].copy()
    window = config["competence_window"]
    stimulus = stimulus.merge(stage, on=["stage", "stage_num"], how="left", suffixes=("", "_stage_table"))
    stimulus["is_competence_window"] = stimulus["stage_num"].between(int(window["stage_min"]), int(window["stage_max"]))
    stimulus["inverse_recovery_probability"] = 1.0 - pd.to_numeric(stimulus["cell_level_recovery_probability_proxy"], errors="coerce")

    metrics = [
        "stimulus_response_norm",
        "fate_displacement_from_baseline",
        "cell_level_feedback_cost_proxy",
        "inverse_recovery_probability",
        "alignment_with_normal_transition",
    ]

    group_rows = []
    test_rows = []
    for metric in metrics:
        inside = pd.to_numeric(stimulus.loc[stimulus["is_competence_window"], metric], errors="coerce").dropna().to_numpy(dtype=float)
        outside = pd.to_numeric(stimulus.loc[~stimulus["is_competence_window"], metric], errors="coerce").dropna().to_numpy(dtype=float)
        for label, values in [("competence_window", inside), ("outside", outside)]:
            ci_low, ci_high = bootstrap_ci(values, rng, reps)
            group_rows.append(
                {
                    "metric": metric,
                    "group": label,
                    "n": int(values.size),
                    "mean": float(values.mean()) if values.size else float("nan"),
                    "median": float(np.median(values)) if values.size else float("nan"),
                    "ci95_low": ci_low,
                    "ci95_high": ci_high,
                }
            )
        diff_boot = []
        if inside.size and outside.size:
            diff_boot = (
                rng.choice(inside, size=(reps, inside.size), replace=True).mean(axis=1)
                - rng.choice(outside, size=(reps, outside.size), replace=True).mean(axis=1)
            )
        test_rows.append(
            {
                "claim": "competence_window_modulates_response",
                "metric": metric,
                "n_inside": int(inside.size),
                "n_outside": int(outside.size),
                "mean_inside": float(inside.mean()) if inside.size else float("nan"),
                "mean_outside": float(outside.mean()) if outside.size else float("nan"),
                "mean_difference_inside_minus_outside": float(inside.mean() - outside.mean()) if inside.size and outside.size else float("nan"),
                "ci95_low": float(np.percentile(diff_boot, 2.5)) if len(diff_boot) else float("nan"),
                "ci95_high": float(np.percentile(diff_boot, 97.5)) if len(diff_boot) else float("nan"),
                "mann_whitney_p": mann_whitney_p(inside, outside),
            }
        )

    group_df = pd.DataFrame(group_rows)
    tests_df = pd.DataFrame(test_rows)
    tests_df["mann_whitney_q_bh"] = p_adjust_bh(tests_df["mann_whitney_p"].tolist())

    corr_rows = []
    for metric in metrics:
        corr_df = stimulus[["stage_vulnerability_score", metric]].replace([np.inf, -np.inf], np.nan).dropna()
        for method in ["pearson", "spearman"]:
            if len(corr_df) >= 3 and corr_df["stage_vulnerability_score"].nunique() > 1 and corr_df[metric].nunique() > 1:
                result = stats.pearsonr(corr_df["stage_vulnerability_score"], corr_df[metric]) if method == "pearson" else stats.spearmanr(corr_df["stage_vulnerability_score"], corr_df[metric])
                r_value = float(result.statistic)
                p_value = float(result.pvalue)
            else:
                r_value = float("nan")
                p_value = float("nan")
            corr_rows.append(
                {
                    "claim": "stage_vulnerability_tracks_response",
                    "metric": metric,
                    "method": method,
                    "n": int(len(corr_df)),
                    "correlation": r_value,
                    "p_value": p_value,
                }
            )
    corr_tests = pd.DataFrame(corr_rows)
    corr_tests["q_value_bh"] = p_adjust_bh(corr_tests["p_value"].tolist())

    stage_summary = (
        stimulus.groupby(["stage", "stage_num", "is_competence_window"], as_index=False)
        .agg(
            mean_response=("stimulus_response_norm", "mean"),
            mean_fate_displacement=("fate_displacement_from_baseline", "mean"),
            mean_recovery_cost=("cell_level_feedback_cost_proxy", "mean"),
            vulnerability_score=("stage_vulnerability_score", "mean"),
            n_tests=("stimulus_response_norm", "size"),
        )
        .sort_values("stage_num")
    )

    group_df.to_csv(paths.tables / "competence_window_group_statistics.csv", index=False)
    tests_df.to_csv(paths.tables / "competence_window_statistical_tests.csv", index=False)
    corr_tests.to_csv(paths.tables / "vulnerability_response_correlation_tests.csv", index=False)
    stage_summary.to_csv(paths.tables / "stage_response_summary.csv", index=False)

    fig, ax1 = plt.subplots(figsize=(8.4, 4.8))
    colors = ["#b4584d" if is_window else "#5d7f9d" for is_window in stage_summary["is_competence_window"]]
    ax1.bar(stage_summary["stage_num"].astype(str), stage_summary["mean_response"], color=colors, alpha=0.9)
    ax1.set_xlabel("Theiler stage")
    ax1.set_ylabel("Mean TF stimulus response")
    ax1.grid(axis="y", alpha=0.25)
    ax2 = ax1.twinx()
    ax2.plot(stage_summary["stage_num"].astype(str), stage_summary["vulnerability_score"], color="#1f1f1f", marker="o", linewidth=1.8)
    ax2.set_ylabel("Stage vulnerability")
    ax1.set_title("Competence-window perturbation response")
    fig.tight_layout()
    fig.savefig(paths.figures / "figure_competence_window_response.png", dpi=240)
    plt.close(fig)

    return {
        "stimulus_rows": int(len(stimulus)),
        "competence_window_rows": int(stimulus["is_competence_window"].sum()),
        "outside_rows": int((~stimulus["is_competence_window"]).sum()),
    }


def fate_recovery_analysis(config: dict[str, Any], tables: dict[str, pd.DataFrame], paths: PipelinePaths) -> dict[str, Any]:
    stimulus = tables["stimulus_response"].copy()
    priority = tables["perturbation_priority"][
        ["tf_name", "best_rescuer_tf", "best_rescue_fraction", "devvcell_priority_score"]
    ].drop_duplicates("tf_name").rename(
        columns={"devvcell_priority_score": "priority_score_from_priority_table"}
    )
    screen = stimulus.merge(priority, on="tf_name", how="left")
    if "devvcell_priority_score" not in screen.columns:
        screen["devvcell_priority_score"] = screen["priority_score_from_priority_table"]
    else:
        screen["devvcell_priority_score"] = pd.to_numeric(
            screen["devvcell_priority_score"], errors="coerce"
        ).fillna(pd.to_numeric(screen["priority_score_from_priority_table"], errors="coerce"))
    screen["response_amplitude"] = pd.to_numeric(screen["stimulus_response_norm"], errors="coerce")
    screen["fate_displacement"] = pd.to_numeric(screen["fate_displacement_from_baseline"], errors="coerce")
    screen["recovery_cost"] = pd.to_numeric(screen["cell_level_feedback_cost_proxy"], errors="coerce")
    screen["recovery_probability"] = pd.to_numeric(screen["cell_level_recovery_probability_proxy"], errors="coerce")
    screen["rescue_candidate"] = screen["best_rescuer_tf"].fillna("not_observed")
    screen["rescue_candidate_status"] = np.where(
        screen["rescue_candidate"].astype(str).eq("not_observed"),
        "not_observed",
        "observed_in_rescue_table",
    )

    screen_cols = [
        "tf_name",
        "devvcell_system",
        "stage",
        "stage_num",
        "target_stage",
        "response_amplitude",
        "fate_displacement",
        "recovery_cost",
        "recovery_probability",
        "alignment_with_normal_transition",
        "devvcell_priority_score",
        "rescue_candidate",
        "best_rescue_fraction",
        "rescue_candidate_status",
        "matched_global_grn_targets",
        "matched_system_targets",
    ]
    screen_out = screen[screen_cols].sort_values(["fate_displacement", "recovery_cost"], ascending=False)
    tf_summary = (
        screen_out.groupby("tf_name", as_index=False)
        .agg(
            n_system_stage_tests=("tf_name", "size"),
            mean_response_amplitude=("response_amplitude", "mean"),
            max_response_amplitude=("response_amplitude", "max"),
            mean_fate_displacement=("fate_displacement", "mean"),
            mean_recovery_cost=("recovery_cost", "mean"),
            mean_recovery_probability=("recovery_probability", "mean"),
            mean_priority_score=("devvcell_priority_score", "mean"),
            rescue_candidate=("rescue_candidate", "first"),
            best_rescue_fraction=("best_rescue_fraction", "max"),
        )
        .sort_values(["mean_fate_displacement", "mean_recovery_cost"], ascending=False)
    )
    screen_out.to_csv(paths.tables / "fate_recovery_virtual_screen.csv", index=False)
    tf_summary.to_csv(paths.tables / "fate_recovery_tf_summary.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    plot_df = tf_summary.head(12).copy()
    max_disp = max(float(plot_df["mean_fate_displacement"].max()), 1e-12)
    scatter = ax.scatter(
        plot_df["mean_response_amplitude"],
        plot_df["mean_recovery_cost"],
        s=80 + 380 * plot_df["mean_fate_displacement"] / max_disp,
        c=plot_df["mean_recovery_probability"],
        cmap="viridis_r",
        edgecolor="white",
        linewidth=0.8,
    )
    for row in plot_df.itertuples(index=False):
        ax.text(row.mean_response_amplitude, row.mean_recovery_cost, row.tf_name, fontsize=8)
    ax.set_xlabel("Mean response amplitude")
    ax.set_ylabel("Mean recovery cost")
    ax.set_title("Fate displacement and recovery cost")
    ax.grid(alpha=0.25)
    fig.colorbar(scatter, ax=ax, fraction=0.04, pad=0.02, label="Recovery probability")
    fig.tight_layout()
    fig.savefig(paths.figures / "figure_fate_recovery_virtual_screen.png", dpi=240)
    plt.close(fig)

    return {
        "screen_rows": int(len(screen_out)),
        "tf_summary_rows": int(len(tf_summary)),
        "top_tf": str(tf_summary.iloc[0]["tf_name"]),
    }


def external_perturbation_analysis(config: dict[str, Any], tables: dict[str, pd.DataFrame], paths: PipelinePaths, rng: np.random.Generator) -> dict[str, Any]:
    reps = int(config["bootstrap_reps"])
    primary = config["primary_models"]["external_perturbation"]
    external = tables["external_condition_metrics"].copy()

    model_summary = (
        external.groupby("model", as_index=False)
        .agg(
            n_conditions=("condition", "nunique"),
            n_pairs_total=("n_pairs", "sum"),
            mean_pair_latent_mse=("pair_latent_mse", "mean"),
            mean_centroid_latent_mse=("centroid_latent_mse", "mean"),
            mean_effect_cosine=("effect_cosine", "mean"),
            mean_rbf_mmd=("rbf_mmd", "mean"),
        )
        .sort_values("mean_pair_latent_mse")
    )

    wide = external.pivot_table(
        index=["condition", "gene", "context", "perturbation_label"],
        columns="model",
        values="pair_latent_mse",
        aggfunc="mean",
    )
    primary_mean = float(model_summary.loc[model_summary["model"] == primary, "mean_pair_latent_mse"].iloc[0])
    paired_rows = []
    for competitor in config["baseline_models"]["external_perturbation"]:
        if competitor not in wide.columns or primary not in wide.columns:
            continue
        diff = (wide[competitor] - wide[primary]).dropna().to_numpy(dtype=float)
        ci_low, ci_high = bootstrap_ci(diff, rng, reps)
        competitor_mean = float(model_summary.loc[model_summary["model"] == competitor, "mean_pair_latent_mse"].iloc[0])
        paired_rows.append(
            {
                "claim": "external_perturbation_calibration",
                "primary_model": primary,
                "competitor_model": competitor,
                "n_paired_conditions": int(diff.size),
                "primary_mean_mse": primary_mean,
                "competitor_mean_mse": competitor_mean,
                "mean_competitor_minus_primary_mse": float(diff.mean()) if diff.size else float("nan"),
                "ci95_low": ci_low,
                "ci95_high": ci_high,
                "relative_improvement_vs_competitor": safe_relative_improvement(competitor_mean, primary_mean),
                "primary_better_conditions": int(np.sum(diff > 0)) if diff.size else 0,
                "wilcoxon_p": paired_wilcoxon(diff),
            }
        )

    paired_df = pd.DataFrame(paired_rows)
    if not paired_df.empty:
        paired_df["wilcoxon_q_bh"] = p_adjust_bh(paired_df["wilcoxon_p"].tolist())

    centered = external.copy()
    centered["condition_mean_pair_latent_mse"] = centered.groupby("condition")["pair_latent_mse"].transform("mean")
    centered["condition_centered_pair_latent_mse"] = centered["pair_latent_mse"] - centered["condition_mean_pair_latent_mse"]
    centered_summary = (
        centered.groupby("model", as_index=False)
        .agg(
            mean_condition_centered_pair_latent_mse=("condition_centered_pair_latent_mse", "mean"),
            mean_pair_latent_mse=("pair_latent_mse", "mean"),
            mean_effect_cosine=("effect_cosine", "mean"),
            n_conditions=("condition", "nunique"),
        )
        .sort_values("mean_pair_latent_mse")
    )

    model_summary.to_csv(paths.tables / "external_perturbation_model_performance.csv", index=False)
    paired_df.to_csv(paths.tables / "external_perturbation_paired_statistical_tests.csv", index=False)
    centered.to_csv(paths.tables / "external_perturbation_condition_centered_metrics.csv", index=False)
    centered_summary.to_csv(paths.tables / "external_perturbation_bias_control_summary.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    plot_df = paired_df.sort_values("mean_competitor_minus_primary_mse", ascending=False)
    ax.bar(plot_df["competitor_model"], plot_df["mean_competitor_minus_primary_mse"], color="#5f8f6b")
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_ylabel("Competitor MSE - primary MSE")
    ax.set_xlabel("External baseline")
    ax.set_title("External Perturb-seq paired calibration evidence")
    ax.grid(axis="y", alpha=0.25)
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(paths.figures / "figure_external_perturbation_paired_tests.png", dpi=240)
    plt.close(fig)

    return {
        "primary_model": primary,
        "primary_mean_mse": primary_mean,
        "n_external_conditions": int(external["condition"].nunique()),
    }


def ablation_analysis(config: dict[str, Any], tables: dict[str, pd.DataFrame], paths: PipelinePaths) -> dict[str, Any]:
    transition_ablation = tables["transition_ablation"].copy()
    stimulus_ablation = tables["stimulus_ablation"].copy()
    transition_ablation.to_csv(paths.tables / "transition_ablation_formal_summary.csv", index=False)
    stimulus_ablation.to_csv(paths.tables / "stimulus_ablation_formal_summary.csv", index=False)
    return {
        "transition_ablation_rows": int(len(transition_ablation)),
        "stimulus_ablation_rows": int(len(stimulus_ablation)),
        "best_transition_ablation": transition_ablation.sort_values("mean_pair_latent_mse").iloc[0].to_dict(),
    }


def build_contract_table(paths: PipelinePaths) -> pd.DataFrame:
    rows = [
        {
            "virtual_cell_component": "normal_next_state",
            "implemented_object": "stage/system-conditioned residual transition operator",
            "formal_output": "predicted latent state at target stage",
            "evidence_table": "transition_model_performance.csv",
            "status": "implemented",
        },
        {
            "virtual_cell_component": "perturbed_state",
            "implemented_object": "TF/GRN stimulus-conditioned latent projection",
            "formal_output": "perturbed latent state proxy",
            "evidence_table": "fate_recovery_virtual_screen.csv",
            "status": "implemented_computational_proxy",
        },
        {
            "virtual_cell_component": "fate_displacement",
            "implemented_object": "distance from normal transition baseline",
            "formal_output": "fate_displacement",
            "evidence_table": "fate_recovery_virtual_screen.csv",
            "status": "implemented_computational_proxy",
        },
        {
            "virtual_cell_component": "recovery_cost",
            "implemented_object": "cell-level feedback cost and recovery probability",
            "formal_output": "recovery_cost; recovery_probability",
            "evidence_table": "fate_recovery_virtual_screen.csv",
            "status": "implemented_computational_proxy",
        },
        {
            "virtual_cell_component": "perturbation_calibration",
            "implemented_object": "external scPerturb guide-transfer benchmark",
            "formal_output": "paired improvement over naive baselines",
            "evidence_table": "external_perturbation_paired_statistical_tests.csv",
            "status": "implemented_first_dataset",
        },
    ]
    out = pd.DataFrame(rows)
    out.to_csv(paths.tables / "virtual_cell_component_contract.csv", index=False)
    return out


def build_claim_gates(config: dict[str, Any], summaries: dict[str, Any], tables: dict[str, pd.DataFrame], paths: PipelinePaths) -> pd.DataFrame:
    thresholds = config["claim_thresholds"]

    transition_tests = pd.read_csv(paths.tables / "transition_paired_statistical_tests.csv")
    external_tests = pd.read_csv(paths.tables / "external_perturbation_paired_statistical_tests.csv")
    competence_tests = pd.read_csv(paths.tables / "competence_window_statistical_tests.csv")

    transition_vs_ridge = transition_tests[transition_tests["competitor_model"] == "ridge"]
    external_vs_identity = external_tests[external_tests["competitor_model"] == "identity"]
    stimulus_rows = summaries["competence"]["stimulus_rows"]
    competence_rows = summaries["competence"]["competence_window_rows"]

    def first_value(frame: pd.DataFrame, column: str, default: float = float("nan")) -> float:
        if frame.empty or column not in frame.columns:
            return default
        return float(frame.iloc[0][column])

    gate_rows = [
        {
            "claim": "stage_conditioned_transition_operator",
            "gate": "minimum paired heldout strata",
            "observed": int(first_value(transition_vs_ridge, "n_paired_strata", 0)),
            "threshold": int(thresholds["minimum_transition_paired_strata"]),
            "pass": int(first_value(transition_vs_ridge, "n_paired_strata", 0)) >= int(thresholds["minimum_transition_paired_strata"]),
            "hard_failure": False,
        },
        {
            "claim": "stage_conditioned_transition_operator",
            "gate": "relative improvement versus ridge",
            "observed": first_value(transition_vs_ridge, "relative_improvement_vs_competitor"),
            "threshold": float(thresholds["minimum_transition_relative_improvement"]),
            "pass": first_value(transition_vs_ridge, "relative_improvement_vs_competitor") >= float(thresholds["minimum_transition_relative_improvement"]),
            "hard_failure": False,
        },
        {
            "claim": "competence_window_response",
            "gate": "minimum stimulus rows",
            "observed": int(stimulus_rows),
            "threshold": int(thresholds["minimum_stimulus_rows"]),
            "pass": int(stimulus_rows) >= int(thresholds["minimum_stimulus_rows"]),
            "hard_failure": False,
        },
        {
            "claim": "competence_window_response",
            "gate": "minimum competence window rows",
            "observed": int(competence_rows),
            "threshold": int(thresholds["minimum_competence_window_rows"]),
            "pass": int(competence_rows) >= int(thresholds["minimum_competence_window_rows"]),
            "hard_failure": False,
        },
        {
            "claim": "external_perturbation_calibration",
            "gate": "minimum external conditions",
            "observed": int(summaries["external"]["n_external_conditions"]),
            "threshold": int(thresholds["minimum_external_conditions"]),
            "pass": int(summaries["external"]["n_external_conditions"]) >= int(thresholds["minimum_external_conditions"]),
            "hard_failure": False,
        },
        {
            "claim": "external_perturbation_calibration",
            "gate": "relative improvement versus identity",
            "observed": first_value(external_vs_identity, "relative_improvement_vs_competitor"),
            "threshold": float(thresholds["minimum_external_relative_improvement_vs_identity"]),
            "pass": first_value(external_vs_identity, "relative_improvement_vs_competitor") >= float(thresholds["minimum_external_relative_improvement_vs_identity"]),
            "hard_failure": False,
        },
        {
            "claim": "competence_window_response",
            "gate": "statistical tests completed",
            "observed": int(len(competence_tests)),
            "threshold": 5,
            "pass": len(competence_tests) >= 5,
            "hard_failure": False,
        },
    ]
    gates = pd.DataFrame(gate_rows)
    gates.to_csv(paths.tables / "claim_gate_matrix.csv", index=False)

    if config.get("failure_policy", {}).get("fail_on_claim_threshold_failure", False):
        failed = gates[~gates["pass"]]
        if not failed.empty:
            raise EvidencePipelineError(f"Claim threshold gate failed: {failed.to_dict(orient='records')}")
    return gates


def write_methods_text(config: dict[str, Any], summaries: dict[str, Any], paths: PipelinePaths, config_path: Path) -> None:
    text = f"""# Paper Methods: DevVCell Formal Evidence Pipeline

The DevVCell evidence pipeline was run with configuration `{rel(config_path)}`. The pipeline validates all required input files, checks table schemas, verifies required primary and baseline models, computes paired statistical comparisons, records full provenance, and writes a manifest linking every claim to output tables and figures.

## Input Validation

For each required result table, the pipeline checks file existence, SHA-256 hash, minimum row count, required columns, and numeric validity for required metrics. Required numeric columns must be present, parseable, non-missing, and finite unless explicitly listed under `schema.nullable_numeric_columns`; nullable columns must still contain at least one finite value and no invalid tokens. Missing inputs, schema errors, missing required models, and invalid required metrics are treated as hard failures according to the configured failure policy.

## Developmental Transition Evidence

Heldout cell-level transition metrics are grouped by model and evaluated using paired system-stage strata. The primary model is `{config['primary_models']['transition']}`. For each baseline, the pipeline computes mean paired MSE difference, bootstrap 95% confidence interval with {config['bootstrap_reps']} resamples, relative improvement, the number of strata in which the primary model is better, and a paired Wilcoxon signed-rank p-value with Benjamini-Hochberg correction.

## Competence-Window Evidence

TF/GRN stimulus-response rows are stratified by the configured competence window, TS{config['competence_window']['stage_min']}--TS{config['competence_window']['stage_max']}. For response amplitude, fate displacement, recovery cost, inverse recovery probability, and alignment with normal transition, the pipeline reports group means, medians, bootstrap confidence intervals, window-versus-outside contrasts, Mann-Whitney tests, and false-discovery-rate adjusted q-values. It also tests Pearson and Spearman correlations between stage vulnerability and stimulus-response metrics.

## Fate and Recovery Virtual Screen

The pipeline constructs a TF-system-stage table containing response amplitude, fate displacement, recovery cost, recovery probability, priority score, and rescue candidate annotations. Rescue candidates are carried forward only as computational hypotheses unless supported by external perturbation or wet-lab rescue assays.

## External Perturbation Calibration

External scPerturb guide-transfer metrics are evaluated by condition-paired comparisons. The primary external model is `{config['primary_models']['external_perturbation']}`. The pipeline computes paired improvements against configured baselines, bootstrap confidence intervals, Wilcoxon signed-rank p-values, effect-cosine summaries, and condition-centered metrics that remove condition-level average difficulty before comparing models.

## Claim Gates

Claim gates are recorded in `tables/claim_gate_matrix.csv`. They are designed to make the current evidence boundary explicit. By default, claim-gate failures are reported but do not abort the pipeline unless `fail_on_claim_threshold_failure` is enabled.

## Interpretation Boundary

The current TF/GRN stimulus, fate-displacement, and recovery-cost outputs are computational proxy readouts. They are suitable for formal hypothesis generation and manuscript-facing evidence tracking, but they do not replace direct CellOT/GEARS/scGen/CPA/foundation-model baselines or stage-dependent perturbation experiments.
"""
    paths.methods.write_text(text, encoding="utf-8")


def build_manifest(
    config: dict[str, Any],
    config_path: Path,
    paths: PipelinePaths,
    validation: pd.DataFrame,
    model_validation: pd.DataFrame,
    gates: pd.DataFrame,
    summaries: dict[str, Any],
    started_utc: str,
) -> dict[str, Any]:
    output_files = []
    for folder in [paths.tables, paths.figures]:
        for path in sorted(folder.glob("*")):
            if path.is_file():
                output_files.append(
                    {
                        "path": rel(path),
                        "size_bytes": int(path.stat().st_size),
                        "sha256": sha256_file(path),
                    }
                )
    output_files.append(
        {
            "path": rel(paths.methods),
            "size_bytes": int(paths.methods.stat().st_size),
            "sha256": sha256_file(paths.methods),
        }
    )

    manifest = {
        "analysis": config["analysis"],
        "started_utc": started_utc,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "command": " ".join(sys.argv),
        "cwd": str(PROJECT_ROOT),
        "config_path": rel(config_path),
        "config_sha256": sha256_file(config_path),
        "input_files": input_file_inventory(config),
        "validation_summary": {
            "schema_pass": bool(validation["schema_pass"].all()),
            "required_models_pass": bool(model_validation["present_in_main_metrics"].fillna(True).all()),
            "claim_gates_pass": bool(gates["pass"].all()),
        },
        "summaries": summaries,
        "claim_gates": gates.to_dict(orient="records"),
        "output_files": output_files,
    }
    write_json(paths.manifest, manifest)
    return manifest


def main() -> None:
    args = parse_args()
    started_utc = datetime.now(timezone.utc).isoformat()
    config_path = resolve_project_path(args.config)
    config = load_json(config_path)
    if args.output_dir is not None:
        config["output_dir"] = args.output_dir
    if args.bootstrap_reps is not None:
        config["bootstrap_reps"] = int(args.bootstrap_reps)
    if args.fail_on_claim_threshold_failure:
        config.setdefault("failure_policy", {})["fail_on_claim_threshold_failure"] = True

    paths = make_paths(config, args.output_dir)
    rng = np.random.default_rng(int(config["seed"]))

    validation, tables = validate_inputs(config, paths)
    model_validation = validate_required_models(config, tables, paths)

    contract = build_contract_table(paths)
    transition = transition_analysis(config, tables, paths, rng)
    competence = competence_analysis(config, tables, paths, rng)
    fate_recovery = fate_recovery_analysis(config, tables, paths)
    external = external_perturbation_analysis(config, tables, paths, rng)
    ablation = ablation_analysis(config, tables, paths)

    summaries = {
        "virtual_cell_contract_rows": int(len(contract)),
        "transition": transition,
        "competence": competence,
        "fate_recovery": fate_recovery,
        "external": external,
        "ablation": ablation,
    }
    gates = build_claim_gates(config, summaries, tables, paths)
    write_methods_text(config, summaries, paths, config_path)
    manifest = build_manifest(config, config_path, paths, validation, model_validation, gates, summaries, started_utc)

    print(json.dumps({"status": "completed", "manifest": rel(paths.manifest), "claim_gates_pass": manifest["validation_summary"]["claim_gates_pass"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
