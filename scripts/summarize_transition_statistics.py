"""Summarize transition benchmark statistics with paired bootstrap intervals."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from devvcell.io import resolve_project_path, write_json  # noqa: E402


plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial Unicode MS", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", default="results/cell_level_v1/tables/cell_level_transition_metrics.csv")
    parser.add_argument("--output-dir", default="results/cell_level_v1")
    parser.add_argument("--reference-model", default="context_residual_mlp")
    parser.add_argument("--bootstrap-reps", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=123)
    return parser.parse_args()


def bootstrap_ci(values: np.ndarray, rng: np.random.Generator, reps: int) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    if len(values) == 0:
        return float("nan"), float("nan")
    samples = rng.choice(values, size=(reps, len(values)), replace=True).mean(axis=1)
    lo, hi = np.percentile(samples, [2.5, 97.5])
    return float(lo), float(hi)


def summarize_models(metrics: pd.DataFrame, rng: np.random.Generator, reps: int) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for model, sub in metrics.groupby("model"):
        row: dict[str, object] = {"model": model, "n_strata": int(len(sub))}
        for metric in ["pair_latent_mse", "centroid_latent_mse", "rbf_mmd"]:
            values = pd.to_numeric(sub[metric], errors="coerce").dropna().to_numpy(dtype=float)
            lo, hi = bootstrap_ci(values, rng, reps)
            row[f"mean_{metric}"] = float(values.mean()) if len(values) else float("nan")
            row[f"ci95_low_{metric}"] = lo
            row[f"ci95_high_{metric}"] = hi
        rows.append(row)
    return pd.DataFrame(rows).sort_values("mean_pair_latent_mse")


def paired_differences(
    metrics: pd.DataFrame,
    reference_model: str,
    rng: np.random.Generator,
    reps: int,
) -> pd.DataFrame:
    key_cols = ["system", "src_stage", "tgt_stage"]
    if reference_model not in set(metrics["model"]):
        reference_model = str(metrics.groupby("model")["pair_latent_mse"].mean().idxmin())

    pivot = metrics.pivot_table(index=key_cols, columns="model", values="pair_latent_mse", aggfunc="mean")
    if reference_model not in pivot.columns:
        raise ValueError(f"Reference model {reference_model!r} is not available.")

    rows: list[dict[str, object]] = []
    for model in sorted(c for c in pivot.columns if c != reference_model):
        paired = pivot[[reference_model, model]].dropna()
        competitor = paired[model].to_numpy(dtype=float)
        reference = paired[reference_model].to_numpy(dtype=float)
        improvement = competitor - reference
        rel_improvement = improvement / np.maximum(competitor, 1e-12)
        lo, hi = bootstrap_ci(improvement, rng, reps)
        rel_lo, rel_hi = bootstrap_ci(rel_improvement, rng, reps)
        rows.append(
            {
                "reference_model": reference_model,
                "competitor_model": model,
                "n_paired_strata": int(len(paired)),
                "mean_competitor_minus_reference_mse": float(improvement.mean()),
                "ci95_low_competitor_minus_reference_mse": lo,
                "ci95_high_competitor_minus_reference_mse": hi,
                "mean_relative_improvement": float(rel_improvement.mean()),
                "ci95_low_relative_improvement": rel_lo,
                "ci95_high_relative_improvement": rel_hi,
                "strata_reference_better": int((improvement > 0).sum()),
            }
        )
    return pd.DataFrame(rows).sort_values("mean_competitor_minus_reference_mse", ascending=False)


def plot_summary(summary: pd.DataFrame, figures_dir: Path) -> None:
    ordered = summary.sort_values("mean_pair_latent_mse")
    yerr = np.vstack(
        [
            ordered["mean_pair_latent_mse"] - ordered["ci95_low_pair_latent_mse"],
            ordered["ci95_high_pair_latent_mse"] - ordered["mean_pair_latent_mse"],
        ]
    )
    fig, ax = plt.subplots(figsize=(7.8, 4.6))
    ax.bar(
        ordered["model"],
        ordered["mean_pair_latent_mse"],
        yerr=yerr,
        capsize=4,
        color=["#2f6f73", "#4e79a7", "#b07aa1", "#6c757d", "#59a14f"][: len(ordered)],
    )
    ax.set_ylabel("Heldout pseudo-pair latent MSE")
    ax.set_xlabel("transition 模型")
    ax.set_title("细胞级 transition bootstrap 95% CI")
    ax.grid(axis="y", alpha=0.25)
    ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(figures_dir / "transition_bootstrap_ci.png", dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    metrics_path = resolve_project_path(args.metrics)
    out_dir = resolve_project_path(args.output_dir)
    tables_dir = out_dir / "tables"
    figures_dir = out_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    metrics = pd.read_csv(metrics_path)
    summary = summarize_models(metrics, rng, args.bootstrap_reps)
    differences = paired_differences(metrics, args.reference_model, rng, args.bootstrap_reps)
    summary.to_csv(tables_dir / "transition_statistical_summary.csv", index=False)
    differences.to_csv(tables_dir / "transition_paired_bootstrap_differences.csv", index=False)
    plot_summary(summary, figures_dir)

    best_model = str(summary.iloc[0]["model"])
    manifest = {
        "analysis": "cell_level_transition_statistical_summary",
        "source_metrics": str(metrics_path.relative_to(PROJECT_ROOT)).replace("\\", "/"),
        "bootstrap_reps": int(args.bootstrap_reps),
        "best_model_by_pair_latent_mse": best_model,
        "model_summary": summary.to_dict(orient="records"),
        "paired_differences": differences.to_dict(orient="records"),
        "outputs": {
            "summary": "results/cell_level_v1/tables/transition_statistical_summary.csv",
            "paired_differences": "results/cell_level_v1/tables/transition_paired_bootstrap_differences.csv",
            "figure": "results/cell_level_v1/figures/transition_bootstrap_ci.png",
        },
    }
    write_json(tables_dir / "transition_statistical_summary.json", manifest)
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
