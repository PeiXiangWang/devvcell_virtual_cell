from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
import pandas as pd

from devguard.io import ensure_dir, write_dataframe, write_manifest
from devguard.robustness import (
    balanced_downsample_class_fractions,
    bootstrap_sample_class_fractions,
    class_fraction_table,
    sample_class_summary,
    summarize_downsample_iterations,
)


DEFAULT_TAL1 = "results/devguard_real/MouseGastrulationData_tal1_chimera_full_integrated_e85_strict/perturbation_classification/cell_normality_classes.csv"
DEFAULT_T = "results/devguard_real/MouseGastrulationData_t_chimera_full_integrated_e85_strict/perturbation_classification/cell_normality_classes.csv"
DEFAULT_CONTROL = "results/devguard_real/MouseGastrulationData_integrated_chimera_controls_e85_strict/heldout_control_classification/heldout_control_normality_classes.csv"
DEFAULT_QUALITY = "results/devguard_real/MouseGastrulationData_integrated_chimera_controls_e85_strict/normality_reference/normality_reference_quality.csv"


def _load_classification(path: str | Path, cohort: str) -> pd.DataFrame:
    frame = pd.read_csv(path, low_memory=False)
    frame["cohort"] = cohort
    if "perturbation_name" not in frame.columns:
        frame["perturbation_name"] = cohort
    frame["perturbation_name"] = frame["perturbation_name"].astype(str).replace({"control": "Control"})
    if cohort == "Control":
        frame["perturbation_name"] = "Control"
    return frame


def _load_quality(path: str | Path, score_method: str) -> pd.DataFrame:
    quality = pd.read_csv(path)
    if "score_method" in quality.columns:
        quality = quality[quality["score_method"].astype(str).eq(score_method)].copy()
    for column in ["calibration_excess_fpr", "high_fpr_flag"]:
        if column not in quality.columns:
            if column == "calibration_excess_fpr":
                quality[column] = quality["heldout_control_fpr"] - quality["alpha"]
            else:
                quality[column] = quality["heldout_control_fpr"] > quality["alpha"] + 0.05
    return quality


def _attach_quality_flags(frame: pd.DataFrame, quality: pd.DataFrame, alpha: float) -> pd.DataFrame:
    out = frame.copy()
    quality_by_group = quality.set_index("reference_group")
    fpr_map = quality_by_group["heldout_control_fpr"].to_dict()
    high_map = quality_by_group["high_fpr_flag"].astype(bool).to_dict()
    lineage_map = quality_by_group["lineage"].to_dict()
    for column in ["reference_current_same", "reference_other_lineage", "assigned_reference_group"]:
        fpr_col = f"{column}_heldout_fpr"
        high_col = f"{column}_high_fpr"
        lineage_col = f"{column}_lineage"
        out[fpr_col] = out[column].map(fpr_map) if column in out.columns else np.nan
        out[high_col] = out[column].map(high_map).astype("boolean").fillna(False).astype(bool) if column in out.columns else False
        out[lineage_col] = out[column].map(lineage_map) if column in out.columns else ""
    current = "reference_current_same_heldout_fpr"
    out["current_reference_missing"] = out[current].isna()
    out["current_reference_high_or_missing"] = out["reference_current_same_high_fpr"] | out["current_reference_missing"]
    fpr = out[current].fillna(np.inf).astype(float)
    out["current_reference_calibration_weight"] = np.where(
        np.isfinite(fpr) & (fpr > 0),
        np.minimum(1.0, alpha / fpr),
        0.0,
    )
    out.loc[fpr <= alpha, "current_reference_calibration_weight"] = 1.0
    return out


def _quality_table(quality: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "reference_group",
        "time_point",
        "lineage",
        "score_method",
        "alpha",
        "n_train",
        "n_calibration",
        "n_test",
        "n_train_units",
        "n_calibration_units",
        "n_test_units",
        "heldout_control_fpr",
        "calibration_excess_fpr",
        "high_fpr_flag",
    ]
    return quality[[column for column in columns if column in quality.columns]].sort_values(
        ["high_fpr_flag", "heldout_control_fpr"], ascending=[False, False]
    )


def _sensitivity_tables(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    scenarios = []
    scenario_frames = {
        "full": (frame, None),
        "drop_current_high_or_missing_fpr_lineages": (frame[~frame["current_reference_high_or_missing"]].copy(), None),
        "drop_assigned_high_fpr_reference": (frame[~frame["assigned_reference_group_high_fpr"]].copy(), None),
        "downweight_current_reference_fpr": (frame, "current_reference_calibration_weight"),
    }
    for scenario, (scenario_frame, weight_col) in scenario_frames.items():
        if scenario_frame.empty:
            continue
        summary = class_fraction_table(
            scenario_frame,
            ["perturbation_name"],
            weight_col=weight_col,
        )
        summary["scenario"] = scenario
        summary["n_input_cells"] = int(scenario_frame.shape[0])
        scenarios.append(summary)
    long = pd.concat(scenarios, axis=0, ignore_index=True) if scenarios else pd.DataFrame()
    pivot = long.pivot_table(
        index=["scenario", "perturbation_name"],
        columns="normality_class",
        values="fraction",
        fill_value=0,
        aggfunc="sum",
    ).reset_index()
    return long, pivot


def _fate_deviation_reference_summary(tal1: pd.DataFrame, quality: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    quality_by_group = quality.set_index("reference_group")
    target_lineage = quality_by_group["lineage"].to_dict()
    target_fpr = quality_by_group["heldout_control_fpr"].to_dict()
    target_high = quality_by_group["high_fpr_flag"].astype(bool).to_dict()
    fate = tal1[tal1["normality_class"].astype(str).eq("fate_deviation")].copy()
    if fate.empty:
        return pd.DataFrame(), pd.DataFrame()
    fate["target_reference_group"] = fate["reference_other_lineage"].fillna(fate["assigned_reference_group"]).astype(str)
    fate["target_reference_lineage"] = fate["target_reference_group"].map(target_lineage).fillna("")
    fate["target_reference_heldout_fpr"] = fate["target_reference_group"].map(target_fpr)
    fate["target_reference_high_fpr"] = fate["target_reference_group"].map(target_high).fillna(False).astype(bool)
    grouped = (
        fate.groupby(["lineage", "target_reference_lineage", "target_reference_group"], dropna=False, observed=True)
        .agg(
            n_cells=("cell_id", "size"),
            mean_p_other_lineage=("p_other_lineage", "mean"),
            median_p_other_lineage=("p_other_lineage", "median"),
            mean_p_current_same=("p_current_same", "mean"),
            median_pvalue_margin=("pvalue_margin", "median"),
            n_samples=("sample_id", "nunique"),
            target_reference_heldout_fpr=("target_reference_heldout_fpr", "first"),
            target_reference_high_fpr=("target_reference_high_fpr", "first"),
        )
        .reset_index()
    )
    grouped["fraction_of_all_tal1_fate_deviation"] = grouped["n_cells"] / max(int(fate.shape[0]), 1)
    grouped["fraction_within_observed_lineage"] = grouped["n_cells"] / grouped.groupby("lineage", dropna=False)["n_cells"].transform("sum")
    grouped = grouped.sort_values(["n_cells", "fraction_within_observed_lineage"], ascending=[False, False]).reset_index(drop=True)

    overall = (
        fate.groupby(["target_reference_lineage", "target_reference_group"], dropna=False, observed=True)
        .agg(
            n_cells=("cell_id", "size"),
            n_observed_lineages=("lineage", "nunique"),
            n_samples=("sample_id", "nunique"),
            mean_p_other_lineage=("p_other_lineage", "mean"),
            target_reference_heldout_fpr=("target_reference_heldout_fpr", "first"),
            target_reference_high_fpr=("target_reference_high_fpr", "first"),
        )
        .reset_index()
    )
    overall["fraction_of_all_tal1_fate_deviation"] = overall["n_cells"] / max(int(fate.shape[0]), 1)
    overall = overall.sort_values("n_cells", ascending=False).reset_index(drop=True)
    return grouped, overall


def analyze_chimera_robustness(
    tal1_csv: str | Path = DEFAULT_TAL1,
    t_csv: str | Path = DEFAULT_T,
    control_csv: str | Path = DEFAULT_CONTROL,
    quality_csv: str | Path = DEFAULT_QUALITY,
    output_dir: str | Path = "results/devguard_real/chimera_robustness",
    *,
    score_method: str = "knn_distance",
    alpha: float = 0.05,
    n_bootstrap: int = 1000,
    n_downsample_iterations: int = 200,
    max_cells_per_stratum: int = 200,
    seed: int = 42,
) -> dict[str, Path]:
    output = ensure_dir(output_dir)
    tal1 = _load_classification(tal1_csv, "Tal1")
    t_frame = _load_classification(t_csv, "T")
    control = _load_classification(control_csv, "Control")
    combined = pd.concat([tal1, t_frame, control], axis=0, ignore_index=True, sort=False)
    combined["perturbation_name"] = combined["perturbation_name"].replace(
        {"Tal1_chimera": "Tal1_chimera", "T_chimera": "T_chimera", "Control": "Control"}
    )
    quality = _load_quality(quality_csv, score_method)
    flagged = _attach_quality_flags(combined, quality, alpha)

    outputs: dict[str, Path] = {}
    outputs["quality_flags"] = write_dataframe(_quality_table(quality), output / "reference_lineage_fpr_flags.csv")
    outputs["sample_class_summary"] = write_dataframe(sample_class_summary(flagged), output / "sample_class_summary.csv")
    outputs["sample_bootstrap"] = write_dataframe(
        bootstrap_sample_class_fractions(flagged, n_bootstrap=n_bootstrap, seed=seed),
        output / "sample_level_bootstrap_class_fractions.csv",
    )
    sensitivity_long, sensitivity_pivot = _sensitivity_tables(flagged)
    outputs["fpr_sensitivity_long"] = write_dataframe(sensitivity_long, output / "high_fpr_lineage_sensitivity_long.csv")
    outputs["fpr_sensitivity_pivot"] = write_dataframe(sensitivity_pivot, output / "high_fpr_lineage_sensitivity_pivot.csv")

    downsample_frames = []
    global_iterations = balanced_downsample_class_fractions(
        flagged,
        strata_cols=[],
        n_iterations=n_downsample_iterations,
        max_cells_per_stratum=None,
        seed=seed,
    )
    global_iterations["scenario"] = "equal_total_cells"
    downsample_frames.append(global_iterations)
    lineage_iterations = balanced_downsample_class_fractions(
        flagged,
        strata_cols=["lineage"],
        n_iterations=n_downsample_iterations,
        max_cells_per_stratum=max_cells_per_stratum,
        seed=seed + 1,
    )
    lineage_iterations["scenario"] = "matched_observed_lineage_counts"
    downsample_frames.append(lineage_iterations)
    downsample_iterations = pd.concat(downsample_frames, axis=0, ignore_index=True)
    outputs["balanced_iterations"] = write_dataframe(
        downsample_iterations,
        output / "balanced_downsample_iterations.csv",
    )
    summary_parts = []
    for scenario, part in downsample_iterations.groupby("scenario", dropna=False, observed=True):
        summary = summarize_downsample_iterations(part)
        summary["scenario"] = scenario
        summary_parts.append(summary)
    downsample_summary = pd.concat(summary_parts, axis=0, ignore_index=True) if summary_parts else pd.DataFrame()
    outputs["balanced_summary"] = write_dataframe(
        downsample_summary,
        output / "balanced_downsample_class_fraction_summary.csv",
    )

    fate_routes, fate_overall = _fate_deviation_reference_summary(tal1, quality)
    outputs["tal1_fate_routes"] = write_dataframe(fate_routes, output / "tal1_fate_deviation_reference_routes.csv")
    outputs["tal1_fate_targets"] = write_dataframe(fate_overall, output / "tal1_fate_deviation_reference_targets.csv")
    outputs["flagged_cells"] = write_dataframe(
        flagged[
            [
                column
                for column in [
                    "cohort",
                    "cell_id",
                    "perturbation_name",
                    "sample_id",
                    "time_point",
                    "lineage",
                    "normality_class",
                    "reference_current_same",
                    "reference_current_same_heldout_fpr",
                    "reference_current_same_high_fpr",
                    "current_reference_missing",
                    "assigned_reference_group",
                    "assigned_reference_group_heldout_fpr",
                    "assigned_reference_group_high_fpr",
                    "current_reference_calibration_weight",
                ]
                if column in flagged.columns
            ]
        ],
        output / "cell_level_calibration_flags.csv",
    )

    write_manifest(
        output / "chimera_robustness_manifest.json",
        name="analyze_chimera_robustness",
        inputs=[str(tal1_csv), str(t_csv), str(control_csv), str(quality_csv)],
        outputs=[str(path) for path in outputs.values()],
        parameters={
            "score_method": score_method,
            "alpha": alpha,
            "n_bootstrap": n_bootstrap,
            "n_downsample_iterations": n_downsample_iterations,
            "max_cells_per_stratum": max_cells_per_stratum,
            "seed": seed,
        },
        metrics={
            "n_tal1_cells": int(tal1.shape[0]),
            "n_t_cells": int(t_frame.shape[0]),
            "n_control_cells": int(control.shape[0]),
            "n_high_fpr_reference_groups": int(quality["high_fpr_flag"].sum()),
        },
    )
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze Tal1/T chimera robustness and calibration sensitivity.")
    parser.add_argument("--tal1-csv", default=DEFAULT_TAL1)
    parser.add_argument("--t-csv", default=DEFAULT_T)
    parser.add_argument("--control-csv", default=DEFAULT_CONTROL)
    parser.add_argument("--quality-csv", default=DEFAULT_QUALITY)
    parser.add_argument("--output-dir", default="results/devguard_real/chimera_robustness")
    parser.add_argument("--score-method", default="knn_distance")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    parser.add_argument("--n-downsample-iterations", type=int, default=200)
    parser.add_argument("--max-cells-per-stratum", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    analyze_chimera_robustness(
        args.tal1_csv,
        args.t_csv,
        args.control_csv,
        args.quality_csv,
        args.output_dir,
        score_method=args.score_method,
        alpha=args.alpha,
        n_bootstrap=args.n_bootstrap,
        n_downsample_iterations=args.n_downsample_iterations,
        max_cells_per_stratum=args.max_cells_per_stratum,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
