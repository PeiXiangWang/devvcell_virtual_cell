from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import joblib
import numpy as np
import pandas as pd

from devguard.conformal import conformal_p_values
from devguard.embedding import apply_batch_centering
from devguard.io import ensure_dir, read_h5ad, write_dataframe, write_manifest
from devguard.normality import score_cells
from devguard.preprocessing import standardize_obs


DEFAULT_H5AD = "data/processed/devguard/GSE212050_strict_sample_13285.h5ad"
DEFAULT_MODEL = "results/devguard_real/GSE212050_strict_sample/normality_reference/devguard_normality_model.joblib"
DEFAULT_HELDOUT = "results/devguard_real/GSE212050_strict_sample/heldout_control_classification/heldout_control_normality_classes.csv"


def _class_false_positive_by_sample(heldout_csv: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(heldout_csv, low_memory=False)
    frame["is_false_positive"] = ~frame["normality_class"].astype(str).eq("within_stage_normal")
    frame["is_abnormal_off_normal"] = frame["normality_class"].astype(str).eq("abnormal_off_normal")
    return (
        frame.groupby(["time_point", "lineage", "sample_id"], dropna=False, observed=True)
        .agg(
            n_cells=("cell_id", "size"),
            false_positive_cells=("is_false_positive", "sum"),
            abnormal_off_normal_cells=("is_abnormal_off_normal", "sum"),
            fate_deviation_cells=("normality_class", lambda values: int((values.astype(str) == "fate_deviation").sum())),
            dominant_class=("normality_class", lambda values: values.astype(str).value_counts().idxmax()),
        )
        .reset_index()
        .assign(
            false_positive_fraction=lambda df: df["false_positive_cells"] / df["n_cells"].where(df["n_cells"] > 0, np.nan),
            abnormal_off_normal_fraction=lambda df: df["abnormal_off_normal_cells"] / df["n_cells"].where(df["n_cells"] > 0, np.nan),
            fate_deviation_fraction=lambda df: df["fate_deviation_cells"] / df["n_cells"].where(df["n_cells"] > 0, np.nan),
        )
        .sort_values(["false_positive_fraction", "n_cells"], ascending=[False, False])
        .reset_index(drop=True)
    )


def _organoid_lineage_composition(obs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_cols = ["time_point", "sample_id"]
    comp = obs.groupby(sample_cols + ["lineage"], dropna=False, observed=True).size().reset_index(name="n_cells")
    comp["fraction"] = comp["n_cells"] / comp.groupby(sample_cols, dropna=False, observed=True)["n_cells"].transform("sum")
    pivot = comp.pivot_table(
        index=sample_cols,
        columns="lineage",
        values="fraction",
        fill_value=0,
        aggfunc="sum",
        observed=True,
    ).reset_index()
    for column in ["mesodermal", "neural", "intermediate"]:
        if column not in pivot.columns:
            pivot[column] = 0.0
    lineage_cols = [column for column in pivot.columns if column not in sample_cols]
    pivot["dominant_lineage"] = pivot[lineage_cols].idxmax(axis=1)
    pivot["dominant_lineage_fraction"] = pivot[lineage_cols].max(axis=1)
    pivot["mesodermal_minus_neural"] = pivot["mesodermal"] - pivot["neural"]
    totals = obs.groupby(sample_cols, dropna=False, observed=True).size().reset_index(name="n_cells")
    bias = pivot.merge(totals, on=sample_cols, how="left").sort_values(["time_point", "mesodermal_minus_neural"])
    return comp.sort_values(sample_cols + ["lineage"]).reset_index(drop=True), bias.reset_index(drop=True)


def _transform_control_embeddings(h5ad_path: str | Path, model_path: str | Path):
    model = joblib.load(model_path)
    adata = read_h5ad(h5ad_path)
    dataset_id = str(adata.obs["dataset_id"].iloc[0]) if "dataset_id" in adata.obs else Path(h5ad_path).stem
    adata = standardize_obs(adata, dataset_id=dataset_id)
    control = adata[adata.obs["is_control"].astype(bool).to_numpy()].copy()
    obs = control.obs.copy().reset_index(drop=True)
    embeddings = model["embedding_model"].transform(control)
    embeddings = apply_batch_centering(
        embeddings,
        obs,
        model.get("batch_centering"),
        fallback_columns=model.get("config", {}).get("embedding", {}).get("batch_center_fallback_columns", ["dataset_id"]),
    )
    return obs, embeddings, model


def _leave_one_organoid_out_fpr(
    obs: pd.DataFrame,
    embeddings: np.ndarray,
    model: dict,
    *,
    score_method: str = "knn_distance",
    alpha: float = 0.05,
    min_query_cells: int = 10,
    min_other_units: int = 4,
    seed: int = 42,
) -> pd.DataFrame:
    config = model.get("config", {})
    grouping = config.get("reference_grouping", {})
    time_col = grouping.get("time_column", "time_point")
    lineage_col = grouping.get("lineage_column", "lineage")
    unit_col = config.get("sample_split", {}).get("unit_column", "sample_id")
    k = int(config.get("knn", {}).get("k", 15))
    regularization = float(config.get("mahalanobis", {}).get("regularization", 0.01))
    train_fraction = float(config.get("splits", {}).get("train_fraction", 0.34))
    calibration_fraction = float(config.get("splits", {}).get("calibration_fraction", 0.33))
    train_share = train_fraction / max(train_fraction + calibration_fraction, 1e-12)
    rng = np.random.default_rng(seed)
    rows = []
    for (time_point, lineage), group in obs.groupby([time_col, lineage_col], dropna=False, observed=True):
        group_positions = group.index.to_numpy(dtype=int)
        units_by_position = group[unit_col].astype("string").fillna("NA").astype(str)
        units = np.asarray(sorted(units_by_position.unique()))
        if units.size < min_other_units + 1:
            continue
        for unit in units:
            query_positions = group_positions[units_by_position.eq(unit).to_numpy()]
            if query_positions.size < min_query_cells:
                continue
            other_units = np.asarray([item for item in units if item != unit])
            if other_units.size < min_other_units:
                continue
            shuffled = rng.permutation(other_units)
            n_train_units = max(1, int(round(shuffled.size * train_share)))
            if n_train_units >= shuffled.size:
                n_train_units = shuffled.size - 1
            train_units = set(shuffled[:n_train_units])
            calibration_units = set(shuffled[n_train_units:])
            train_positions = group_positions[units_by_position.isin(train_units).to_numpy()]
            calibration_positions = group_positions[units_by_position.isin(calibration_units).to_numpy()]
            if min(train_positions.size, calibration_positions.size) == 0:
                continue
            calibration_scores = score_cells(
                embeddings[calibration_positions],
                embeddings[train_positions],
                method=score_method,
                k=k,
                regularization=regularization,
            )
            query_scores = score_cells(
                embeddings[query_positions],
                embeddings[train_positions],
                method=score_method,
                k=k,
                regularization=regularization,
            )
            p_values = conformal_p_values(calibration_scores, query_scores)
            rows.append(
                {
                    "time_point": time_point,
                    "lineage": lineage,
                    "sample_id": unit,
                    "score_method": score_method,
                    "alpha": alpha,
                    "n_query_cells": int(query_positions.size),
                    "n_train_cells": int(train_positions.size),
                    "n_calibration_cells": int(calibration_positions.size),
                    "n_other_units": int(other_units.size),
                    "n_train_units": int(len(train_units)),
                    "n_calibration_units": int(len(calibration_units)),
                    "leave_one_organoid_fpr": float((p_values < alpha).mean()),
                    "median_p_value": float(np.median(p_values)),
                    "p10": float(np.quantile(p_values, 0.10)),
                    "p90": float(np.quantile(p_values, 0.90)),
                }
            )
    return pd.DataFrame(rows).sort_values(["leave_one_organoid_fpr", "n_query_cells"], ascending=[False, False]).reset_index(drop=True)


def analyze_gse212050_organoid_heterogeneity(
    input_h5ad: str | Path = DEFAULT_H5AD,
    model_path: str | Path = DEFAULT_MODEL,
    heldout_csv: str | Path = DEFAULT_HELDOUT,
    output_dir: str | Path = "results/devguard_real/GSE212050_strict_sample/organoid_heterogeneity_control",
    *,
    score_method: str = "knn_distance",
    alpha: float = 0.05,
    seed: int = 42,
) -> dict[str, Path]:
    output = ensure_dir(output_dir)
    obs, embeddings, model = _transform_control_embeddings(input_h5ad, model_path)
    composition, bias = _organoid_lineage_composition(obs)
    loo = _leave_one_organoid_out_fpr(
        obs,
        embeddings,
        model,
        score_method=score_method,
        alpha=alpha,
        seed=seed,
    )
    outputs: dict[str, Path] = {}
    outputs["composition"] = write_dataframe(composition, output / "organoid_lineage_composition.csv")
    outputs["bias"] = write_dataframe(bias, output / "organoid_lineage_bias_summary.csv")
    if Path(heldout_csv).exists():
        outputs["heldout_sample_fpr"] = write_dataframe(
            _class_false_positive_by_sample(heldout_csv),
            output / "heldout_false_positive_by_organoid.csv",
        )
    outputs["leave_one_organoid"] = write_dataframe(loo, output / "leave_one_organoid_out_fpr.csv")
    write_manifest(
        output / "gse212050_organoid_heterogeneity_manifest.json",
        name="analyze_gse212050_organoid_heterogeneity",
        inputs=[str(input_h5ad), str(model_path), str(heldout_csv)],
        outputs=[str(path) for path in outputs.values()],
        parameters={"score_method": score_method, "alpha": alpha, "seed": seed},
        metrics={
            "n_cells": int(obs.shape[0]),
            "n_organoid_units": int(obs["sample_id"].nunique(dropna=False)),
            "n_leave_one_organoid_tests": int(loo.shape[0]),
            "mean_leave_one_organoid_fpr": float(loo["leave_one_organoid_fpr"].mean()) if not loo.empty else float("nan"),
            "max_leave_one_organoid_fpr": float(loo["leave_one_organoid_fpr"].max()) if not loo.empty else float("nan"),
        },
    )
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description="Control DevGuard calibration against GSE212050 inter-organoid heterogeneity.")
    parser.add_argument("--input-h5ad", default=DEFAULT_H5AD)
    parser.add_argument("--model-path", default=DEFAULT_MODEL)
    parser.add_argument("--heldout-csv", default=DEFAULT_HELDOUT)
    parser.add_argument("--output-dir", default="results/devguard_real/GSE212050_strict_sample/organoid_heterogeneity_control")
    parser.add_argument("--score-method", default="knn_distance")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    analyze_gse212050_organoid_heterogeneity(
        args.input_h5ad,
        args.model_path,
        args.heldout_csv,
        args.output_dir,
        score_method=args.score_method,
        alpha=args.alpha,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
