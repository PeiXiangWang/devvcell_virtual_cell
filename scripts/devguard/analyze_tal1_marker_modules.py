from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import anndata as ad
import pandas as pd

from devguard.io import ensure_dir, write_dataframe, write_manifest
from devguard.markers import MARKER_MODULES, log_normalized_expression, score_marker_modules


DEFAULT_TAL1_H5AD = "data/processed/devguard/MouseGastrulationData_tal1_chimera_full.h5ad"
DEFAULT_CONTROL_H5AD = "data/processed/devguard/MouseGastrulationData_integrated_chimera_controls.h5ad"
DEFAULT_TAL1_CLASSES = "results/devguard_real/MouseGastrulationData_tal1_chimera_full_integrated_e85_strict/perturbation_classification/cell_normality_classes.csv"
DEFAULT_CONTROL_CLASSES = "results/devguard_real/MouseGastrulationData_integrated_chimera_controls_e85_strict/heldout_control_classification/heldout_control_normality_classes.csv"


def _subset_by_cell_ids(adata, cell_ids: pd.Series):
    ids = cell_ids.astype(str).tolist()
    if "cell_id" not in adata.obs.columns:
        raise ValueError("AnnData obs is missing cell_id.")
    obs = adata.obs.copy()
    obs["_devguard_cell_id"] = obs["cell_id"].astype(str)
    if not obs["_devguard_cell_id"].is_unique:
        raise ValueError("cell_id values are not unique; cannot align classification table.")
    adata.obs["_devguard_cell_id"] = obs["_devguard_cell_id"]
    adata.obs_names = adata.obs["_devguard_cell_id"].astype(str)
    missing = sorted(set(ids) - set(adata.obs_names.astype(str)))
    if missing:
        preview = ", ".join(missing[:5])
        raise ValueError(f"Missing {len(missing)} classified cells in AnnData: {preview}")
    return adata[ids].copy()


def _load_scored_cohort(h5ad_path: str | Path, classes_csv: str | Path, cohort: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    classes = pd.read_csv(classes_csv, low_memory=False)
    adata = ad.read_h5ad(h5ad_path)
    subset = _subset_by_cell_ids(adata, classes["cell_id"])
    scores, module_gene_table = score_marker_modules(subset, MARKER_MODULES)
    marker_genes = sorted({gene for genes in MARKER_MODULES.values() for gene in genes})
    marker_expression, _ = log_normalized_expression(subset, marker_genes)

    metadata_cols = [
        "cell_id",
        "dataset_id",
        "sample_id",
        "time_point",
        "lineage",
        "cell_type",
        "perturbation_name",
        "normality_class",
        "assigned_reference_group",
        "reference_other_lineage",
    ]
    metadata = classes[[column for column in metadata_cols if column in classes.columns]].copy()
    metadata.index = subset.obs_names.astype(str)
    scored = pd.concat([metadata, scores], axis=1)
    scored["cohort"] = cohort
    marker_expression = pd.concat([metadata, marker_expression], axis=1)
    marker_expression["cohort"] = cohort
    module_gene_table["cohort"] = cohort
    return scored, marker_expression, module_gene_table


def _summarize_scores(frame: pd.DataFrame, group_cols: list[str], module_cols: list[str]) -> pd.DataFrame:
    return (
        frame.groupby(group_cols, dropna=False, observed=True)
        .agg(n_cells=("cell_id", "size"), **{f"{module}_mean": (module, "mean") for module in module_cols})
        .reset_index()
    )


def _module_effects(frame: pd.DataFrame, module_cols: list[str]) -> pd.DataFrame:
    baseline = frame[
        frame["cohort"].eq("Control") & frame["normality_class"].astype(str).eq("within_stage_normal")
    ]
    baseline_means = baseline[module_cols].mean()
    rows = []
    for (cohort, normality_class), group in frame.groupby(["cohort", "normality_class"], dropna=False, observed=True):
        for module in module_cols:
            rows.append(
                {
                    "cohort": cohort,
                    "normality_class": normality_class,
                    "module": module,
                    "n_cells": int(group.shape[0]),
                    "mean_score": float(group[module].mean()),
                    "control_within_stage_mean": float(baseline_means[module]),
                    "delta_vs_control_within_stage": float(group[module].mean() - baseline_means[module]),
                }
            )
    return pd.DataFrame(rows).sort_values(["cohort", "normality_class", "delta_vs_control_within_stage"], ascending=[True, True, False])


def analyze_tal1_marker_modules(
    tal1_h5ad: str | Path = DEFAULT_TAL1_H5AD,
    control_h5ad: str | Path = DEFAULT_CONTROL_H5AD,
    tal1_classes: str | Path = DEFAULT_TAL1_CLASSES,
    control_classes: str | Path = DEFAULT_CONTROL_CLASSES,
    output_dir: str | Path = "results/devguard_real/MouseGastrulationData_tal1_chimera_full_integrated_e85_strict/marker_modules",
) -> dict[str, Path]:
    output = ensure_dir(output_dir)
    tal1_scores, tal1_markers, tal1_gene_table = _load_scored_cohort(tal1_h5ad, tal1_classes, "Tal1_chimera")
    control_scores, control_markers, control_gene_table = _load_scored_cohort(control_h5ad, control_classes, "Control")
    scores = pd.concat([tal1_scores, control_scores], axis=0, ignore_index=True, sort=False)
    marker_expression = pd.concat([tal1_markers, control_markers], axis=0, ignore_index=True, sort=False)
    module_gene_table = pd.concat([tal1_gene_table, control_gene_table], axis=0, ignore_index=True, sort=False)
    module_cols = list(MARKER_MODULES)
    marker_cols = sorted({gene for genes in MARKER_MODULES.values() for gene in genes if gene in marker_expression.columns})

    outputs: dict[str, Path] = {}
    outputs["module_gene_table"] = write_dataframe(module_gene_table, output / "marker_module_gene_availability.csv")
    outputs["per_cell_module_scores"] = write_dataframe(scores, output / "cell_marker_module_scores.csv")
    outputs["module_by_class"] = write_dataframe(
        _summarize_scores(scores, ["cohort", "normality_class"], module_cols),
        output / "module_scores_by_class.csv",
    )
    outputs["module_by_lineage_class"] = write_dataframe(
        _summarize_scores(scores, ["cohort", "lineage", "normality_class"], module_cols),
        output / "module_scores_by_lineage_class.csv",
    )
    outputs["module_effects"] = write_dataframe(
        _module_effects(scores, module_cols),
        output / "module_score_effects_vs_control_within_stage.csv",
    )
    marker_summary = (
        marker_expression.groupby(["cohort", "normality_class"], dropna=False, observed=True)
        .agg(n_cells=("cell_id", "size"), **{f"{gene}_mean": (gene, "mean") for gene in marker_cols})
        .reset_index()
    )
    outputs["marker_expression"] = write_dataframe(marker_summary, output / "marker_expression_by_class.csv")
    write_manifest(
        output / "tal1_marker_module_manifest.json",
        name="analyze_tal1_marker_modules",
        inputs=[str(tal1_h5ad), str(control_h5ad), str(tal1_classes), str(control_classes)],
        outputs=[str(path) for path in outputs.values()],
        parameters={"modules": MARKER_MODULES},
        metrics={
            "n_tal1_cells": int(tal1_scores.shape[0]),
            "n_control_cells": int(control_scores.shape[0]),
            "n_modules": len(module_cols),
            "n_marker_genes_available": len(marker_cols),
        },
    )
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description="Score Tal1 chimera cells by lineage marker modules.")
    parser.add_argument("--tal1-h5ad", default=DEFAULT_TAL1_H5AD)
    parser.add_argument("--control-h5ad", default=DEFAULT_CONTROL_H5AD)
    parser.add_argument("--tal1-classes", default=DEFAULT_TAL1_CLASSES)
    parser.add_argument("--control-classes", default=DEFAULT_CONTROL_CLASSES)
    parser.add_argument("--output-dir", default="results/devguard_real/MouseGastrulationData_tal1_chimera_full_integrated_e85_strict/marker_modules")
    args = parser.parse_args()
    analyze_tal1_marker_modules(
        args.tal1_h5ad,
        args.control_h5ad,
        args.tal1_classes,
        args.control_classes,
        args.output_dir,
    )


if __name__ == "__main__":
    main()
