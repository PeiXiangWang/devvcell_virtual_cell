from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401

from devspectrum.io import ensure_dir, load_config, read_h5ad, write_dataframe, write_manifest
from devspectrum.timeseries import aggregate_condition_timeseries, build_stage_lineage_module_timeseries, make_quick_timeseries


def build_timeseries_from_config(config_path: str | Path, *, quick: bool = False) -> Path:
    config = load_config(config_path)
    output_root = Path(config.get("output_dir", "results/devspectrum"))
    output_dir = ensure_dir(output_root / "timeseries")
    if quick:
        timeseries = make_quick_timeseries()
        coverage = timeseries[["module_name", "n_genes_in_module", "module_gene_coverage", "feature_type"]].drop_duplicates()
    else:
        input_h5ad = Path(config["input_h5ad"])
        adata = read_h5ad(input_h5ad)
        dataset_id = str(adata.obs["dataset_id"].iloc[0]) if "dataset_id" in adata.obs else input_h5ad.stem
        features = config.get("features", {})
        timeseries, coverage = build_stage_lineage_module_timeseries(
            adata,
            dataset_id=dataset_id,
            module_registry=config.get("modules", {}).get("module_registry"),
            gene_symbol_map=config.get("modules", {}).get("gene_symbol_map"),
            time_column=config.get("time_column", "time_point"),
            time_numeric_column=config.get("time_numeric_column", "time_numeric"),
            lineage_column=config.get("lineage_column", "lineage"),
            sample_column=config.get("sample_column", "sample_id"),
            condition_column=config.get("condition_column", "condition"),
            min_cells_per_group=int(config.get("min_cells_per_group", 50)),
            include_module_scores=bool(features.get("module_scores", True)),
            include_latent_dimensions=bool(features.get("latent_dimensions", False)),
            latent_dim=int(features.get("latent_dim", 5)),
            seed=int(config.get("seed", 42)),
        )
    aggregate = aggregate_condition_timeseries(timeseries)
    timeseries_path = output_dir / "stage_lineage_module_timeseries.csv"
    aggregate_path = output_dir / "condition_stage_lineage_module_timeseries.csv"
    coverage_path = output_dir / "module_gene_coverage.csv"
    write_dataframe(timeseries, timeseries_path)
    write_dataframe(aggregate, aggregate_path)
    write_dataframe(coverage, coverage_path)
    write_manifest(
        output_dir / "timeseries_manifest.json",
        name="build_stage_lineage_module_timeseries",
        inputs=[] if quick else [str(config.get("input_h5ad"))],
        outputs=[str(timeseries_path), str(aggregate_path), str(coverage_path)],
        parameters={**config, "quick": quick},
        metrics={
            "n_rows": int(timeseries.shape[0]),
            "n_condition_rows": int(aggregate.shape[0]),
            "n_lineages": int(timeseries["lineage"].nunique()),
            "n_modules": int(timeseries["module_name"].nunique()),
            "n_time_points": int(timeseries["time_point"].nunique()),
        },
    )
    return timeseries_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Build DevSpectrum stage-lineage-module time-series features.")
    parser.add_argument("--config", default="config/devspectrum/devspectrum_gse212050_control.json")
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    build_timeseries_from_config(args.config, quick=args.quick)


if __name__ == "__main__":
    main()
