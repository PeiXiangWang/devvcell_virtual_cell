from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import pandas as pd

from devguard.io import ensure_dir, read_h5ad, write_dataframe, write_manifest
from devguard.markers import MARKER_MODULES, assign_top_module, score_marker_modules


DEFAULT_INPUT = "data/processed/devguard/GSE123187_tomo_3files.h5ad"


def map_gse123187_spatial_tomo_lineages(
    input_h5ad: str | Path = DEFAULT_INPUT,
    output_dir: str | Path = "results/devguard_real/GSE123187_tomo_spatial_lineage_mapping",
    *,
    mode_filter: str = "tomo_seq",
) -> dict[str, Path]:
    output = ensure_dir(output_dir)
    adata = read_h5ad(input_h5ad)
    if mode_filter and "gse123187_mode" in adata.obs.columns:
        mask = adata.obs["gse123187_mode"].astype(str).eq(mode_filter).to_numpy()
        adata = adata[mask].copy()
    scores, module_gene_table = score_marker_modules(adata, MARKER_MODULES)
    obs_cols = [
        "cell_id",
        "sample_id",
        "time_point",
        "batch",
        "gse123187_mode",
        "source_member",
        "tomo_position",
        "tomo_axis_fraction",
        "tomo_axis_bin",
    ]
    obs = adata.obs[[column for column in obs_cols if column in adata.obs.columns]].copy()
    obs.index = adata.obs_names.astype(str)
    segment_scores = pd.concat([obs, scores], axis=1)
    segment_scores["assigned_marker_module"] = assign_top_module(scores).to_numpy()
    module_cols = list(MARKER_MODULES)

    group_cols = [column for column in ["source_member", "sample_id", "tomo_axis_bin", "assigned_marker_module"] if column in segment_scores.columns]
    axis_counts = segment_scores.groupby(group_cols, dropna=False, observed=True).size().reset_index(name="n_segments")
    total_cols = [column for column in ["source_member", "sample_id", "tomo_axis_bin"] if column in segment_scores.columns]
    axis_counts["fraction_within_axis_bin"] = axis_counts["n_segments"] / axis_counts.groupby(total_cols, dropna=False, observed=True)[
        "n_segments"
    ].transform("sum")
    mean_scores = (
        segment_scores.groupby(total_cols, dropna=False, observed=True)[module_cols]
        .mean()
        .reset_index()
    )
    axis_summary = axis_counts.merge(mean_scores, on=total_cols, how="left").sort_values(
        total_cols + ["fraction_within_axis_bin"], ascending=[True] * len(total_cols) + [False]
    )
    source_summary = (
        segment_scores.groupby(["source_member", "assigned_marker_module"], dropna=False, observed=True)
        .size()
        .reset_index(name="n_segments")
    )
    source_summary["fraction_within_source"] = source_summary["n_segments"] / source_summary.groupby("source_member", dropna=False, observed=True)[
        "n_segments"
    ].transform("sum")
    outputs: dict[str, Path] = {}
    outputs["module_gene_table"] = write_dataframe(module_gene_table, output / "gse123187_marker_module_gene_availability.csv")
    outputs["segment_scores"] = write_dataframe(segment_scores, output / "tomo_segment_marker_module_scores.csv")
    outputs["axis_summary"] = write_dataframe(axis_summary, output / "tomo_axis_lineage_mapping.csv")
    outputs["source_summary"] = write_dataframe(
        source_summary.sort_values(["source_member", "fraction_within_source"], ascending=[True, False]),
        output / "tomo_source_lineage_summary.csv",
    )
    write_manifest(
        output / "gse123187_tomo_lineage_mapping_manifest.json",
        name="map_gse123187_spatial_tomo_lineages",
        inputs=[str(input_h5ad)],
        outputs=[str(path) for path in outputs.values()],
        parameters={"mode_filter": mode_filter, "modules": MARKER_MODULES},
        metrics={
            "n_segments": int(segment_scores.shape[0]),
            "n_sources": int(segment_scores["source_member"].nunique(dropna=False)) if "source_member" in segment_scores.columns else 0,
            "n_modules": len(module_cols),
        },
    )
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description="Map GSE123187 tomo-seq segments to marker-module lineages.")
    parser.add_argument("--input-h5ad", default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", default="results/devguard_real/GSE123187_tomo_spatial_lineage_mapping")
    parser.add_argument("--mode-filter", default="tomo_seq")
    args = parser.parse_args()
    map_gse123187_spatial_tomo_lineages(args.input_h5ad, args.output_dir, mode_filter=args.mode_filter)


if __name__ == "__main__":
    main()
