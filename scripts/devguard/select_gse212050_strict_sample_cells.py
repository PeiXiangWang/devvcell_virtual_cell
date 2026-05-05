from __future__ import annotations

import argparse
import math
from pathlib import Path

import _bootstrap  # noqa: F401
import pandas as pd

from devguard.io import ensure_dir, write_dataframe, write_manifest


def _infer_gse212050_metadata(obs: pd.DataFrame) -> pd.DataFrame:
    time_point = obs["timepoint.demultiplexed"].fillna(obs["timepoint"]).astype(str)
    cell_type = (
        obs["celltype.mapped.extended"]
        .fillna(obs["celltype.mapped.original"])
        .fillna(obs["cluster"])
        .astype(str)
    )
    lineage = obs["gastr_type"].fillna("").astype(str)
    lineage = lineage.mask(lineage.isin(["", "NA", "nan"]), cell_type)
    sample = obs["sample"].astype(str)
    multi_class = obs["MULTI_class"].fillna("").astype(str)
    sample_id = multi_class.where(multi_class.str.startswith("Bar"), sample)
    return pd.DataFrame(
        {
            "cell_barcode": obs["cell_barcode"].astype(str),
            "time_point": time_point,
            "lineage": lineage,
            "sample_id": sample_id,
        }
    )


def select_strict_cells(
    obs_csv: str | Path,
    output_cells: str | Path,
    *,
    summary_csv: str | Path,
    min_cells_per_group: int = 200,
    min_units_per_group: int = 8,
    max_cells_per_group: int = 1500,
    seed: int = 42,
) -> Path:
    obs = pd.read_csv(obs_csv, low_memory=False)
    metadata = _infer_gse212050_metadata(obs)
    group_summary = (
        metadata.groupby(["time_point", "lineage"], observed=True)
        .agg(n_cells=("cell_barcode", "size"), n_units=("sample_id", "nunique"))
        .reset_index()
    )
    keep_groups = group_summary[
        (group_summary["n_cells"] >= min_cells_per_group) & (group_summary["n_units"] >= min_units_per_group)
    ][["time_point", "lineage"]]

    selected: list[str] = []
    for _, group_row in keep_groups.sort_values(["time_point", "lineage"]).iterrows():
        mask = (metadata["time_point"] == group_row["time_point"]) & (metadata["lineage"] == group_row["lineage"])
        group = metadata.loc[mask]
        units = sorted(group["sample_id"].unique())
        per_unit_cap = max(1, math.ceil(max_cells_per_group / len(units)))
        picks: list[str] = []
        for unit in units:
            unit_cells = group.loc[group["sample_id"] == unit, "cell_barcode"]
            n_take = min(len(unit_cells), per_unit_cap)
            picks.extend(unit_cells.sample(n=n_take, random_state=seed).tolist())
        if len(picks) > max_cells_per_group:
            picks = pd.Series(picks).sample(n=max_cells_per_group, random_state=seed).tolist()
        selected.extend(picks)

    output_path = Path(output_cells)
    ensure_dir(output_path.parent)
    output_path.write_text("\n".join(selected) + "\n", encoding="utf-8")

    selected_metadata = metadata[metadata["cell_barcode"].isin(selected)]
    selected_summary = (
        selected_metadata.groupby(["time_point", "lineage"], observed=True)
        .agg(n_cells=("cell_barcode", "size"), n_units=("sample_id", "nunique"))
        .reset_index()
    )
    summary_path = Path(summary_csv)
    write_dataframe(selected_summary, summary_path)
    write_manifest(
        summary_path.with_suffix(".manifest.json"),
        name="select_gse212050_strict_sample_cells",
        inputs=[str(obs_csv)],
        outputs=[str(output_path), str(summary_path)],
        parameters={
            "min_cells_per_group": min_cells_per_group,
            "min_units_per_group": min_units_per_group,
            "max_cells_per_group": max_cells_per_group,
            "seed": seed,
        },
        metrics={"n_cells": len(selected), "n_groups": int(selected_summary.shape[0])},
    )
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Select GSE212050 cells for strict sample-level DevGuard calibration.")
    parser.add_argument("--obs-csv", default="data/external/GSE212050/components/obs.csv")
    parser.add_argument("--output-cells", default="data/external/GSE212050/gse212050_strict_sample_cells.txt")
    parser.add_argument("--summary-csv", default="results/devguard/dataset_metadata/GSE212050_strict_sample_selection.csv")
    parser.add_argument("--min-cells-per-group", type=int, default=200)
    parser.add_argument("--min-units-per-group", type=int, default=8)
    parser.add_argument("--max-cells-per-group", type=int, default=1500)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    select_strict_cells(
        args.obs_csv,
        args.output_cells,
        summary_csv=args.summary_csv,
        min_cells_per_group=args.min_cells_per_group,
        min_units_per_group=args.min_units_per_group,
        max_cells_per_group=args.max_cells_per_group,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
