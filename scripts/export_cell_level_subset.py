"""Export a reproducible cell-level H5AD subset for DevVCell training.

The full source matrix is large. This script samples a balanced subset by
developmental stage and broad biological system, then writes a compact H5AD
that downstream experiments can load repeatedly.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from devvcell.io import load_json, resolve_project_path  # noqa: E402
from devvcell.stages import stage_number  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        default="config/cell_level_baseline.json",
        help="Path to the cell-level experiment config.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite the processed subset if it already exists.",
    )
    return parser.parse_args()


def category_pattern_mask(series: pd.Series, patterns: list[str]) -> pd.Series:
    """Return a boolean mask without expanding categorical values unnecessarily."""

    regex = re.compile("|".join(f"(?:{p})" for p in patterns), flags=re.IGNORECASE)
    if isinstance(series.dtype, pd.CategoricalDtype):
        matching = [cat for cat in series.cat.categories.astype(str) if regex.search(cat)]
        return series.isin(matching)
    return series.astype(str).str.contains(regex, na=False)


def stage_mask(series: pd.Series, stages: set[int]) -> pd.Series:
    if isinstance(series.dtype, pd.CategoricalDtype):
        matching = [cat for cat in series.cat.categories if stage_number(cat) in stages]
        return series.isin(matching)
    return series.map(stage_number).isin(stages)


def sample_indices(
    obs: pd.DataFrame,
    config: dict,
    rng: np.random.Generator,
) -> tuple[np.ndarray, pd.DataFrame]:
    stage_values = set(int(s) for s in config["stages"])
    in_stage_window = stage_mask(obs["development_stage"], stage_values)
    max_cells = int(config["max_cells_per_stage_system"])

    selected: list[int] = []
    manifest_rows: list[dict[str, object]] = []

    for system_name, spec in config["systems"].items():
        system_mask = category_pattern_mask(obs[spec["obs_column"]], spec["patterns"])
        eligible = obs.index[in_stage_window & system_mask]
        if len(eligible) == 0:
            manifest_rows.append(
                {
                    "system": system_name,
                    "stage": "all",
                    "available_cells": 0,
                    "sampled_cells": 0,
                }
            )
            continue

        eligible_obs = obs.loc[eligible, ["development_stage"]].copy()
        eligible_obs["stage_num"] = eligible_obs["development_stage"].map(stage_number)
        for stage_num in sorted(stage_values):
            stage_index = eligible_obs.index[eligible_obs["stage_num"] == stage_num].to_numpy()
            available = int(len(stage_index))
            if available == 0:
                sampled = np.array([], dtype=object)
            elif available <= max_cells:
                sampled = stage_index
            else:
                sampled = rng.choice(stage_index, size=max_cells, replace=False)
            selected.extend(obs.index.get_indexer(sampled).tolist())
            manifest_rows.append(
                {
                    "system": system_name,
                    "stage": f"Theiler stage {stage_num}",
                    "available_cells": available,
                    "sampled_cells": int(len(sampled)),
                }
            )

    selected_array = np.array(sorted(set(i for i in selected if i >= 0)), dtype=np.int64)
    return selected_array, pd.DataFrame(manifest_rows)


def assign_system_labels(obs: pd.DataFrame, config: dict) -> pd.Series:
    labels = pd.Series("unassigned", index=obs.index, dtype="object")
    for system_name, spec in config["systems"].items():
        labels.loc[category_pattern_mask(obs[spec["obs_column"]], spec["patterns"])] = system_name
    return labels


def main() -> None:
    args = parse_args()
    config = load_json(args.config)
    raw_path = resolve_project_path(config["raw_h5ad_path"])
    out_path = resolve_project_path(config["processed_subset_path"])
    manifest_path = out_path.with_suffix(".manifest.csv")

    if out_path.exists() and not args.force:
        raise FileExistsError(f"{out_path} already exists. Pass --force to overwrite it.")

    rng = np.random.default_rng(int(config["seed"]))
    print(f"Opening backed H5AD: {raw_path}")
    adata = ad.read_h5ad(raw_path, backed="r")
    obs = adata.obs.copy()

    selected_indices, manifest = sample_indices(obs, config, rng)
    if len(selected_indices) == 0:
        raise RuntimeError("No cells were selected. Check stages and system patterns.")

    print(f"Loading selected cells into memory: {len(selected_indices):,}")
    subset = adata[selected_indices, :].to_memory()
    subset.obs["devvcell_system"] = assign_system_labels(subset.obs, config).astype("category")
    subset.obs["stage_num"] = subset.obs["development_stage"].map(stage_number).astype(int)
    subset.uns["devvcell_subset_config"] = json.dumps(config, ensure_ascii=False)
    subset.uns["devvcell_source_h5ad"] = str(raw_path)
    subset.uns["devvcell_selected_cells"] = int(subset.n_obs)
    subset.uns["devvcell_selected_genes"] = int(subset.n_vars)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing subset: {out_path}")
    subset.write_h5ad(out_path, compression="gzip")
    manifest.to_csv(manifest_path, index=False)

    adata.file.close()
    print(f"Done. Cells={subset.n_obs:,}; genes={subset.n_vars:,}; manifest={manifest_path}")


if __name__ == "__main__":
    main()
