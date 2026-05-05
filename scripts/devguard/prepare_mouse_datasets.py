from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401

from devguard.fixtures import create_quick_fixture
from devguard.io import dataset_registry_frame, ensure_dir, load_json, write_dataframe, write_h5ad, write_manifest
from devguard.preprocessing import metadata_summary, standardize_obs, validate_obs_schema


def prepare_quick_dataset(config: dict, output_path: str | Path | None = None) -> Path:
    registry = dataset_registry_frame(config)
    row = registry.loc[registry["dataset_id"].eq("devguard_quick_mouse")].iloc[0]
    output = Path(output_path or row["local_path"])
    metadata_output = Path(row["metadata_output"])
    adata = create_quick_fixture()
    adata = standardize_obs(adata, dataset_id="devguard_quick_mouse")
    missing = validate_obs_schema(adata.obs)
    if missing:
        raise ValueError(f"Quick fixture missing DevGuard obs columns: {missing}")
    write_h5ad(adata, output)
    write_dataframe(metadata_summary(adata.obs), metadata_output)
    write_manifest(
        "results/devguard/dataset_metadata/devguard_quick_mouse_manifest.json",
        name="prepare_quick_dataset",
        inputs=["generated_quick_fixture"],
        outputs=[str(output), str(metadata_output)],
        parameters={"software_validation_only": True},
        metrics={"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars)},
    )
    return output


def prepare_registry(config_path: str | Path, quick_fixture: bool = False) -> None:
    config = load_json(config_path)
    ensure_dir("results/devguard/dataset_metadata")
    registry = dataset_registry_frame(config)
    write_dataframe(registry, "results/devguard/dataset_metadata/dataset_registry.csv")
    if quick_fixture:
        prepare_quick_dataset(config)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare mouse datasets for DevGuard.")
    parser.add_argument("--config", default="config/devguard/datasets_mouse.json")
    parser.add_argument("--quick-fixture", action="store_true", help="Generate the local software validation fixture.")
    args = parser.parse_args()
    prepare_registry(args.config, quick_fixture=args.quick_fixture)


if __name__ == "__main__":
    main()
