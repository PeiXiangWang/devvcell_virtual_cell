from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import anndata as ad

from devguard.io import write_h5ad, write_manifest


def combine_control_h5ads(inputs: list[str | Path], output: str | Path, *, dataset_id: str) -> Path:
    adatas = []
    input_paths = [Path(path) for path in inputs]
    for path in input_paths:
        current = ad.read_h5ad(path)
        source_dataset = str(current.obs["dataset_id"].iloc[0]) if "dataset_id" in current.obs else path.stem
        if "is_control" in current.obs:
            current = current[current.obs["is_control"].astype(bool).to_numpy()].copy()
        current.obs["source_dataset_id"] = source_dataset
        current.obs["sample_id"] = source_dataset + "__" + current.obs["sample_id"].astype(str)
        current.obs["cell_id"] = source_dataset + "__" + current.obs["cell_id"].astype(str)
        current.obs["dataset_id"] = dataset_id
        current.obs_names = current.obs["cell_id"].astype(str)
        adatas.append(current)
    combined = ad.concat(adatas, join="inner", merge="same", label="source_h5ad", keys=[path.stem for path in input_paths])
    combined.obs["dataset_id"] = dataset_id
    combined.obs["condition"] = "control"
    combined.obs["perturbation_name"] = "control"
    combined.obs["perturbation_type"] = "none"
    combined.obs["is_control"] = True
    combined.obs["is_perturbed"] = False
    output_path = write_h5ad(combined, output)
    write_manifest(
        Path("results/devguard/dataset_metadata") / f"{dataset_id}_combine_manifest.json",
        name="combine_control_h5ads",
        inputs=[str(path) for path in input_paths],
        outputs=[str(output_path)],
        parameters={"dataset_id": dataset_id},
        metrics={"n_cells": int(combined.n_obs), "n_genes": int(combined.n_vars), "n_inputs": len(input_paths)},
    )
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Combine control cells from multiple DevGuard H5AD files.")
    parser.add_argument("--input", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--dataset-id", required=True)
    args = parser.parse_args()
    combine_control_h5ads(args.input, args.output, dataset_id=args.dataset_id)


if __name__ == "__main__":
    main()
