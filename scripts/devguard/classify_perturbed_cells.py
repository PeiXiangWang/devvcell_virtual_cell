from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import joblib

from devguard.classification import classify_cells_against_reference, summarize_classes
from devguard.io import ensure_dir, load_json, read_h5ad, write_dataframe, write_manifest
from devguard.preprocessing import standardize_obs


def classify_perturbed_cells(config_path: str | Path) -> Path:
    config = load_json(config_path)
    input_h5ad = Path(config["input_h5ad"])
    model_path = Path(config["model_path"])
    output_dir = ensure_dir(config.get("output_dir", "results/devguard/perturbation_classification"))
    model = joblib.load(model_path)
    adata = read_h5ad(input_h5ad)
    dataset_id = str(adata.obs["dataset_id"].iloc[0]) if "dataset_id" in adata.obs else input_h5ad.stem
    adata = standardize_obs(adata, dataset_id=dataset_id)
    perturbed = adata[adata.obs["is_perturbed"].astype(bool).to_numpy()].copy()
    if perturbed.n_obs == 0:
        raise ValueError("No perturbed cells found for classification.")
    embeddings = model["embedding_model"].transform(perturbed)
    obs = perturbed.obs.copy().reset_index(drop=True)
    score_method = config.get("score_method", "knn_distance")
    results = classify_cells_against_reference(
        embeddings,
        obs,
        model["groups"],
        score_method=score_method,
        alpha=float(config.get("alpha", model.get("alpha", 0.05))),
        k=int(model["config"].get("knn", {}).get("k", 15)),
        regularization=float(model["config"].get("mahalanobis", {}).get("regularization", 0.01)),
    )
    cell_path = output_dir / "cell_normality_classes.csv"
    condition_path = output_dir / "condition_class_summary.csv"
    lineage_path = output_dir / "lineage_class_summary.csv"
    write_dataframe(results, cell_path)
    condition_cols = ["dataset_id", "condition", "perturbation_name", "time_point", "lineage"]
    write_dataframe(summarize_classes(results, condition_cols), condition_path)
    write_dataframe(summarize_classes(results, ["lineage", "perturbation_name"]), lineage_path)
    write_manifest(
        output_dir / "perturbation_classification_manifest.json",
        name="classify_perturbed_cells",
        inputs=[str(input_h5ad), str(model_path)],
        outputs=[str(cell_path), str(condition_path), str(lineage_path)],
        parameters=config,
        metrics={"n_perturbed_cells": int(results.shape[0])},
    )
    return cell_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Classify perturbed cells with DevGuard.")
    parser.add_argument("--config", default="config/devguard/perturbation_tests.json")
    args = parser.parse_args()
    classify_perturbed_cells(args.config)


if __name__ == "__main__":
    main()
