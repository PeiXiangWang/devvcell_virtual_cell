from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import joblib
import numpy as np

from devguard.classification import classify_cells_against_reference, summarize_classes
from devguard.embedding import apply_batch_centering
from devguard.io import ensure_dir, read_h5ad, write_dataframe, write_manifest
from devguard.preprocessing import standardize_obs


def classify_reference_test_cells(
    model_path: str | Path,
    output_dir: str | Path,
    *,
    score_method: str = "knn_distance",
) -> Path:
    model = joblib.load(model_path)
    input_h5ad = Path(model["input_h5ad"])
    adata = read_h5ad(input_h5ad)
    dataset_id = str(adata.obs["dataset_id"].iloc[0]) if "dataset_id" in adata.obs else input_h5ad.stem
    adata = standardize_obs(adata, dataset_id=dataset_id)
    control = adata[adata.obs["is_control"].astype(bool).to_numpy()].copy()
    test_indices = np.unique(np.concatenate([group.test_indices for group in model["groups"].values()]))
    test = control[test_indices].copy()
    obs = test.obs.copy().reset_index(drop=True)
    embeddings = model["embedding_model"].transform(test)
    embeddings = apply_batch_centering(
        embeddings,
        obs,
        model.get("batch_centering"),
        fallback_columns=model.get("config", {}).get("embedding", {}).get("batch_center_fallback_columns", ["dataset_id"]),
    )
    config = model.get("config", {})
    results = classify_cells_against_reference(
        embeddings,
        obs,
        model["groups"],
        score_method=score_method,
        alpha=float(model.get("alpha", 0.05)),
        ambiguous_margin=float(config.get("ambiguous_margin", 0.02)),
        k=int(config.get("knn", {}).get("k", 15)),
        regularization=float(config.get("mahalanobis", {}).get("regularization", 0.01)),
    )
    output = ensure_dir(output_dir)
    cell_path = output / "heldout_control_normality_classes.csv"
    summary_path = output / "heldout_control_class_summary.csv"
    write_dataframe(results, cell_path)
    write_dataframe(summarize_classes(results, ["dataset_id", "time_point", "lineage"]), summary_path)
    write_manifest(
        output / "heldout_control_classification_manifest.json",
        name="classify_reference_test_cells",
        inputs=[str(model_path), str(input_h5ad)],
        outputs=[str(cell_path), str(summary_path)],
        parameters={"score_method": score_method},
        metrics={"n_heldout_control_cells": int(results.shape[0])},
    )
    return cell_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Classify heldout control reference-test cells against a DevGuard model.")
    parser.add_argument("--model-path", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--score-method", default="knn_distance")
    args = parser.parse_args()
    classify_reference_test_cells(args.model_path, args.output_dir, score_method=args.score_method)


if __name__ == "__main__":
    main()
