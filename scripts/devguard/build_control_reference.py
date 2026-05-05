from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import joblib

from devguard.embedding import SVDEmbeddingModel, apply_batch_centering, fit_batch_centering
from devguard.io import ensure_dir, load_json, read_h5ad, write_dataframe, write_manifest
from devguard.normality import build_normality_groups, centroids_frame, quality_frame, reference_cells_frame, thresholds_frame
from devguard.preprocessing import standardize_obs


def build_control_reference(config_path: str | Path) -> Path:
    config = load_json(config_path)
    input_h5ad = Path(config["input_h5ad"])
    output_dir = ensure_dir(config.get("output_dir", "results/devguard/normality_reference"))
    model_output = Path(config.get("model_output", output_dir / "devguard_normality_model.joblib"))
    adata = read_h5ad(input_h5ad)
    dataset_id = str(adata.obs["dataset_id"].iloc[0]) if "dataset_id" in adata.obs else input_h5ad.stem
    adata = standardize_obs(adata, dataset_id=dataset_id)
    control_mask = adata.obs["is_control"].astype(bool).to_numpy()
    control = adata[control_mask].copy()
    control_obs = control.obs.copy().reset_index(drop=True)

    embedding_cfg = config.get("embedding", {})
    embedding_model = SVDEmbeddingModel(
        n_hvg=int(embedding_cfg.get("n_hvg", 3000)),
        latent_dim=int(embedding_cfg.get("latent_dim", 50)),
        normalize_total=embedding_cfg.get("normalize_total", 10000),
        log1p=bool(embedding_cfg.get("log1p", True)),
        random_state=int(config.get("seed", 42)),
    )
    control_embeddings = embedding_model.fit_transform(control)
    batch_centering = None
    batch_center_column = embedding_cfg.get("batch_center_column")
    if batch_center_column:
        if batch_center_column not in control_obs.columns:
            raise ValueError(f"batch_center_column not found in obs: {batch_center_column}")
        batch_centering = fit_batch_centering(control_embeddings, control_obs, column=batch_center_column)
        control_embeddings = apply_batch_centering(control_embeddings, control_obs, batch_centering)

    grouping = config.get("reference_grouping", {})
    splits = config.get("splits", {})
    sample_split = config.get("sample_split", {})
    groups = build_normality_groups(
        control_embeddings,
        control_obs,
        time_column=grouping.get("time_column", "time_point"),
        time_numeric_column=grouping.get("time_numeric_column", "time_numeric"),
        lineage_column=grouping.get("lineage_column", "lineage"),
        min_cells_per_group=int(grouping.get("min_cells_per_group", 30)),
        train_fraction=float(splits.get("train_fraction", 0.6)),
        calibration_fraction=float(splits.get("calibration_fraction", 0.2)),
        split_strategy=sample_split.get("strategy", "cell"),
        split_unit_column=sample_split.get("unit_column", "sample_id"),
        allow_cell_split_fallback=bool(sample_split.get("allow_cell_fallback", True)),
        min_units_per_group=int(sample_split.get("min_units_per_group", 0)),
        seed=int(config.get("seed", 42)),
        score_methods=config.get("score_methods", ["knn_distance", "mahalanobis"]),
        k=int(config.get("knn", {}).get("k", 15)),
        regularization=float(config.get("mahalanobis", {}).get("regularization", 0.01)),
    )
    if not groups:
        raise ValueError("No normality groups were built. Check min_cells_per_group and control annotations.")

    alpha = float(config.get("alpha", 0.05))
    outputs = [
        output_dir / "reference_cells.parquet",
        output_dir / "stage_lineage_centroids.csv",
        output_dir / "conformal_thresholds.csv",
        output_dir / "normality_reference_quality.csv",
        model_output,
    ]
    write_dataframe(reference_cells_frame(groups, control_obs), outputs[0])
    write_dataframe(centroids_frame(groups), outputs[1])
    write_dataframe(thresholds_frame(groups, alpha), outputs[2])
    quality = quality_frame(groups, alpha)
    write_dataframe(quality, outputs[3])

    ensure_dir(model_output.parent)
    joblib.dump(
        {
            "embedding_model": embedding_model,
            "groups": groups,
            "config": config,
            "input_h5ad": str(input_h5ad),
            "score_methods": config.get("score_methods", ["knn_distance", "mahalanobis"]),
            "alpha": alpha,
            "batch_centering": batch_centering,
        },
        model_output,
    )
    write_manifest(
        output_dir / "normality_reference_manifest.json",
        name="build_control_reference",
        inputs=[str(input_h5ad)],
        outputs=[str(path) for path in outputs],
        parameters=config,
        metrics={
            "n_control_cells": int(control.n_obs),
            "n_reference_groups": int(len(groups)),
            "max_heldout_control_fpr": float(quality["heldout_control_fpr"].max()),
        },
    )
    return model_output


def main() -> None:
    parser = argparse.ArgumentParser(description="Build DevGuard control normality references.")
    parser.add_argument("--config", default="config/devguard/normality_model.json")
    args = parser.parse_args()
    build_control_reference(args.config)


if __name__ == "__main__":
    main()
