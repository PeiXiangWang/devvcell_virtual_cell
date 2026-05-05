"""Build stage-lineage-module time-series tables from single-cell data."""

from __future__ import annotations

import numpy as np
import pandas as pd

from devguard.embedding import SVDEmbeddingModel
from devguard.preprocessing import standardize_obs
from devspectrum.modules import load_module_registry, module_long_table, score_modules


def _metadata_frame(adata, columns: list[str]) -> pd.DataFrame:
    metadata = adata.obs.copy()
    for column in columns:
        if column not in metadata.columns:
            metadata[column] = "NA"
    out = metadata[columns].copy()
    out.index = adata.obs_names.astype(str)
    return out


def _score_latent_dimensions(
    adata,
    metadata: pd.DataFrame,
    *,
    n_hvg: int = 2000,
    latent_dim: int = 5,
    seed: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    model = SVDEmbeddingModel(n_hvg=n_hvg, latent_dim=latent_dim, random_state=seed)
    embeddings = model.fit_transform(adata)
    scores = pd.DataFrame(
        embeddings[:, : min(latent_dim, embeddings.shape[1])],
        index=adata.obs_names.astype(str),
        columns=[f"latent_{i + 1}" for i in range(min(latent_dim, embeddings.shape[1]))],
    )
    coverage = pd.DataFrame(
        [
            {
                "module_name": column,
                "requested_genes": "",
                "available_genes": "",
                "resolved_var_names": "",
                "n_requested_genes": 0,
                "n_genes_in_module": 0,
                "module_gene_coverage": 1.0,
            }
            for column in scores.columns
        ]
    )
    long = module_long_table(scores, metadata, coverage)
    long["feature_type"] = "latent_dimension"
    long = long.rename(columns={"module_score": "feature_value"})
    return long, coverage


def build_stage_lineage_module_timeseries(
    adata,
    *,
    dataset_id: str,
    module_registry: str | None = None,
    gene_symbol_map: str | None = None,
    time_column: str = "time_point",
    time_numeric_column: str = "time_numeric",
    lineage_column: str = "lineage",
    sample_column: str = "sample_id",
    condition_column: str = "condition",
    min_cells_per_group: int = 20,
    include_module_scores: bool = True,
    include_latent_dimensions: bool = False,
    latent_dim: int = 5,
    seed: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate cells into sample-level developmental time-series features."""

    adata = standardize_obs(adata, dataset_id=dataset_id)
    metadata_columns = [
        "cell_id",
        "dataset_id",
        condition_column,
        sample_column,
        time_column,
        time_numeric_column,
        lineage_column,
    ]
    metadata = _metadata_frame(adata, metadata_columns)
    metadata = metadata.rename(
        columns={
            condition_column: "condition",
            sample_column: "sample_id",
            time_column: "time_point",
            time_numeric_column: "time_numeric",
            lineage_column: "lineage",
        }
    )
    metadata["time_numeric"] = pd.to_numeric(metadata["time_numeric"], errors="coerce")

    long_frames = []
    coverage_frames = []
    if include_module_scores:
        modules = load_module_registry(module_registry)
        scores, coverage = score_modules(adata, modules, gene_symbol_map=gene_symbol_map)
        long = module_long_table(scores, metadata, coverage).rename(columns={"module_score": "feature_value"})
        long["feature_type"] = "module_score"
        long_frames.append(long)
        coverage_frames.append(coverage.assign(feature_type="module_score"))
    if include_latent_dimensions:
        latent_long, latent_coverage = _score_latent_dimensions(
            adata,
            metadata,
            latent_dim=latent_dim,
            seed=seed,
        )
        long_frames.append(latent_long)
        coverage_frames.append(latent_coverage.assign(feature_type="latent_dimension"))

    if not long_frames:
        raise ValueError("No DevSpectrum features requested.")
    cell_features = pd.concat(long_frames, axis=0, ignore_index=True, sort=False)
    cell_features = cell_features[pd.to_numeric(cell_features["feature_value"], errors="coerce").notna()].copy()
    if cell_features.empty:
        raise ValueError("All requested DevSpectrum feature values are missing.")
    group_cols = ["dataset_id", "condition", "sample_id", "time_point", "time_numeric", "lineage", "module_name", "feature_type"]
    grouped = (
        cell_features.groupby(group_cols, dropna=False, observed=True)
        .agg(
            feature_value=("feature_value", "mean"),
            feature_std=("feature_value", "std"),
            n_cells=("cell_id", "size"),
            n_genes_in_module=("n_genes_in_module", "first"),
            module_gene_coverage=("module_gene_coverage", "first"),
        )
        .reset_index()
    )
    grouped = grouped[grouped["n_cells"] >= int(min_cells_per_group)].copy()
    grouped["feature_std"] = grouped["feature_std"].fillna(0.0)
    grouped["time_numeric"] = pd.to_numeric(grouped["time_numeric"], errors="coerce")
    coverage_table = pd.concat(coverage_frames, axis=0, ignore_index=True, sort=False).drop_duplicates(
        ["feature_type", "module_name"]
    )
    return grouped.sort_values(group_cols).reset_index(drop=True), coverage_table.reset_index(drop=True)


def aggregate_condition_timeseries(timeseries: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["dataset_id", "condition", "time_point", "time_numeric", "lineage", "module_name", "feature_type"]
    return (
        timeseries.groupby(group_cols, dropna=False, observed=True)
        .agg(
            feature_value=("feature_value", "mean"),
            feature_std=("feature_value", "std"),
            n_cells=("n_cells", "sum"),
            n_samples=("sample_id", "nunique"),
            n_genes_in_module=("n_genes_in_module", "first"),
            module_gene_coverage=("module_gene_coverage", "first"),
        )
        .reset_index()
        .assign(feature_std=lambda df: df["feature_std"].fillna(0.0))
        .sort_values(group_cols)
        .reset_index(drop=True)
    )


def make_quick_timeseries() -> pd.DataFrame:
    """Small deterministic fixture used by quick-mode tests and demos."""

    rng = np.random.default_rng(42)
    rows = []
    times = np.asarray([3.0, 3.5, 4.0, 4.5, 5.0])
    for lineage in ["mesodermal", "neural", "intermediate"]:
        for sample_idx in range(4):
            for time in times:
                for module_name in ["mesoderm", "neural", "stress"]:
                    trend = {"mesoderm": time - 3.0, "neural": 5.0 - time, "stress": 0.2 * np.sin(time * np.pi)}[module_name]
                    lineage_shift = 0.6 if module_name.startswith(lineage[:5]) else 0.0
                    rows.append(
                        {
                            "dataset_id": "devspectrum_quick",
                            "condition": "control",
                            "sample_id": f"s{sample_idx}",
                            "time_point": f"d{time:g}",
                            "time_numeric": time,
                            "lineage": lineage,
                            "module_name": module_name,
                            "feature_type": "module_score",
                            "feature_value": float(trend + lineage_shift + rng.normal(0, 0.03)),
                            "feature_std": 0.03,
                            "n_cells": 50,
                            "n_genes_in_module": 5,
                            "module_gene_coverage": 1.0,
                        }
                    )
    return pd.DataFrame(rows)
