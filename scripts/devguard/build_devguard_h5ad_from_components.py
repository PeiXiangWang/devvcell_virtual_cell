from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse

from devguard.io import write_dataframe, write_h5ad, write_manifest
from devguard.preprocessing import metadata_summary, parse_time_numeric, standardize_obs


def _read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False)


def _sanitize_obs(obs: pd.DataFrame) -> pd.DataFrame:
    obs = obs.copy()
    index_col = "cell_barcode" if "cell_barcode" in obs.columns else ("cell" if "cell" in obs.columns else None)
    if index_col:
        obs.index = obs[index_col].astype(str)
    else:
        obs.index = [f"cell_{i}" for i in range(len(obs))]
    obs.index.name = None
    for column in obs.columns:
        if pd.api.types.is_object_dtype(obs[column]):
            obs[column] = obs[column].fillna("").astype(str)
    return obs


def _sanitize_var(var: pd.DataFrame) -> pd.DataFrame:
    var = var.copy()
    index_col = "gene_name" if "gene_name" in var.columns else ("SYMBOL" if "SYMBOL" in var.columns else "gene_id")
    var.index = var[index_col].fillna(var.get("gene_id", "")).astype(str)
    if not var.index.is_unique:
        counts: dict[str, int] = {}
        unique_names = []
        for name in var.index:
            count = counts.get(name, 0)
            unique_names.append(name if count == 0 else f"{name}-{count}")
            counts[name] = count + 1
        var.index = unique_names
    var.index.name = None
    return var


def _column(obs: pd.DataFrame, *names: str, default: object = "NA") -> pd.Series:
    for name in names:
        if name in obs.columns:
            return obs[name]
    return pd.Series(default, index=obs.index)


def infer_public_obs(obs: pd.DataFrame, *, dataset_id: str, dataset_kind: str, perturbation_name: str) -> pd.DataFrame:
    obs = obs.copy()
    if dataset_kind == "gse212050":
        obs["cell_id"] = _column(obs, "cell", default=obs.index.astype(str)).astype(str)
        obs["dataset_id"] = dataset_id
        obs["species"] = "Mus musculus"
        obs["system"] = "gastruloid"
        obs["time_point"] = _column(obs, "timepoint.demultiplexed", "timepoint", "day", default="unknown").astype(str)
        time_numeric = _column(obs, "day", default=obs["time_point"]).map(parse_time_numeric)
        obs["time_numeric"] = time_numeric.where(~time_numeric.isna(), obs["time_point"].map(parse_time_numeric))
        obs["condition"] = "control"
        obs["perturbation_name"] = "control"
        obs["perturbation_type"] = "none"
        sample = _column(obs, "sample", default=obs.index.astype(str)).astype(str)
        multi_class = _column(obs, "MULTI_class", default="").astype(str)
        obs["sample_id"] = np.where(multi_class.str.startswith("Bar"), multi_class, sample)
        obs["batch"] = _column(obs, "batch", "lane", default="NA").astype(str)
        obs["cell_type"] = _column(obs, "celltype.mapped.extended", "celltype.mapped.original", "cluster", default="NA").astype(str)
        lineage = _column(obs, "gastr_type", default=obs["cell_type"]).astype("string")
        lineage = lineage.mask(lineage.isna() | lineage.eq("") | lineage.eq("NA"), obs["cell_type"].astype("string"))
        obs["lineage"] = lineage.astype(str)
    elif dataset_kind in {"embryo_atlas", "mgd_atlas"}:
        obs["cell_id"] = _column(obs, "cell", "cell_barcode", default=obs.index.astype(str)).astype(str)
        obs["dataset_id"] = dataset_id
        obs["species"] = "Mus musculus"
        obs["system"] = "embryo"
        obs["time_point"] = _column(obs, "stage", default="unknown").astype(str)
        obs["time_numeric"] = obs["time_point"].map(parse_time_numeric)
        obs["condition"] = "control"
        obs["perturbation_name"] = "control"
        obs["perturbation_type"] = "none"
        obs["sample_id"] = _column(obs, "sample", "pool", default=obs.index.astype(str)).astype(str)
        obs["batch"] = _column(obs, "sequencing.batch", default="NA").astype(str)
        obs["cell_type"] = _column(obs, "celltype", "celltype.mapped", "cluster", default="NA").astype(str)
        obs["lineage"] = obs["cell_type"]
    elif dataset_kind in {"tal1_chimera", "t_chimera", "wt_chimera", "mgd_chimera"}:
        obs["cell_id"] = _column(obs, "cell", "cell_barcode", default=obs.index.astype(str)).astype(str)
        obs["dataset_id"] = dataset_id
        obs["species"] = "Mus musculus"
        obs["system"] = "embryo_chimera"
        obs["time_point"] = _column(obs, "stage", "theiler", default="unknown").astype(str)
        obs["time_numeric"] = obs["time_point"].map(parse_time_numeric)
        tomato = _column(obs, "tomato", default=False)
        tomato_bool = pd.Series(tomato, index=obs.index).astype(str).str.lower().isin(["true", "1", "yes", "tomato"])
        is_real_ko = dataset_kind in {"tal1_chimera", "t_chimera", "mgd_chimera"} and perturbation_name != "control"
        obs["condition"] = np.where(tomato_bool & is_real_ko, "perturbation", "control")
        obs["perturbation_name"] = np.where(obs["condition"].eq("perturbation"), perturbation_name, "control")
        obs["perturbation_type"] = np.where(obs["condition"].eq("perturbation"), "genetic", "none")
        obs["sample_id"] = _column(obs, "sample", "pool", default=obs.index.astype(str)).astype(str)
        obs["batch"] = _column(obs, "sequencing.batch", default="NA").astype(str)
        obs["cell_type"] = _column(obs, "celltype", "celltype.mapped", "cluster", default="NA").astype(str)
        obs["lineage"] = obs["cell_type"]
    elif dataset_kind == "gse123187":
        obs["cell_id"] = _column(obs, "cell", "cell_barcode", default=obs.index.astype(str)).astype(str)
        obs["dataset_id"] = dataset_id
        obs["species"] = "Mus musculus"
        obs["system"] = "gastruloid"
        obs["time_point"] = _column(obs, "time_point", "sample", default="unknown").astype(str)
        obs["time_numeric"] = obs["time_point"].map(parse_time_numeric)
        obs["condition"] = "control"
        obs["perturbation_name"] = "control"
        obs["perturbation_type"] = "none"
        obs["sample_id"] = _column(obs, "sample", default=obs.index.astype(str)).astype(str)
        obs["batch"] = _column(obs, "run", default="NA").astype(str)
        obs["cell_type"] = _column(obs, "celltype", default="NA").astype(str)
        obs["lineage"] = obs["cell_type"]
    else:
        raise ValueError(f"Unsupported dataset_kind: {dataset_kind}")

    for column in ["dose", "duration"]:
        if column not in obs.columns:
            obs[column] = "NA"
    obs["is_control"] = obs["condition"].astype(str).str.lower().eq("control")
    obs["is_perturbed"] = ~obs["is_control"]
    return obs


def build_h5ad_from_components(
    components: str | Path,
    output: str | Path,
    *,
    dataset_id: str,
    dataset_kind: str,
    perturbation_name: str = "control",
) -> Path:
    component_dir = Path(components)
    matrix = io.mmread(component_dir / "matrix.mtx").tocsr()
    obs = _sanitize_obs(_read_table(component_dir / "obs.csv"))
    var = _sanitize_var(_read_table(component_dir / "var.csv"))
    X = sparse.csr_matrix(matrix.T)
    if X.shape[0] != len(obs) or X.shape[1] != len(var):
        raise ValueError(f"Shape mismatch: X={X.shape}, obs={len(obs)}, var={len(var)}")
    obs = infer_public_obs(obs, dataset_id=dataset_id, dataset_kind=dataset_kind, perturbation_name=perturbation_name)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names_make_unique()
    adata = standardize_obs(adata, dataset_id=dataset_id)
    output_path = write_h5ad(adata, output)
    metadata_path = Path("results/devguard/dataset_metadata") / f"{dataset_id}_metadata.csv"
    write_dataframe(metadata_summary(adata.obs), metadata_path)
    write_manifest(
        Path("results/devguard/dataset_metadata") / f"{dataset_id}_h5ad_manifest.json",
        name="build_devguard_h5ad_from_components",
        inputs=[str(component_dir)],
        outputs=[str(output_path), str(metadata_path)],
        parameters={"dataset_id": dataset_id, "dataset_kind": dataset_kind, "perturbation_name": perturbation_name},
        metrics={"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars)},
    )
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert exported public dataset components into DevGuard H5AD.")
    parser.add_argument("--components", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--dataset-kind", required=True)
    parser.add_argument("--perturbation-name", default="control")
    args = parser.parse_args()
    build_h5ad_from_components(
        args.components,
        args.output,
        dataset_id=args.dataset_id,
        dataset_kind=args.dataset_kind,
        perturbation_name=args.perturbation_name,
    )


if __name__ == "__main__":
    main()
