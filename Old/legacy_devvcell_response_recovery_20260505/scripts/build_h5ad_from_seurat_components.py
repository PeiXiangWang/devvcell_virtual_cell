"""Build an H5AD file from components exported from a Seurat RDS object."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy import io, sparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--components", default="data/external/response_transfer/GSE208369_components")
    parser.add_argument("--output", default="data/external/response_transfer/primary_perturbation.h5ad")
    return parser.parse_args()


def sanitize_obs(obs: pd.DataFrame) -> pd.DataFrame:
    obs = obs.copy()
    if "cell_barcode" in obs.columns:
        obs.index = obs["cell_barcode"].astype(str)
    else:
        obs.index = [f"cell_{i}" for i in range(len(obs))]
    obs.index.name = None
    for column in obs.columns:
        if pd.api.types.is_object_dtype(obs[column]):
            obs[column] = obs[column].fillna("").astype(str)
    return obs


def sanitize_var(var: pd.DataFrame) -> pd.DataFrame:
    var = var.copy()
    index_col = "gene_name" if "gene_name" in var.columns else "gene_id"
    var.index = var[index_col].astype(str)
    var.index.name = None
    for column in var.columns:
        if pd.api.types.is_object_dtype(var[column]):
            var[column] = var[column].fillna("").astype(str)
    return var


def main() -> None:
    args = parse_args()
    component_dir = Path(args.components)
    output = Path(args.output)
    matrix = io.mmread(component_dir / "matrix.mtx").tocsr()
    obs = sanitize_obs(pd.read_csv(component_dir / "obs.csv", low_memory=False))
    var = sanitize_var(pd.read_csv(component_dir / "var.csv", low_memory=False))

    # Matrix Market export from Seurat is genes x cells; AnnData expects cells x genes.
    X = sparse.csr_matrix(matrix.T)
    if X.shape[0] != len(obs) or X.shape[1] != len(var):
        raise ValueError(f"Shape mismatch: X={X.shape}, obs={len(obs)}, var={len(var)}")

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.uns["source"] = {
        "component_dir": str(component_dir),
        "conversion": "Seurat RDS -> Matrix Market/CSV -> H5AD",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output, compression="gzip")
    print(f"Wrote H5AD: {output}")
    print(f"Shape: {adata.n_obs} cells x {adata.n_vars} genes")


if __name__ == "__main__":
    main()
