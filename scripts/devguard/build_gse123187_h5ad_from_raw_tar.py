from __future__ import annotations

import argparse
import gzip
import re
import tarfile
from pathlib import Path

import _bootstrap  # noqa: F401
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from devguard.io import write_dataframe, write_h5ad, write_manifest
from devguard.preprocessing import metadata_summary, parse_time_numeric, standardize_obs


def _sample_id_from_name(name: str) -> str:
    return Path(name).name.split("_", 1)[0]


def _time_point_from_name(name: str) -> str:
    match = re.search(r"-(\d+(?:\.\d+)?)dAA-", name)
    if match:
        return f"{match.group(1)}dAA"
    return "unknown"


def _run_from_name(name: str) -> str:
    match = re.search(r"-(Run[^_]+)_", name)
    return match.group(1) if match else "NA"


def _mode_from_name(name: str) -> str:
    lower = name.lower()
    if "tomo" in lower:
        return "tomo_seq"
    if "single" in lower:
        return "single_gastruloid"
    return "unknown"


def _read_count_member(tar: tarfile.TarFile, name: str) -> tuple[sparse.csr_matrix, pd.DataFrame, pd.DataFrame]:
    member = tar.extractfile(name)
    if member is None:
        raise ValueError(f"Cannot read tar member: {name}")
    with gzip.GzipFile(fileobj=member) as handle:
        frame = pd.read_csv(handle, sep="\t", index_col=0)
    matrix = sparse.csr_matrix(frame.to_numpy(dtype=np.float32).T)
    gene_parts = frame.index.astype(str).str.split("_", n=2, expand=True)
    var = pd.DataFrame(
        {
            "gene_id": gene_parts.get_level_values(0),
            "gene_symbol": gene_parts.get_level_values(1) if gene_parts.nlevels > 1 else frame.index.astype(str),
            "raw_feature_id": frame.index.astype(str),
        },
        index=frame.index.astype(str),
    )
    sample_id = _sample_id_from_name(name)
    time_point = _time_point_from_name(name)
    obs = pd.DataFrame(
        {
            "cell_id": [f"{sample_id}_{cell}" for cell in frame.columns.astype(str)],
            "dataset_id": "GSE123187",
            "species": "Mus musculus",
            "system": "gastruloid",
            "time_point": time_point,
            "time_numeric": parse_time_numeric(time_point),
            "condition": "control",
            "perturbation_name": "control",
            "perturbation_type": "none",
            "dose": "NA",
            "duration": "NA",
            "sample_id": sample_id,
            "batch": _run_from_name(name),
            "cell_type": "NA",
            "lineage": "NA",
            "is_control": True,
            "is_perturbed": False,
            "gse123187_mode": _mode_from_name(name),
            "source_member": name,
            "tomo_position": list(frame.columns.astype(str)),
        }
    )
    obs.index = obs["cell_id"].astype(str)
    return matrix, obs, var


def build_gse123187_h5ad_from_raw_tar(
    raw_tar: str | Path,
    output: str | Path,
    *,
    max_files: int | None = None,
    count_suffix: str = ".coutb.tsv.gz",
) -> Path:
    raw_tar = Path(raw_tar)
    matrices = []
    obs_frames = []
    var_reference: pd.DataFrame | None = None
    selected_names = []
    with tarfile.open(raw_tar) as tar:
        names = [name for name in tar.getnames() if name.endswith(count_suffix)]
        names = sorted(names)
        if max_files is not None:
            names = names[:max_files]
        if not names:
            raise ValueError(f"No members ending with {count_suffix} found in {raw_tar}")
        for name in names:
            matrix, obs, var = _read_count_member(tar, name)
            if var_reference is None:
                var_reference = var
            elif not var.index.equals(var_reference.index):
                matrix = sparse.csr_matrix(
                    pd.DataFrame.sparse.from_spmatrix(matrix.T, index=var.index).reindex(var_reference.index).fillna(0).sparse.to_coo().T
                )
            matrices.append(matrix)
            obs_frames.append(obs)
            selected_names.append(name)
    X = sparse.vstack(matrices, format="csr")
    obs = pd.concat(obs_frames, axis=0)
    assert var_reference is not None
    adata = ad.AnnData(X=X, obs=obs, var=var_reference)
    adata.var_names_make_unique()
    adata = standardize_obs(adata, dataset_id="GSE123187")
    output_path = write_h5ad(adata, output)
    metadata_path = Path("results/devguard/dataset_metadata/GSE123187_metadata.csv")
    write_dataframe(metadata_summary(adata.obs), metadata_path)
    write_manifest(
        Path("results/devguard/dataset_metadata/GSE123187_h5ad_manifest.json"),
        name="build_gse123187_h5ad_from_raw_tar",
        inputs=[str(raw_tar)],
        outputs=[str(output_path), str(metadata_path)],
        parameters={"max_files": max_files, "count_suffix": count_suffix},
        metrics={"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars), "n_files": len(selected_names)},
    )
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a DevGuard H5AD from GSE123187 RAW tar count matrices.")
    parser.add_argument("--raw-tar", default="data/external/GSE123187/GSE123187_RAW.tar")
    parser.add_argument("--output", default="data/processed/devguard/GSE123187_preview.h5ad")
    parser.add_argument("--max-files", type=int, default=4)
    parser.add_argument("--count-suffix", default=".coutb.tsv.gz")
    args = parser.parse_args()
    build_gse123187_h5ad_from_raw_tar(
        args.raw_tar,
        args.output,
        max_files=args.max_files,
        count_suffix=args.count_suffix,
    )


if __name__ == "__main__":
    main()
