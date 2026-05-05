"""Embedding model used by the DevGuard MVP."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import sparse
from sklearn.decomposition import TruncatedSVD


def _as_matrix(adata):
    matrix = adata.X
    if sparse.issparse(matrix):
        return matrix.tocsr().astype(float)
    return np.asarray(matrix, dtype=float)


def _normalize_total_log1p(matrix, normalize_total: float | None, log1p: bool):
    if normalize_total:
        if sparse.issparse(matrix):
            row_sums = np.asarray(matrix.sum(axis=1)).ravel()
            scale = np.divide(
                normalize_total,
                row_sums,
                out=np.zeros_like(row_sums, dtype=float),
                where=row_sums > 0,
            )
            matrix = matrix.multiply(scale[:, None]).tocsr()
        else:
            row_sums = matrix.sum(axis=1)
            scale = np.divide(
                normalize_total,
                row_sums,
                out=np.zeros_like(row_sums, dtype=float),
                where=row_sums > 0,
            )
            matrix = matrix * scale[:, None]
    if log1p:
        if sparse.issparse(matrix):
            matrix = matrix.copy()
            matrix.data = np.log1p(matrix.data)
        else:
            matrix = np.log1p(matrix)
    return matrix


def _feature_variance(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        mean = np.asarray(matrix.mean(axis=0)).ravel()
        sq_mean = np.asarray(matrix.multiply(matrix).mean(axis=0)).ravel()
        return np.maximum(sq_mean - mean**2, 0.0)
    return np.var(matrix, axis=0)


@dataclass
class SVDEmbeddingModel:
    """Small, auditable count-to-latent embedding model."""

    n_hvg: int = 3000
    latent_dim: int = 50
    normalize_total: float | None = 10000
    log1p: bool = True
    random_state: int = 42
    selected_genes_: list[str] = field(default_factory=list)
    svd_: TruncatedSVD | None = None

    def fit_transform(self, adata) -> np.ndarray:
        matrix = _normalize_total_log1p(_as_matrix(adata), self.normalize_total, self.log1p)
        var = _feature_variance(matrix)
        n_features = matrix.shape[1]
        n_select = min(max(1, self.n_hvg), n_features)
        selected_idx = np.argsort(var)[::-1][:n_select]
        selected_idx = np.sort(selected_idx)
        self.selected_genes_ = [str(adata.var_names[i]) for i in selected_idx]
        selected = matrix[:, selected_idx]
        n_components = min(self.latent_dim, max(1, selected.shape[0] - 1), max(1, selected.shape[1] - 1))
        self.svd_ = TruncatedSVD(n_components=n_components, random_state=self.random_state)
        return np.asarray(self.svd_.fit_transform(selected), dtype=float)

    def transform(self, adata) -> np.ndarray:
        if self.svd_ is None or not self.selected_genes_:
            raise RuntimeError("SVDEmbeddingModel must be fitted before transform().")
        var_index = {str(name): idx for idx, name in enumerate(adata.var_names)}
        missing = [gene for gene in self.selected_genes_ if gene not in var_index]
        if missing:
            preview = ", ".join(missing[:5])
            raise ValueError(f"Input data is missing {len(missing)} selected genes: {preview}")
        selected_idx = [var_index[gene] for gene in self.selected_genes_]
        matrix = _normalize_total_log1p(_as_matrix(adata), self.normalize_total, self.log1p)
        selected = matrix[:, selected_idx]
        return np.asarray(self.svd_.transform(selected), dtype=float)


def fit_batch_centering(embeddings: np.ndarray, obs, *, column: str) -> dict:
    labels = obs[column].astype("string").fillna("NA").astype(str)
    global_center = np.asarray(embeddings.mean(axis=0), dtype=float)
    centers = {}
    for label in sorted(labels.unique()):
        mask = labels.eq(label).to_numpy()
        centers[label] = np.asarray(embeddings[mask].mean(axis=0), dtype=float)
    return {"column": column, "global_center": global_center, "centers": centers}


def apply_batch_centering(
    embeddings: np.ndarray,
    obs,
    state: dict | None,
    *,
    fallback_columns: list[str] | tuple[str, ...] = (),
) -> np.ndarray:
    if not state:
        return embeddings
    columns = [state["column"], *fallback_columns]
    labels = None
    for column in columns:
        if column in obs.columns:
            labels = obs[column].astype("string").fillna("NA").astype(str)
            break
    if labels is None:
        return embeddings - state["global_center"]
    centered = np.asarray(embeddings, dtype=float).copy()
    centers = state["centers"]
    global_center = state["global_center"]
    for label in sorted(labels.unique()):
        mask = labels.eq(label).to_numpy()
        centered[mask] = centered[mask] - centers.get(label, global_center)
    return centered
