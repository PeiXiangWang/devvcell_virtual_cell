"""Marker-gene module scoring helpers for DevGuard validation analyses."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import sparse


MARKER_MODULES: dict[str, list[str]] = {
    "haematoendothelial": ["Tal1", "Lmo2", "Kdr", "Etv2", "Pecam1", "Cdh5", "Runx1", "Gata2"],
    "erythroid": ["Gata1", "Klf1", "Hba-x", "Hbb-bh1", "Hbb-y", "Alas2", "Nfe2"],
    "endothelial": ["Kdr", "Pecam1", "Cdh5", "Sox7", "Sox17", "Tek", "Egfl7", "Flt1"],
    "primitive_streak_mesoderm": ["T", "Eomes", "Mixl1", "Mesp1", "Msgn1", "Tbx6"],
    "allantois_exe_mesoderm": ["Tbx4", "Hand1", "Cdx2", "Bmp4", "Pdgfra", "T", "Prdm1"],
    "cardiac_mesoderm": ["Mesp1", "Nkx2-5", "Tbx5", "Tnnt2", "Myh6", "Hand2", "Isl1"],
    "paraxial_somite": ["Msgn1", "Tbx6", "Meox1", "Mesp2", "Uncx", "Pax3"],
    "neural_ectoderm": ["Sox2", "Sox1", "Pax6", "Pou3f1", "Nes", "Otx2", "Six3"],
    "endoderm": ["Foxa2", "Sox17", "Gata4", "Gata6", "Apoa1", "Hhex"],
    "pluripotency_epiblast": ["Pou5f1", "Nanog", "Fgf5", "Utf1"],
    "wnt_response": ["Axin2", "Sp5", "Lef1", "Tcf7", "Wnt3", "Cdx2", "T"],
}


def _as_array(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix)


def gene_lookup(var: pd.DataFrame, var_names: Sequence[object]) -> dict[str, str]:
    """Build a case-insensitive gene symbol/id lookup to AnnData var names."""

    lookup: dict[str, str] = {}
    candidate_columns = ["SYMBOL", "gene_symbol", "gene_name", "gene_id", "ENSEMBL", "raw_feature_id"]
    for name in var_names:
        text = str(name)
        lookup.setdefault(text.lower(), text)
    for column in candidate_columns:
        if column not in var.columns:
            continue
        for var_name, value in zip(var_names, var[column].astype(str), strict=False):
            if value and value.lower() not in {"", "na", "nan", "none"}:
                lookup.setdefault(value.lower(), str(var_name))
    return lookup


def resolve_genes(adata, genes: Sequence[str]) -> dict[str, str]:
    lookup = gene_lookup(adata.var, adata.var_names)
    resolved: dict[str, str] = {}
    for gene in genes:
        match = lookup.get(str(gene).lower())
        if match is not None:
            resolved[str(gene)] = match
    return resolved


def log_normalized_expression(
    adata,
    genes: Sequence[str],
    *,
    normalize_total: float = 10000.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return log-normalized expression for available genes and a gene-presence table."""

    unique_genes = list(dict.fromkeys(str(gene) for gene in genes))
    resolved = resolve_genes(adata, unique_genes)
    presence = pd.DataFrame(
        [
            {
                "requested_gene": gene,
                "resolved_var_name": resolved.get(gene, ""),
                "available": gene in resolved,
            }
            for gene in unique_genes
        ]
    )
    if not resolved:
        return pd.DataFrame(index=adata.obs_names.astype(str)), presence

    matrix = adata.X
    totals = np.asarray(matrix.sum(axis=1)).ravel() if sparse.issparse(matrix) else np.asarray(matrix).sum(axis=1)
    scale = np.divide(
        normalize_total,
        totals,
        out=np.zeros_like(totals, dtype=float),
        where=totals > 0,
    )
    selected_names = [resolved[gene] for gene in unique_genes if gene in resolved]
    selected_labels = [gene for gene in unique_genes if gene in resolved]
    selected = adata[:, selected_names].X
    values = _as_array(selected).astype(float)
    values = np.log1p(values * scale[:, None])
    expression = pd.DataFrame(values, index=adata.obs_names.astype(str), columns=selected_labels)
    return expression, presence


def score_marker_modules(
    adata,
    modules: Mapping[str, Sequence[str]] | None = None,
    *,
    normalize_total: float = 10000.0,
    min_available_genes: int = 1,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score each cell by averaging log-normalized marker expression per module."""

    module_map = dict(modules or MARKER_MODULES)
    all_genes = [gene for genes in module_map.values() for gene in genes]
    expression, presence = log_normalized_expression(adata, all_genes, normalize_total=normalize_total)
    available_genes = set(presence.loc[presence["available"], "requested_gene"].astype(str))
    resolved_by_gene = presence.set_index("requested_gene")["resolved_var_name"].to_dict()
    scores = pd.DataFrame(index=adata.obs_names.astype(str))
    module_rows = []
    for module, genes in module_map.items():
        available = [gene for gene in genes if gene in available_genes]
        module_rows.append(
            {
                "module": module,
                "requested_genes": ";".join(genes),
                "available_genes": ";".join(available),
                "resolved_var_names": ";".join(str(resolved_by_gene.get(gene, "")) for gene in available),
                "n_requested_genes": len(genes),
                "n_available_genes": len(available),
            }
        )
        if len(available) >= min_available_genes:
            scores[module] = expression[available].mean(axis=1)
        else:
            scores[module] = np.nan
    return scores, pd.DataFrame(module_rows)


def assign_top_module(scores: pd.DataFrame) -> pd.Series:
    if scores.empty:
        return pd.Series(dtype=str)
    return scores.idxmax(axis=1).fillna("unassigned")
