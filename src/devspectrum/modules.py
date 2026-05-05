"""Module registry and module-score calculation for DevSpectrum."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import pandas as pd
import numpy as np

from devguard.markers import log_normalized_expression
from devguard.io import load_json


DEFAULT_MODULES: dict[str, list[str]] = {
    "mesoderm": ["T", "Eomes", "Mixl1", "Mesp1", "Msgn1", "Tbx6", "Mesp2", "Pdgfra"],
    "neural": ["Sox2", "Sox1", "Pax6", "Pou3f1", "Nes", "Otx2", "Six3", "Zic2"],
    "intermediate": ["Cdx2", "Fgf8", "Wnt3a", "Tbx6", "Pax3", "Meox1", "Nkx1-2"],
    "stress": ["Fos", "Jun", "Atf3", "Ddit3", "Hspa1a", "Hspa1b", "Hsp90aa1"],
    "cell_cycle": ["Mki67", "Top2a", "Pcna", "Cdk1", "Ccnb1", "Ccnb2", "Mcm5"],
    "apoptosis": ["Bax", "Bcl2l11", "Casp3", "Casp7", "Trp53", "Pmaip1", "Bbc3"],
    "haematoendothelial": ["Tal1", "Lmo2", "Kdr", "Etv2", "Pecam1", "Cdh5", "Runx1", "Gata2"],
    "erythroid": ["Gata1", "Klf1", "Hba-x", "Hbb-bh1", "Hbb-y", "Alas2", "Nfe2"],
    "endothelial": ["Kdr", "Pecam1", "Cdh5", "Sox7", "Sox17", "Tek", "Egfl7", "Flt1"],
    "cardiac": ["Mesp1", "Nkx2-5", "Tbx5", "Tnnt2", "Myh6", "Hand2", "Isl1"],
    "extraembryonic_mesoderm": ["Tbx4", "Hand1", "Cdx2", "Bmp4", "Pdgfra", "T", "Prdm1"],
    "paraxial_mesoderm": ["Msgn1", "Tbx6", "Meox1", "Mesp2", "Uncx", "Pax3"],
    "wnt_response": ["Axin2", "Sp5", "Lef1", "Tcf7", "Wnt3", "Cdx2", "T"],
}


def load_module_registry(path: str | Path | None = None) -> dict[str, list[str]]:
    if path is None:
        return dict(DEFAULT_MODULES)
    payload = load_json(path)
    modules = payload.get("modules", payload)
    return {str(name): [str(gene) for gene in genes] for name, genes in modules.items()}


def apply_gene_symbol_map(adata, map_path: str | Path | None = None) -> None:
    """Attach gene symbols to ``adata.var`` when source H5AD only has Ensembl IDs."""

    if map_path is None:
        return
    path = Path(map_path)
    if not path.exists():
        return
    frame = pd.read_csv(path, compression="infer")
    if frame.shape[1] == 1 and "," in frame.columns[0]:
        frame = pd.read_csv(path, compression="infer", sep=",")
    if {"ens_id", "symbol"}.issubset(frame.columns):
        mapping = frame.dropna(subset=["ens_id", "symbol"]).set_index("ens_id")["symbol"].astype(str).to_dict()
    elif {"gene_id", "gene_name"}.issubset(frame.columns):
        mapping = frame.dropna(subset=["gene_id", "gene_name"]).set_index("gene_id")["gene_name"].astype(str).to_dict()
    else:
        return
    var = adata.var.copy()
    ids = var["gene_id"].astype(str) if "gene_id" in var.columns else pd.Series(adata.var_names.astype(str), index=var.index)
    symbols = ids.map(mapping)
    if symbols.notna().any():
        var["SYMBOL"] = symbols.fillna(ids).astype(str)
        var["gene_symbol"] = var["SYMBOL"]
        adata.var = var


def score_modules(
    adata,
    modules: Mapping[str, Sequence[str]] | None = None,
    *,
    normalize_total: float = 10000.0,
    gene_symbol_map: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return per-cell module scores and per-module gene coverage."""

    apply_gene_symbol_map(adata, gene_symbol_map)
    module_map = dict(modules or DEFAULT_MODULES)
    all_genes = [gene for genes in module_map.values() for gene in genes]
    expression, presence = log_normalized_expression(adata, all_genes, normalize_total=normalize_total)
    available = set(presence.loc[presence["available"], "requested_gene"].astype(str))
    resolved = presence.set_index("requested_gene")["resolved_var_name"].to_dict()

    scores = pd.DataFrame(index=adata.obs_names.astype(str))
    rows = []
    for module_name, genes in module_map.items():
        present = [gene for gene in genes if gene in available]
        scores[module_name] = expression[present].mean(axis=1) if present else np.nan
        rows.append(
            {
                "module_name": module_name,
                "requested_genes": ";".join(genes),
                "available_genes": ";".join(present),
                "resolved_var_names": ";".join(str(resolved.get(gene, "")) for gene in present),
                "n_requested_genes": int(len(genes)),
                "n_genes_in_module": int(len(present)),
                "module_gene_coverage": float(len(present) / max(len(genes), 1)),
            }
        )
    return scores.astype(float), pd.DataFrame(rows)


def module_long_table(scores: pd.DataFrame, metadata: pd.DataFrame, coverage: pd.DataFrame) -> pd.DataFrame:
    """Convert wide module scores to a long table with coverage fields attached."""

    score_frame = scores.copy()
    score_frame.index = metadata.index
    long = (
        score_frame.reset_index(names="_cell_index")
        .melt(id_vars="_cell_index", var_name="module_name", value_name="module_score")
        .merge(metadata.reset_index(names="_cell_index"), on="_cell_index", how="left")
    )
    return long.merge(coverage, on="module_name", how="left")
