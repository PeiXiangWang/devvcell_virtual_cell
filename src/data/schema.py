from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd


TIME_PATTERNS = ("time", "day", "stage", "pseudotime", "somite", "theiler")
CELLTYPE_PATTERNS = ("cell_type", "celltype", "cluster", "fate", "lineage", "annotation")
BARCODE_PATTERNS = ("barcode", "clone", "lineage_tracing", "lineage_barcode")
CONDITION_PATTERNS = ("batch", "donor", "condition", "treatment", "perturb", "stim", "disease", "sample")
SPATIAL_PATTERNS = ("spatial", "x_centroid", "y_centroid", "x_coord", "y_coord", "tomo_axis", "position")
MULTIMODAL_PATTERNS = ("atac", "protein", "adt", "cite", "multiome", "peak", "modality")
CELL_CYCLE_MARKERS = {
    "human": ["MKI67", "TOP2A", "PCNA", "MCM5", "TYMS", "UBE2C"],
    "mouse": ["Mki67", "Top2a", "Pcna", "Mcm5", "Tyms", "Ube2c"],
}
LR_GENES = [
    "Fgf8",
    "Fgfr1",
    "Wnt3",
    "Fzd1",
    "Kitl",
    "Kit",
    "Cxcl12",
    "Cxcr4",
    "Vegfa",
    "Kdr",
    "Dll1",
    "Notch1",
    "FGF8",
    "FGFR1",
    "WNT3",
    "FZD1",
    "KITLG",
    "KIT",
    "CXCL12",
    "CXCR4",
    "VEGFA",
    "KDR",
    "DLL1",
    "NOTCH1",
]


def _matching_columns(columns: Iterable[str], patterns: Iterable[str]) -> list[str]:
    hits = []
    for col in columns:
        low = str(col).lower()
        if any(pattern in low for pattern in patterns):
            hits.append(str(col))
    return hits


def gene_symbols(adata: ad.AnnData) -> list[str]:
    for col in ("gene_short_name", "gene_name", "SYMBOL", "symbol", "feature_name", "gene_symbol"):
        if col in adata.var:
            vals = adata.var[col].astype(str).tolist()
            if vals:
                return vals
    return [str(v) for v in adata.var_names]


def inspect_h5ad(path: Path) -> dict[str, Any]:
    out: dict[str, Any] = {"path": str(path), "format": "h5ad", "readable": False}
    try:
        data = ad.read_h5ad(path, backed="r")
        obs_cols = [str(c) for c in data.obs.columns]
        var_cols = [str(c) for c in data.var.columns]
        genes = gene_symbols(data)
        lower_genes = {g.lower() for g in genes[: min(len(genes), 50000)]}
        cycle = [g for geneset in CELL_CYCLE_MARKERS.values() for g in geneset if g.lower() in lower_genes]
        lr = [g for g in LR_GENES if g.lower() in lower_genes]
        out.update(
            {
                "readable": True,
                "shape": [int(data.n_obs), int(data.n_vars)],
                "obs_columns": obs_cols,
                "var_columns": var_cols,
                "obsm_keys": [str(k) for k in data.obsm.keys()],
                "layers": [str(k) for k in data.layers.keys()],
                "uns_keys": [str(k) for k in data.uns.keys()],
                "expression_matrix": "X",
                "time_fields": _matching_columns(obs_cols, TIME_PATTERNS),
                "cell_type_fields": _matching_columns(obs_cols, CELLTYPE_PATTERNS),
                "lineage_barcode_fields": _matching_columns(obs_cols, BARCODE_PATTERNS),
                "condition_fields": _matching_columns(obs_cols, CONDITION_PATTERNS),
                "spatial_fields": _matching_columns(obs_cols, SPATIAL_PATTERNS)
                + [k for k in data.obsm.keys() if "spatial" in str(k).lower()],
                "multimodal_fields": _matching_columns(obs_cols + var_cols + [str(k) for k in data.obsm.keys()], MULTIMODAL_PATTERNS),
                "cell_cycle_marker_hits": cycle,
                "ligand_receptor_gene_hits": lr,
                "has_lineage_tracing": bool(_matching_columns(obs_cols, ("clone", "barcode", "lineage_tracing"))),
                "has_perturbation_labels": bool(_matching_columns(obs_cols, ("perturb", "treatment", "stim", "condition"))),
            }
        )
        data.file.close()
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"
    return out


def inspect_mtx_bundle(path: Path) -> dict[str, Any]:
    out: dict[str, Any] = {"path": str(path), "format": "mtx", "readable": False}
    try:
        shape = None
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if line.startswith("%"):
                    continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    shape = [int(parts[0]), int(parts[1])]
                    break
        if shape is None:
            raise ValueError("Could not parse MatrixMarket dimensions")
        obs_path = path.with_name("obs.csv")
        var_path = path.with_name("var.csv")
        obs_cols: list[str] = []
        var_cols: list[str] = []
        if obs_path.exists():
            obs_cols = [str(c) for c in pd.read_csv(obs_path, nrows=5).columns]
        if var_path.exists():
            var_cols = [str(c) for c in pd.read_csv(var_path, nrows=5).columns]
        out.update(
            {
                "readable": True,
                "shape": shape,
                "obs_columns": obs_cols,
                "var_columns": var_cols,
                "expression_matrix": "matrix.mtx",
                "time_fields": _matching_columns(obs_cols, TIME_PATTERNS),
                "cell_type_fields": _matching_columns(obs_cols, CELLTYPE_PATTERNS),
                "lineage_barcode_fields": _matching_columns(obs_cols, BARCODE_PATTERNS),
                "condition_fields": _matching_columns(obs_cols, CONDITION_PATTERNS),
                "spatial_fields": _matching_columns(obs_cols, SPATIAL_PATTERNS),
                "multimodal_fields": _matching_columns(obs_cols + var_cols, MULTIMODAL_PATTERNS),
            }
        )
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"
    return out


def summarize_records(records: list[dict[str, Any]]) -> dict[str, Any]:
    readable = [r for r in records if r.get("readable")]
    return {
        "n_files_scanned": len(records),
        "n_readable_expression_objects": len(readable),
        "formats": dict(pd.Series([r.get("format", "unknown") for r in records]).value_counts()) if records else {},
        "time_capable_files": [r["path"] for r in readable if r.get("time_fields")],
        "perturbation_capable_files": [r["path"] for r in readable if r.get("has_perturbation_labels") or r.get("condition_fields")],
        "lineage_or_celltype_files": [r["path"] for r in readable if r.get("cell_type_fields")],
        "spatial_candidate_files": [r["path"] for r in readable if r.get("spatial_fields")],
        "ligand_receptor_candidate_files": [r["path"] for r in readable if r.get("ligand_receptor_gene_hits")],
    }


def markdown_audit(records: list[dict[str, Any]], summary: dict[str, Any]) -> str:
    lines = [
        "# SwarmLineage-OT Data Audit",
        "",
        "This audit was generated from local files only. It records fields that are available for modelling; it does not certify biological validity.",
        "",
        "## Summary",
        "",
        f"- Files scanned: {summary['n_files_scanned']}",
        f"- Readable expression objects: {summary['n_readable_expression_objects']}",
        f"- Time/stage-capable files: {len(summary['time_capable_files'])}",
        f"- Perturbation/condition-capable files: {len(summary['perturbation_capable_files'])}",
        f"- Spatial candidate files: {len(summary['spatial_candidate_files'])}",
        f"- Ligand-receptor candidate files: {len(summary['ligand_receptor_candidate_files'])}",
        "",
        "## File-Level Schema",
        "",
    ]
    for r in records:
        lines += [
            f"### `{r['path']}`",
            "",
            f"- format: {r.get('format')}",
            f"- readable: {r.get('readable')}",
        ]
        if r.get("shape"):
            lines.append(f"- shape: {r['shape'][0]} cells/rows × {r['shape'][1]} genes/features")
        if r.get("error"):
            lines.append(f"- error: {r['error']}")
        for key, label in [
            ("time_fields", "time/day/stage/pseudotime"),
            ("cell_type_fields", "cell type/cluster/fate/lineage"),
            ("lineage_barcode_fields", "lineage/barcode"),
            ("condition_fields", "batch/donor/condition/treatment/perturbation"),
            ("spatial_fields", "spatial coordinates"),
            ("multimodal_fields", "multimodal fields"),
            ("cell_cycle_marker_hits", "cell-cycle/proliferation markers"),
            ("ligand_receptor_gene_hits", "ligand-receptor genes"),
        ]:
            vals = r.get(key) or []
            lines.append(f"- {label}: {', '.join(map(str, vals[:20])) if vals else 'not detected'}")
        lines.append(f"- lineage-tracing evidence likely present: {bool(r.get('has_lineage_tracing'))}")
        lines.append(f"- perturb-seq/treatment labels likely present: {bool(r.get('has_perturbation_labels'))}")
        lines.append("")
    lines += [
        "## Main Modelling Input Decision",
        "",
        "`data/processed/cell_level_subset_v1.h5ad` is used as the default quick real-data input because it has eight ordered developmental stages, lineage/cell-type labels, and a manageable 19,156 × 3,000 matrix. The 12.6GB `data/scLine_pro.h5ad` is retained as a larger extension target.",
        "",
        "If no real time field exists, stage or pseudotime is treated as an exploratory fallback. Main claims in generated manuscripts are gated accordingly.",
        "",
    ]
    return "\n".join(lines)
