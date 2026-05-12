from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

LR_PAIRS = [
    ("Fgf8", "Fgfr1"),
    ("Wnt3", "Fzd1"),
    ("Kitl", "Kit"),
    ("Cxcl12", "Cxcr4"),
    ("Vegfa", "Kdr"),
    ("Dll1", "Notch1"),
    ("FGF8", "FGFR1"),
    ("WNT3", "FZD1"),
    ("KITLG", "KIT"),
    ("CXCL12", "CXCR4"),
    ("VEGFA", "KDR"),
    ("DLL1", "NOTCH1"),
]


def _gene_index(var: pd.DataFrame) -> dict[str, int]:
    names = [str(v) for v in var.index]
    for col in ("swarm_gene_symbol", "gene_short_name", "gene_name", "SYMBOL", "feature_name"):
        if col in var:
            names += [str(v) for v in var[col].tolist()]
    out = {}
    for i, name in enumerate([str(v) for v in var.index]):
        out[name.lower()] = i
    for col in ("swarm_gene_symbol", "gene_short_name", "gene_name", "SYMBOL", "feature_name"):
        if col in var:
            for i, name in enumerate(var[col].astype(str)):
                out.setdefault(name.lower(), i)
    return out


def cci_context(adata, labels: np.ndarray) -> tuple[np.ndarray, list[tuple[str, str]]]:
    index = _gene_index(adata.var)
    pairs = [(l, r) for l, r in LR_PAIRS if l.lower() in index and r.lower() in index]
    if not pairs:
        return np.zeros(adata.n_obs), []
    x = adata.X
    scores = np.zeros(adata.n_obs, dtype=float)
    for ligand, receptor in pairs:
        li, ri = index[ligand.lower()], index[receptor.lower()]
        lig = np.asarray(x[:, li].todense()).ravel() if sparse.issparse(x) else np.asarray(x[:, li]).ravel()
        rec = np.asarray(x[:, ri].todense()).ravel() if sparse.issparse(x) else np.asarray(x[:, ri]).ravel()
        scores += np.sqrt(np.maximum(lig, 0) * np.maximum(rec, 0))
    scores /= max(len(pairs), 1)
    lo, hi = np.quantile(scores, [0.05, 0.95]) if np.any(scores) else (0.0, 1.0)
    scores = np.clip((scores - lo) / max(hi - lo, 1e-8), 0.0, 1.0)
    return scores, pairs


def cci_delta(z: np.ndarray, labels: np.ndarray, cci_score: np.ndarray, strength: float = 0.08) -> np.ndarray:
    centers = {lab: z[labels == lab].mean(axis=0) for lab in sorted(set(labels))}
    global_center = z.mean(axis=0)
    out = np.zeros_like(z)
    for i, lab in enumerate(labels):
        out[i] = strength * cci_score[i] * (0.65 * (centers[lab] - z[i]) + 0.35 * (global_center - z[i]))
    return out

