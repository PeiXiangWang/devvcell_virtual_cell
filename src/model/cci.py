from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neighbors import NearestNeighbors

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


def sender_receiver_graph(adata, z: np.ndarray, labels: np.ndarray, k: int = 15) -> tuple[pd.DataFrame, np.ndarray, list[tuple[str, str]]]:
    """Build a LIANA/CellPhoneDB-style sender-receiver graph from LR expression."""
    scores, pairs = cci_context(adata, labels)
    if z.shape[0] <= 2 or not pairs:
        return pd.DataFrame(), scores, pairs
    k = min(k, z.shape[0] - 1)
    neigh = NearestNeighbors(n_neighbors=k + 1).fit(z).kneighbors(z, return_distance=False)[:, 1:]
    rows = []
    node_signal = np.zeros(z.shape[0], dtype=float)
    for i in range(z.shape[0]):
        recv = str(labels[i])
        senders = labels[neigh[i]].astype(str)
        sender_score = scores[neigh[i]]
        for sender in np.unique(senders):
            mask = senders == sender
            weight = float(sender_score[mask].mean() * scores[i])
            if weight > 0:
                rows.append({"sender": sender, "receiver": recv, "weight": weight, "n_neighbors": int(mask.sum())})
                node_signal[i] += weight
    if node_signal.max() > 0:
        node_signal = node_signal / node_signal.max()
    graph = pd.DataFrame(rows)
    if not graph.empty:
        graph = graph.groupby(["sender", "receiver"], as_index=False).agg(weight=("weight", "mean"), n_neighbors=("n_neighbors", "sum"))
    return graph, node_signal, pairs


def cci_branch_delta(z: np.ndarray, labels: np.ndarray, cci_signal: np.ndarray, fate_bias: np.ndarray | None = None, strength: float = 0.04) -> np.ndarray:
    centers = {lab: z[labels == lab].mean(axis=0) for lab in sorted(set(labels))}
    out = np.zeros_like(z)
    for i, lab in enumerate(labels):
        target = centers.get(str(lab), z.mean(axis=0))
        out[i] = strength * cci_signal[i] * (target - z[i])
    if fate_bias is not None and fate_bias.size:
        out *= (1.0 + np.nan_to_num(fate_bias[:, None], nan=0.0))
    return out
