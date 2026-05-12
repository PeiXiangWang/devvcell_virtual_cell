from __future__ import annotations

import numpy as np
import pandas as pd


def coupling_lineage_edges(npz: dict, obs: pd.DataFrame, cell_type_key: str) -> pd.DataFrame:
    src = npz["source_indices"]
    tgt = npz["target_indices"]
    plan = npz["plan"]
    src_types = obs.iloc[src][cell_type_key].astype(str).to_numpy()
    tgt_types = obs.iloc[tgt][cell_type_key].astype(str).to_numpy()
    rows = []
    for s_type in sorted(set(src_types)):
        s_mask = src_types == s_type
        for t_type in sorted(set(tgt_types)):
            mass = float(plan[np.ix_(s_mask, tgt_types == t_type)].sum())
            if mass > 0:
                rows.append({"source_lineage": s_type, "target_lineage": t_type, "mass": mass})
    return pd.DataFrame(rows)


def js_divergence(p: np.ndarray, q: np.ndarray) -> float:
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    p = p / max(p.sum(), 1e-12)
    q = q / max(q.sum(), 1e-12)
    m = 0.5 * (p + q)
    kl_pm = np.sum(np.where(p > 0, p * np.log(np.maximum(p, 1e-12) / np.maximum(m, 1e-12)), 0.0))
    kl_qm = np.sum(np.where(q > 0, q * np.log(np.maximum(q, 1e-12) / np.maximum(m, 1e-12)), 0.0))
    return float(0.5 * (kl_pm + kl_qm))

