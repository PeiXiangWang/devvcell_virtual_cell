from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def make_synthetic_lineage_adata(
    n_times: int = 5,
    cells_per_time: int = 90,
    n_genes: int = 120,
    seed: int = 17,
) -> ad.AnnData:
    """Create a small nontrivial time-series AnnData for smoke tests.

    The fixture has two branching lineages, stage-specific growth, cell-cycle
    markers and ligand/receptor genes. It is synthetic and only used for CI and
    leakage checks.
    """
    rng = np.random.default_rng(seed)
    genes = [f"Gene{i}" for i in range(n_genes - 10)] + [
        "Mki67",
        "Top2a",
        "Pcna",
        "Fgf8",
        "Fgfr1",
        "Cxcl12",
        "Cxcr4",
        "Wnt3",
        "Fzd1",
        "Kdr",
    ]
    rows = []
    mats = []
    for t in range(n_times):
        time_value = 12.0 + t
        p_neural = 0.25 + 0.1 * t
        for i in range(cells_per_time):
            lineage = "neural" if rng.random() < p_neural else "mesoderm_muscle"
            base = rng.gamma(shape=1.2, scale=1.0, size=n_genes)
            trend = np.zeros(n_genes)
            trend[:20] += t * 0.35
            if lineage == "neural":
                trend[20:45] += 2.2 + 0.2 * t
                trend[-6] += 1.2
                trend[-5] += 1.0
            else:
                trend[45:70] += 2.0 + 0.15 * t
                trend[-7] += 1.3
                trend[-4] += 1.0
            cycle = rng.random() < (0.25 + 0.05 * (lineage == "mesoderm_muscle"))
            if cycle:
                trend[-10:-7] += 2.0
            mu = np.exp(np.log1p(base + trend))
            counts = rng.poisson(np.clip(mu, 0.01, 50.0))
            mats.append(counts)
            rows.append(
                {
                    "cell_id": f"synthetic_t{t}_{i}",
                    "time_point": f"T{time_value:g}",
                    "time_numeric": time_value,
                    "stage_num": time_value,
                    "lineage": lineage,
                    "cell_type": lineage,
                    "condition": "control",
                    "batch": f"batch_{i % 3}",
                    "cell_cycle_score_seed": float(cycle),
                    "is_perturbed": False,
                }
            )
    adata = ad.AnnData(sparse.csr_matrix(np.vstack(mats).astype(np.float32)))
    adata.obs = pd.DataFrame(rows).set_index("cell_id")
    adata.var = pd.DataFrame({"swarm_gene_symbol": genes}, index=genes)
    return adata


def write_synthetic_fixture(path: str | Path, seed: int = 17) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    make_synthetic_lineage_adata(seed=seed).write_h5ad(out)
    return out

