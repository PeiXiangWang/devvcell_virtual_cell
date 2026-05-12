from __future__ import annotations

import numpy as np

from src.data.fixtures import make_synthetic_lineage_adata
from src.model.simulator import VARIANTS
from src.ot_teacher.run_moscot import compute_pair_coupling


def test_ablation_contract_separates_intrinsic_and_teacher_variants():
    variants = {v.name: v for v in VARIANTS}
    assert variants["M1_intrinsic_neural"].model_role == "intrinsic"
    assert variants["M1_intrinsic_neural"].flags.use_teacher is False
    assert variants["M2_ot_teacher_force"].model_role == "teacher"
    assert variants["M2_ot_teacher_force"].flags.use_teacher is True
    assert variants["M0_linear_label_interpolation"].baseline == "linear"
    assert variants["M0b_ot_interpolation"].baseline == "ot"
    assert variants["M9_full_memory"].flags.use_memory is True


def test_toy_teacher_backend_pair_has_valid_transition_rows():
    adata = make_synthetic_lineage_adata(n_times=2, cells_per_time=20, n_genes=40, seed=1)
    rng = np.random.default_rng(1)
    adata.obsm["X_pca"] = rng.normal(size=(adata.n_obs, 6))
    result = compute_pair_coupling(
        z=adata.obsm["X_pca"],
        obs=adata.obs,
        t0=12.0,
        t1=13.0,
        time_key="time_numeric",
        cell_type_key="lineage",
        max_cells=15,
        epsilon=0.1,
        seed=1,
    )
    assert result["method"] == "numpy_sinkhorn"
    row_sums = result["transition"].sum(axis=1)
    assert np.allclose(row_sums, 1.0, atol=1e-4)
    assert result["barycentric"].shape[0] == result["source_indices"].shape[0]

