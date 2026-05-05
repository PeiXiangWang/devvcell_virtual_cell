import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from devguard.markers import log_normalized_expression, score_marker_modules
from devguard.robustness import balanced_downsample_class_fractions, bootstrap_sample_class_fractions


def test_marker_scoring_resolves_symbols_from_var_columns():
    matrix = sparse.csr_matrix([[10, 0, 5], [0, 20, 5]], dtype=float)
    obs = pd.DataFrame(index=["cell1", "cell2"])
    var = pd.DataFrame({"SYMBOL": ["Tal1", "Kdr", "Sox2"]}, index=["ENSMUSG1", "ENSMUSG2", "ENSMUSG3"])
    adata = ad.AnnData(X=matrix, obs=obs, var=var)

    expression, presence = log_normalized_expression(adata, ["Tal1", "Kdr", "Missing"])
    assert expression.columns.tolist() == ["Tal1", "Kdr"]
    assert not bool(presence.set_index("requested_gene").loc["Missing", "available"])

    scores, module_table = score_marker_modules(adata, {"vascular": ["Tal1", "Kdr"]})
    assert "vascular" in scores.columns
    assert module_table.loc[0, "n_available_genes"] == 2


def test_sample_bootstrap_and_balanced_downsample_return_class_fractions():
    frame = pd.DataFrame(
        {
            "perturbation_name": ["A"] * 4 + ["B"] * 4,
            "sample_id": ["a1", "a1", "a2", "a2", "b1", "b1", "b2", "b2"],
            "lineage": ["x", "y", "x", "y"] * 2,
            "normality_class": [
                "within_stage_normal",
                "fate_deviation",
                "within_stage_normal",
                "abnormal_off_normal",
                "within_stage_normal",
                "within_stage_normal",
                "fate_deviation",
                "abnormal_off_normal",
            ],
        }
    )
    boot = bootstrap_sample_class_fractions(frame, n_bootstrap=10, seed=1)
    assert set(boot["perturbation_name"]) == {"A", "B"}
    assert boot["bootstrap_ci_lower"].le(boot["bootstrap_ci_upper"]).all()

    downsampled = balanced_downsample_class_fractions(
        frame,
        strata_cols=["lineage"],
        n_iterations=5,
        max_cells_per_stratum=1,
        seed=1,
    )
    assert downsampled["n_strata"].max() == 2
    assert downsampled.groupby(["iteration", "perturbation_name"])["fraction"].sum().round(6).eq(1.0).all()
