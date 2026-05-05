import pandas as pd

from devspectrum.devguard_link import correlate_dti_with_spectral_residuals, failure_mode_spectral_signature
from devspectrum.rescue import rank_rescue_candidates


def test_devguard_link_and_rescue_candidate_tables():
    residuals = pd.DataFrame(
        {
            "cohort": ["Tal1_chimera"] * 4,
            "lineage": ["a", "a", "b", "b"],
            "module_name": ["m1", "m2", "m1", "m2"],
            "normality_class": ["abnormal_off_normal", "fate_deviation"] * 2,
            "spectral_residual": [2.0, 0.5, 1.0, -0.2],
            "absolute_spectral_residual": [2.0, 0.5, 1.0, 0.2],
            "n_cells": [20, 20, 20, 20],
        }
    )
    signature = failure_mode_spectral_signature(residuals)
    assert not signature.empty
    candidates, report = rank_rescue_candidates(residuals)
    assert candidates.iloc[0]["module_name"] == "m1"
    assert "in silico" in report

    residual_summary = pd.DataFrame(
        {
            "cohort": ["Tal1_chimera", "Tal1_chimera", "Tal1_chimera"],
            "lineage": ["a", "b", "c"],
            "global_spectral_distance": [3.0, 2.0, 1.0],
            "mean_absolute_residual": [2.0, 1.0, 0.5],
            "max_absolute_residual": [3.0, 2.0, 1.0],
        }
    )
    dti = pd.DataFrame({"lineage": ["a", "b", "c"], "DTI": [-1.0, 0.0, 1.0]})
    corr = correlate_dti_with_spectral_residuals(residual_summary, {"Tal1_chimera": dti})
    assert corr["n_lineages"].min() == 3
