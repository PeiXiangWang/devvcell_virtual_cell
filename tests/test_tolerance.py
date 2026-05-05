import pandas as pd

from devguard.tolerance import compute_developmental_tolerance_index


def test_developmental_tolerance_index_formula():
    frame = pd.DataFrame(
        {
            "time_point": ["E7.5"] * 10,
            "lineage": ["mesoderm"] * 10,
            "perturbation_name": ["FGF8_KO"] * 10,
            "normality_class": [
                "within_stage_normal",
                "within_stage_normal",
                "within_stage_normal",
                "within_stage_normal",
                "within_stage_normal",
                "developmental_delay",
                "developmental_delay",
                "developmental_acceleration",
                "fate_deviation",
                "abnormal_off_normal",
            ],
        }
    )
    dti = compute_developmental_tolerance_index(frame)
    assert dti.loc[0, "R_within_stage_normal"] == 0.5
    assert abs(dti.loc[0, "DTI"]) < 1e-12
