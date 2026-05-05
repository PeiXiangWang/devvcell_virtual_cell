import numpy as np
import pandas as pd

from devguard.normality import build_normality_groups
from devguard.tolerance import bootstrap_dti_ci


def test_normality_groups_can_split_by_sample_id():
    rows = []
    embeddings = []
    for sample_idx in range(6):
        for cell_idx in range(6):
            rows.append(
                {
                    "time_point": "E7.5",
                    "time_numeric": 7.5,
                    "lineage": "mesoderm",
                    "is_control": True,
                    "sample_id": f"embryo_{sample_idx}",
                }
            )
            embeddings.append([sample_idx, cell_idx / 10])
    obs = pd.DataFrame(rows)
    groups = build_normality_groups(
        np.asarray(embeddings, dtype=float),
        obs,
        min_cells_per_group=12,
        split_strategy="sample",
        split_unit_column="sample_id",
        score_methods=["knn_distance"],
        k=3,
    )
    group = groups["E7.5__mesoderm"]
    assert group.split_strategy == "sample"
    assert set(group.train_units).isdisjoint(group.calibration_units)
    assert set(group.train_units).isdisjoint(group.test_units)
    assert set(group.calibration_units).isdisjoint(group.test_units)


def test_sample_split_indices_are_global_positions_with_string_index():
    rows = []
    embeddings = []
    for cell_idx in range(12):
        rows.append(
            {
                "time_point": "E6.5",
                "time_numeric": 6.5,
                "lineage": "epiblast",
                "is_control": True,
                "sample_id": f"other_{cell_idx // 4}",
            }
        )
        embeddings.append([0.0, float(cell_idx)])
    for sample_idx in range(6):
        for cell_idx in range(6):
            rows.append(
                {
                    "time_point": "E7.5",
                    "time_numeric": 7.5,
                    "lineage": "mesoderm",
                    "is_control": True,
                    "sample_id": f"embryo_{sample_idx}",
                }
            )
            embeddings.append([10.0 + sample_idx, cell_idx / 10])
    obs = pd.DataFrame(rows)
    obs.index = [f"cell_{idx}" for idx in range(obs.shape[0])]

    groups = build_normality_groups(
        np.asarray(embeddings, dtype=float),
        obs,
        min_cells_per_group=12,
        split_strategy="sample",
        split_unit_column="sample_id",
        score_methods=["knn_distance"],
        k=3,
    )

    group = groups["E7.5__mesoderm"]
    for units, indices in [
        (group.train_units, group.train_indices),
        (group.calibration_units, group.calibration_indices),
        (group.test_units, group.test_indices),
    ]:
        observed_units = set(obs.iloc[indices]["sample_id"].astype(str))
        assert observed_units == set(units)


def test_sample_split_can_require_minimum_units():
    rows = []
    embeddings = []
    for sample_idx in range(2):
        for cell_idx in range(6):
            rows.append(
                {
                    "time_point": "E6.5",
                    "time_numeric": 6.5,
                    "lineage": "epiblast",
                    "is_control": True,
                    "sample_id": f"low_{sample_idx}",
                }
            )
            embeddings.append([sample_idx, cell_idx / 10])
    for sample_idx in range(6):
        for cell_idx in range(6):
            rows.append(
                {
                    "time_point": "E7.5",
                    "time_numeric": 7.5,
                    "lineage": "mesoderm",
                    "is_control": True,
                    "sample_id": f"embryo_{sample_idx}",
                }
            )
            embeddings.append([10 + sample_idx, cell_idx / 10])

    groups = build_normality_groups(
        np.asarray(embeddings, dtype=float),
        pd.DataFrame(rows),
        min_cells_per_group=12,
        split_strategy="sample",
        split_unit_column="sample_id",
        min_units_per_group=3,
        allow_cell_split_fallback=False,
        score_methods=["knn_distance"],
        k=3,
    )

    assert "E6.5__epiblast" not in groups
    assert groups["E7.5__mesoderm"].split_strategy == "sample"


def test_bootstrap_dti_ci_returns_sample_level_interval():
    frame = pd.DataFrame(
        {
            "time_point": ["E7.5"] * 8,
            "lineage": ["mesoderm"] * 8,
            "perturbation_name": ["FGF8_KO"] * 8,
            "sample_id": ["s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4"],
            "normality_class": [
                "within_stage_normal",
                "within_stage_normal",
                "developmental_delay",
                "developmental_delay",
                "fate_deviation",
                "fate_deviation",
                "abnormal_off_normal",
                "abnormal_off_normal",
            ],
        }
    )
    ci = bootstrap_dti_ci(frame, n_bootstrap=20, seed=1)
    assert ci.loc[0, "bootstrap_level"] == "sample"
    assert ci.loc[0, "DTI_ci_lower"] <= ci.loc[0, "DTI_ci_upper"]
