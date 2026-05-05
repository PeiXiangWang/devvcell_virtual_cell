"""Deterministic software-validation fixtures for DevGuard."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd


def _group_mean(time_index: int, lineage_index: int, n_genes: int) -> np.ndarray:
    mean = np.full(n_genes, 4.0)
    mean[time_index * 10 : time_index * 10 + 10] += 22.0
    mean[30 + lineage_index * 10 : 30 + lineage_index * 10 + 10] += 28.0
    mean[60:70] += (time_index + 1) * 3.0
    return mean


def create_quick_fixture(seed: int = 42) -> ad.AnnData:
    """Create a deterministic mouse fixture for software validation only."""

    _ = np.random.default_rng(seed)
    n_genes = 90
    genes = [f"Gene{i:03d}" for i in range(n_genes)]
    times = [("E6.5", 6.5), ("E7.5", 7.5), ("E8.5", 8.5)]
    lineages = ["epiblast", "mesoderm", "neural"]
    matrices = []
    rows = []

    for time_index, (time_point, time_numeric) in enumerate(times):
        for lineage_index, lineage in enumerate(lineages):
            mean = _group_mean(time_index, lineage_index, n_genes)
            for i in range(150):
                matrices.append(mean.copy())
                rows.append(
                    {
                        "cell_id": f"ctrl_{time_point}_{lineage}_{i}",
                        "dataset_id": "devguard_quick_mouse",
                        "species": "Mus musculus",
                        "system": "gastruloid",
                        "time_point": time_point,
                        "time_numeric": time_numeric,
                        "condition": "control",
                        "perturbation_name": "control",
                        "perturbation_type": "none",
                        "dose": "NA",
                        "duration": "NA",
                        "sample_id": f"organoid_{i % 5}",
                        "batch": "quick_batch",
                        "cell_type": lineage,
                        "lineage": lineage,
                        "is_control": True,
                        "is_perturbed": False,
                    }
                )

    perturbation_specs = [
        ("FGF8_KO", "within_stage_normal_like", 1, 1, "mesoderm"),
        ("FGF8_KO", "delay_like", 0, 1, "mesoderm"),
        ("FGF8_KO", "acceleration_like", 2, 1, "mesoderm"),
        ("FGF8_KO", "fate_deviation_like", 1, 2, "mesoderm"),
        ("FGF8_KO", "abnormal_like", 1, 1, "mesoderm"),
    ]
    for perturbation_name, subtype, ref_time_index, ref_lineage_index, observed_lineage in perturbation_specs:
        for i in range(24):
            if subtype == "abnormal_like":
                mean = np.full(n_genes, 4.0)
                mean[75:90] += 55.0
            else:
                mean = _group_mean(ref_time_index, ref_lineage_index, n_genes)
            matrices.append(mean.copy())
            rows.append(
                {
                    "cell_id": f"pert_{subtype}_{i}",
                    "dataset_id": "devguard_quick_mouse",
                    "species": "Mus musculus",
                    "system": "gastruloid",
                    "time_point": "E7.5",
                    "time_numeric": 7.5,
                    "condition": "perturbation",
                    "perturbation_name": perturbation_name,
                    "perturbation_type": "genetic",
                    "dose": "NA",
                    "duration": "NA",
                    "sample_id": f"perturb_organoid_{i % 3}",
                    "batch": "quick_batch",
                    "cell_type": observed_lineage,
                    "lineage": observed_lineage,
                    "is_control": False,
                    "is_perturbed": True,
                    "quick_fixture_expected_pattern": subtype,
                }
            )

    obs = pd.DataFrame(rows).set_index("cell_id", drop=False)
    return ad.AnnData(X=np.vstack(matrices), obs=obs, var=pd.DataFrame(index=genes))
