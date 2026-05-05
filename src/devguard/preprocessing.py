"""Metadata standardization for DevGuard datasets."""

from __future__ import annotations

import re
from typing import Mapping

import numpy as np
import pandas as pd

OBS_SCHEMA_COLUMNS = [
    "cell_id",
    "dataset_id",
    "species",
    "system",
    "time_point",
    "time_numeric",
    "condition",
    "perturbation_name",
    "perturbation_type",
    "dose",
    "duration",
    "sample_id",
    "batch",
    "cell_type",
    "lineage",
    "is_control",
    "is_perturbed",
]

CONTROL_LABELS = {"control", "normal", "untreated", "vehicle", "wildtype", "wt"}


def parse_time_numeric(value: object) -> float:
    """Extract a numeric developmental time from labels such as E7.5 or day4."""

    if value is None or (isinstance(value, float) and np.isnan(value)):
        return float("nan")
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value)
    match = re.search(r"-?\d+(?:\.\d+)?", str(value))
    if match is None:
        return float("nan")
    return float(match.group(0))


def _coalesce_column(obs: pd.DataFrame, target: str, source: str) -> None:
    if target not in obs.columns and source in obs.columns:
        obs[target] = obs[source]


def standardize_obs(
    adata,
    *,
    dataset_id: str,
    species: str = "Mus musculus",
    system: str = "gastruloid",
    column_map: Mapping[str, str] | None = None,
) -> object:
    """Return a copy with the required DevGuard obs schema.

    ``column_map`` maps DevGuard target names to source obs column names.
    """

    result = adata.copy()
    obs = result.obs.copy()
    for target, source in (column_map or {}).items():
        _coalesce_column(obs, target, source)

    if "cell_id" not in obs.columns:
        obs["cell_id"] = obs.index.astype(str)
    if "dataset_id" not in obs.columns:
        obs["dataset_id"] = dataset_id
    if "species" not in obs.columns:
        obs["species"] = species
    if "system" not in obs.columns:
        obs["system"] = system
    if "time_point" not in obs.columns:
        obs["time_point"] = "unknown_time"
    if "time_numeric" not in obs.columns:
        obs["time_numeric"] = obs["time_point"].map(parse_time_numeric)
    else:
        obs["time_numeric"] = obs["time_numeric"].map(parse_time_numeric)
    if "condition" not in obs.columns:
        obs["condition"] = "control"
    if "perturbation_name" not in obs.columns:
        obs["perturbation_name"] = np.where(
            obs["condition"].astype(str).str.lower().isin(CONTROL_LABELS),
            "control",
            "unspecified_perturbation",
        )
    if "perturbation_type" not in obs.columns:
        obs["perturbation_type"] = np.where(
            obs["perturbation_name"].astype(str).str.lower().eq("control"),
            "none",
            "unspecified",
        )
    for column in ["dose", "duration", "sample_id", "batch", "cell_type", "lineage"]:
        if column not in obs.columns:
            obs[column] = "NA"

    condition_lower = obs["condition"].astype(str).str.lower()
    if "is_control" not in obs.columns:
        obs["is_control"] = condition_lower.isin(CONTROL_LABELS)
    else:
        obs["is_control"] = obs["is_control"].astype(bool)
    if "is_perturbed" not in obs.columns:
        obs["is_perturbed"] = ~obs["is_control"]
    else:
        obs["is_perturbed"] = obs["is_perturbed"].astype(bool)

    for column in OBS_SCHEMA_COLUMNS:
        if column not in obs.columns:
            obs[column] = "NA"

    result.obs = obs
    return result


def validate_obs_schema(obs: pd.DataFrame) -> list[str]:
    """Return missing required DevGuard obs columns."""

    return [column for column in OBS_SCHEMA_COLUMNS if column not in obs.columns]


def metadata_summary(obs: pd.DataFrame) -> pd.DataFrame:
    """Summarize cell counts by dataset, condition, time and lineage."""

    group_cols = ["dataset_id", "species", "system", "condition", "time_point", "lineage"]
    available = [column for column in group_cols if column in obs.columns]
    return obs.groupby(available, dropna=False, observed=True).size().reset_index(name="n_cells")
