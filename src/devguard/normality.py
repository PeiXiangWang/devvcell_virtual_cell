"""Stage-lineage normality reference construction."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from sklearn.covariance import EmpiricalCovariance
from sklearn.neighbors import NearestNeighbors

from devguard.conformal import conformal_p_values, conformal_threshold, false_positive_rate

SCORE_METHODS = ("knn_distance", "mahalanobis")


def make_group_id(time_point: object, lineage: object) -> str:
    raw = f"{time_point}__{lineage}"
    return (
        str(raw)
        .replace(" ", "_")
        .replace("/", "_")
        .replace("\\", "_")
        .replace(":", "_")
    )


def knn_distance_scores(query_z: np.ndarray, reference_z: np.ndarray, k: int = 15) -> np.ndarray:
    if reference_z.shape[0] == 0:
        raise ValueError("KNN reference set is empty.")
    n_neighbors = min(max(1, k), reference_z.shape[0])
    nn = NearestNeighbors(n_neighbors=n_neighbors)
    nn.fit(reference_z)
    distances, _ = nn.kneighbors(query_z)
    return distances.mean(axis=1)


def mahalanobis_scores(query_z: np.ndarray, reference_z: np.ndarray, regularization: float = 0.01) -> np.ndarray:
    if reference_z.shape[0] < 2:
        raise ValueError("Mahalanobis reference set requires at least two cells.")
    center = reference_z.mean(axis=0)
    cov = EmpiricalCovariance().fit(reference_z).covariance_
    cov = cov + np.eye(cov.shape[0]) * regularization
    inv_cov = np.linalg.pinv(cov)
    delta = query_z - center
    return np.einsum("ij,jk,ik->i", delta, inv_cov, delta)


def score_cells(
    query_z: np.ndarray,
    reference_z: np.ndarray,
    *,
    method: str,
    k: int = 15,
    regularization: float = 0.01,
) -> np.ndarray:
    if method == "knn_distance":
        return knn_distance_scores(query_z, reference_z, k=k)
    if method == "mahalanobis":
        return mahalanobis_scores(query_z, reference_z, regularization=regularization)
    raise ValueError(f"Unknown score method: {method}")


@dataclass
class NormalityGroup:
    group_id: str
    time_point: str
    time_numeric: float
    lineage: str
    train_indices: np.ndarray
    calibration_indices: np.ndarray
    test_indices: np.ndarray
    train_embeddings: np.ndarray
    split_strategy: str = "cell"
    split_unit_column: str = ""
    train_units: list[str] = field(default_factory=list)
    calibration_units: list[str] = field(default_factory=list)
    test_units: list[str] = field(default_factory=list)
    calibration_scores: dict[str, np.ndarray] = field(default_factory=dict)
    test_scores: dict[str, np.ndarray] = field(default_factory=dict)

    @property
    def n_train(self) -> int:
        return int(self.train_indices.size)

    @property
    def n_calibration(self) -> int:
        return int(self.calibration_indices.size)

    @property
    def n_test(self) -> int:
        return int(self.test_indices.size)


def _split_indices(indices: np.ndarray, train_fraction: float, calibration_fraction: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = indices.size
    n_train = max(2, int(round(n * train_fraction)))
    n_calibration = max(1, int(round(n * calibration_fraction)))
    if n_train + n_calibration >= n:
        n_calibration = max(1, min(n_calibration, n - 3))
        n_train = max(2, n - n_calibration - 1)
    n_test = n - n_train - n_calibration
    if n_test < 1:
        raise ValueError("Not enough cells to create train/calibration/test splits.")
    return indices[:n_train], indices[n_train : n_train + n_calibration], indices[n_train + n_calibration :]


def _split_by_units(
    frame: pd.DataFrame,
    rng: np.random.Generator,
    *,
    split_unit_column: str,
    train_fraction: float,
    calibration_fraction: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str], list[str], list[str]]:
    units = frame[split_unit_column].astype("string").fillna("NA").astype(str)
    unique_units = np.asarray(sorted(units.unique()))
    if unique_units.size < 3:
        raise ValueError(
            f"Need at least 3 unique {split_unit_column} values for sample-level split; "
            f"found {unique_units.size}."
        )
    unique_units = rng.permutation(unique_units)
    n_units = unique_units.size
    n_train = max(1, int(round(n_units * train_fraction)))
    n_calibration = max(1, int(round(n_units * calibration_fraction)))
    if n_train + n_calibration >= n_units:
        n_calibration = max(1, min(n_calibration, n_units - 2))
        n_train = max(1, n_units - n_calibration - 1)
    train_units = unique_units[:n_train]
    calibration_units = unique_units[n_train : n_train + n_calibration]
    test_units = unique_units[n_train + n_calibration :]

    def collect(selected_units: np.ndarray) -> np.ndarray:
        selected = frame.loc[units.isin(selected_units), "_devguard_position"].to_numpy()
        return selected.astype(int)

    train_idx = collect(train_units)
    cal_idx = collect(calibration_units)
    test_idx = collect(test_units)
    if min(train_idx.size, cal_idx.size, test_idx.size) == 0:
        raise ValueError(f"Sample-level split for {split_unit_column} produced an empty split.")
    return train_idx, cal_idx, test_idx, list(train_units), list(calibration_units), list(test_units)


def build_normality_groups(
    embeddings: np.ndarray,
    obs: pd.DataFrame,
    *,
    time_column: str = "time_point",
    time_numeric_column: str = "time_numeric",
    lineage_column: str = "lineage",
    condition_column: str = "is_control",
    min_cells_per_group: int = 30,
    train_fraction: float = 0.6,
    calibration_fraction: float = 0.2,
    split_strategy: str = "cell",
    split_unit_column: str = "sample_id",
    allow_cell_split_fallback: bool = True,
    min_units_per_group: int = 0,
    seed: int = 42,
    score_methods: list[str] | tuple[str, ...] = SCORE_METHODS,
    k: int = 15,
    regularization: float = 0.01,
) -> dict[str, NormalityGroup]:
    """Build normality reference groups from control cells only."""

    rng = np.random.default_rng(seed)
    if condition_column in obs.columns:
        control_mask = obs[condition_column].astype(bool).to_numpy()
    else:
        control_mask = obs["condition"].astype(str).str.lower().eq("control").to_numpy()

    groups: dict[str, NormalityGroup] = {}
    control_obs = obs.loc[control_mask].copy()
    control_obs["_devguard_position"] = np.flatnonzero(control_mask)
    grouped = control_obs.groupby([time_column, lineage_column], dropna=False, observed=True)
    for (time_point, lineage), frame in grouped:
        if frame.shape[0] < min_cells_per_group:
            continue
        if split_strategy == "sample" and min_units_per_group > 0:
            if split_unit_column not in frame.columns:
                if not allow_cell_split_fallback:
                    raise ValueError(f"Missing split unit column for sample-level split: {split_unit_column}")
            else:
                n_units = frame[split_unit_column].astype("string").fillna("NA").nunique()
                if n_units < min_units_per_group:
                    continue
        resolved_strategy = split_strategy
        train_units: list[str] = []
        calibration_units: list[str] = []
        test_units: list[str] = []
        if split_strategy == "sample":
            try:
                train_idx, cal_idx, test_idx, train_units, calibration_units, test_units = _split_by_units(
                    frame,
                    rng,
                    split_unit_column=split_unit_column,
                    train_fraction=train_fraction,
                    calibration_fraction=calibration_fraction,
                )
            except ValueError:
                if not allow_cell_split_fallback:
                    raise
                resolved_strategy = "cell_fallback"
                indices = frame.index.to_numpy()
                if not np.issubdtype(indices.dtype, np.integer):
                    indices = obs.index.get_indexer(frame.index)
                indices = rng.permutation(indices.astype(int))
                train_idx, cal_idx, test_idx = _split_indices(indices, train_fraction, calibration_fraction)
        elif split_strategy == "cell":
            indices = frame.index.to_numpy()
            if not np.issubdtype(indices.dtype, np.integer):
                indices = obs.index.get_indexer(frame.index)
            indices = rng.permutation(indices.astype(int))
            train_idx, cal_idx, test_idx = _split_indices(indices, train_fraction, calibration_fraction)
        else:
            raise ValueError("split_strategy must be 'cell' or 'sample'.")
        group_id = make_group_id(time_point, lineage)
        train_z = embeddings[train_idx]
        group = NormalityGroup(
            group_id=group_id,
            time_point=str(time_point),
            time_numeric=float(pd.to_numeric(frame[time_numeric_column], errors="coerce").median()),
            lineage=str(lineage),
            train_indices=train_idx,
            calibration_indices=cal_idx,
            test_indices=test_idx,
            train_embeddings=train_z,
            split_strategy=resolved_strategy,
            split_unit_column=split_unit_column if resolved_strategy.startswith("sample") else "",
            train_units=train_units,
            calibration_units=calibration_units,
            test_units=test_units,
        )
        for method in score_methods:
            group.calibration_scores[method] = score_cells(
                embeddings[cal_idx],
                train_z,
                method=method,
                k=k,
                regularization=regularization,
            )
            group.test_scores[method] = score_cells(
                embeddings[test_idx],
                train_z,
                method=method,
                k=k,
                regularization=regularization,
            )
        groups[group_id] = group
    return groups


def reference_cells_frame(groups: dict[str, NormalityGroup], obs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group in groups.values():
        for split_name, indices in [
            ("train", group.train_indices),
            ("calibration", group.calibration_indices),
            ("test", group.test_indices),
        ]:
            split_obs = obs.iloc[indices].copy()
            split_obs["reference_group"] = group.group_id
            split_obs["reference_split"] = split_name
            rows.append(split_obs)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, axis=0)


def centroids_frame(groups: dict[str, NormalityGroup]) -> pd.DataFrame:
    rows = []
    for group in groups.values():
        center = group.train_embeddings.mean(axis=0)
        row = {
            "reference_group": group.group_id,
            "time_point": group.time_point,
            "time_numeric": group.time_numeric,
            "lineage": group.lineage,
            "n_train": group.n_train,
            "n_calibration": group.n_calibration,
            "n_test": group.n_test,
            "split_strategy": group.split_strategy,
            "split_unit_column": group.split_unit_column,
            "n_train_units": len(group.train_units),
            "n_calibration_units": len(group.calibration_units),
            "n_test_units": len(group.test_units),
        }
        row.update({f"z{i + 1}": value for i, value in enumerate(center)})
        rows.append(row)
    return pd.DataFrame(rows)


def thresholds_frame(groups: dict[str, NormalityGroup], alpha: float) -> pd.DataFrame:
    rows = []
    for group in groups.values():
        for method, scores in group.calibration_scores.items():
            rows.append(
                {
                    "reference_group": group.group_id,
                    "time_point": group.time_point,
                    "time_numeric": group.time_numeric,
                    "lineage": group.lineage,
                    "score_method": method,
                    "alpha": alpha,
                    "threshold": conformal_threshold(scores, alpha),
                    "n_train": group.n_train,
                    "n_calibration": group.n_calibration,
                    "n_test": group.n_test,
                    "split_strategy": group.split_strategy,
                    "split_unit_column": group.split_unit_column,
                    "n_train_units": len(group.train_units),
                    "n_calibration_units": len(group.calibration_units),
                    "n_test_units": len(group.test_units),
                }
            )
    return pd.DataFrame(rows)


def quality_frame(groups: dict[str, NormalityGroup], alpha: float) -> pd.DataFrame:
    rows = []
    for group in groups.values():
        for method, cal_scores in group.calibration_scores.items():
            test_scores = group.test_scores[method]
            p_values = conformal_p_values(cal_scores, test_scores)
            rows.append(
                {
                    "reference_group": group.group_id,
                    "time_point": group.time_point,
                    "time_numeric": group.time_numeric,
                    "lineage": group.lineage,
                    "score_method": method,
                    "alpha": alpha,
                    "n_train": group.n_train,
                    "n_calibration": group.n_calibration,
                    "n_test": group.n_test,
                    "split_strategy": group.split_strategy,
                    "split_unit_column": group.split_unit_column,
                    "n_train_units": len(group.train_units),
                    "n_calibration_units": len(group.calibration_units),
                    "n_test_units": len(group.test_units),
                    "low_heldout_flag": bool(group.n_test < 20 or (group.split_strategy == "sample" and len(group.test_units) < 2)),
                    "heldout_control_fpr": false_positive_rate(cal_scores, test_scores, alpha),
                    "median_heldout_p_value": float(np.median(p_values)) if p_values.size else float("nan"),
                }
            )
    return pd.DataFrame(rows)
