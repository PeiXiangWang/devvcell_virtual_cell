"""Response-recovery metrics and classification rules."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from devvcell.stages import stage_number
from devvcell.tables import latent_columns, numeric_matrix


RESPONSE_CLASSES = (
    "reversible_response",
    "developmental_delay",
    "fate_deflection",
    "off_manifold_collapse",
)


@dataclass(frozen=True)
class ClassificationThresholds:
    delay_margin: float = 0.05
    fate_margin: float = 0.05
    off_manifold_quantile: float = 0.95
    off_manifold_multiplier: float = 1.25
    recovery_cost_high_quantile: float = 0.75


def _thresholds(config: dict | None) -> ClassificationThresholds:
    cfg = (config or {}).get("classification_thresholds", config or {})
    return ClassificationThresholds(
        delay_margin=float(cfg.get("delay_margin", 0.05)),
        fate_margin=float(cfg.get("fate_margin", 0.05)),
        off_manifold_quantile=float(cfg.get("off_manifold_quantile", 0.95)),
        off_manifold_multiplier=float(cfg.get("off_manifold_multiplier", 1.25)),
        recovery_cost_high_quantile=float(cfg.get("recovery_cost_high_quantile", 0.75)),
    )


def _nearest_centroid_distances(values: np.ndarray) -> np.ndarray:
    if len(values) < 2:
        return np.ones(len(values), dtype=float)
    diff = values[:, None, :] - values[None, :, :]
    dist = np.linalg.norm(diff, axis=2)
    dist[dist == 0] = np.nan
    nearest = np.nanmin(dist, axis=1)
    nearest[~np.isfinite(nearest)] = np.nanmedian(nearest[np.isfinite(nearest)]) if np.isfinite(nearest).any() else 1.0
    return nearest


def _safe_min_distance(point: np.ndarray, pool: np.ndarray) -> float:
    if len(pool) == 0:
        return float("inf")
    return float(np.linalg.norm(pool - point[None, :], axis=1).min())


def classify_from_latent_tables(
    transferred: pd.DataFrame,
    centroids: pd.DataFrame,
    config: dict | None = None,
) -> pd.DataFrame:
    """Classify perturbed stage/cell-type states against a normal manifold.

    Required columns in both tables are ``stage`` and ``cell_type``. Centroid
    latent dimensions use ``latent_`` columns. Transferred responses may either
    contain ``response_latent_`` columns or already contain perturbed
    ``perturbed_latent_`` columns.
    """

    thresholds = _thresholds(config)
    centroid_cols = latent_columns(centroids, prefixes=("latent_",))
    if not centroid_cols:
        raise ValueError("Centroid table must contain latent_* columns.")

    response_cols = latent_columns(transferred, prefixes=("response_latent_",))
    perturbed_cols = latent_columns(transferred, prefixes=("perturbed_latent_",))
    if not response_cols and not perturbed_cols:
        raise ValueError("Transferred table must contain response_latent_* or perturbed_latent_* columns.")

    centroid_keyed = centroids.copy()
    centroid_keyed["stage_num"] = centroid_keyed["stage"].map(stage_number).astype(int)
    centroid_keyed["cell_type"] = centroid_keyed["cell_type"].astype(str)
    centroid_values = numeric_matrix(centroid_keyed, centroid_cols)
    nearest_normal = _nearest_centroid_distances(centroid_values)
    off_threshold = float(np.nanquantile(nearest_normal, thresholds.off_manifold_quantile) * thresholds.off_manifold_multiplier)
    if not np.isfinite(off_threshold) or off_threshold <= 0:
        off_threshold = 1.0

    centroid_lookup = {
        (str(row.stage), str(row.cell_type)): centroid_values[idx]
        for idx, row in centroid_keyed[["stage", "cell_type"]].iterrows()
    }

    rows: list[dict[str, object]] = []
    for _, row in transferred.iterrows():
        stage = str(row["stage"])
        cell_type = str(row["cell_type"])
        key = (stage, cell_type)
        if key not in centroid_lookup:
            continue

        normal = centroid_lookup[key]
        stage_num = int(stage_number(stage))
        if response_cols:
            response = row[response_cols].astype(float).to_numpy()
            perturbed = normal + response
        else:
            perturbed = row[perturbed_cols].astype(float).to_numpy()
            response = perturbed - normal

        same_cell = centroid_keyed["cell_type"].astype(str) == cell_type
        earlier = centroid_keyed[same_cell & (centroid_keyed["stage_num"] < stage_num)]
        future_same = centroid_keyed[same_cell & (centroid_keyed["stage_num"] >= stage_num)]
        future_other = centroid_keyed[(~same_cell) & (centroid_keyed["stage_num"] >= stage_num)]

        current_distance = float(np.linalg.norm(perturbed - normal))
        early_distance = _safe_min_distance(perturbed, numeric_matrix(earlier, centroid_cols)) if len(earlier) else float("inf")
        future_same_distance = _safe_min_distance(perturbed, numeric_matrix(future_same, centroid_cols))
        future_other_distance = _safe_min_distance(perturbed, numeric_matrix(future_other, centroid_cols))
        global_distance = _safe_min_distance(perturbed, centroid_values)

        delay_score = early_distance - current_distance
        fate_deflection_index = future_same_distance - future_other_distance
        recovery_cost = min(current_distance, future_same_distance)

        if global_distance > off_threshold:
            response_class = "off_manifold_collapse"
        elif np.isfinite(delay_score) and delay_score < -thresholds.delay_margin:
            response_class = "developmental_delay"
        elif np.isfinite(fate_deflection_index) and fate_deflection_index > thresholds.fate_margin:
            response_class = "fate_deflection"
        else:
            response_class = "reversible_response"

        rows.append(
            {
                **{col: row[col] for col in transferred.columns if not col.startswith(("response_latent_", "perturbed_latent_"))},
                "stage_num": stage_num,
                "response_amplitude": float(np.linalg.norm(response)),
                "recovery_cost": float(recovery_cost),
                "developmental_delay_score": float(delay_score) if np.isfinite(delay_score) else np.nan,
                "fate_deflection_index": float(fate_deflection_index) if np.isfinite(fate_deflection_index) else np.nan,
                "off_manifold_score": float(global_distance),
                "off_manifold_threshold": off_threshold,
                "response_recovery_class": response_class,
            }
        )

    return pd.DataFrame(rows)


def summarize_response_recovery(classes: pd.DataFrame) -> pd.DataFrame:
    if classes.empty:
        return pd.DataFrame()
    summary = (
        classes.groupby(["stage", "stage_num", "response_recovery_class"], as_index=False)
        .agg(
            n_cases=("response_recovery_class", "size"),
            mean_response_amplitude=("response_amplitude", "mean"),
            mean_recovery_cost=("recovery_cost", "mean"),
            mean_off_manifold_score=("off_manifold_score", "mean"),
        )
        .sort_values(["stage_num", "response_recovery_class"])
    )
    totals = summary.groupby(["stage", "stage_num"], as_index=False)["n_cases"].sum().rename(columns={"n_cases": "stage_cases"})
    summary = summary.merge(totals, on=["stage", "stage_num"], how="left")
    summary["class_fraction"] = summary["n_cases"] / summary["stage_cases"].replace(0, np.nan)
    return summary
