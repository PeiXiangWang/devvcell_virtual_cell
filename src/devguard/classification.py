"""Perturbed-cell classification against DevGuard normality groups."""

from __future__ import annotations

import numpy as np
import pandas as pd

from devguard.conformal import conformal_p_value
from devguard.normality import NormalityGroup, score_cells

CLASS_PRIORITY = [
    "within_stage_normal",
    "developmental_delay",
    "developmental_acceleration",
    "fate_deviation",
    "abnormal_off_normal",
]


def classify_cell_from_pvalues(
    *,
    p_current_same: float,
    p_early_same: float,
    p_late_same: float,
    p_other_lineage: float,
    p_any_normal: float,
    alpha: float = 0.05,
) -> str:
    if p_current_same >= alpha:
        return "within_stage_normal"
    if p_early_same >= alpha:
        return "developmental_delay"
    if p_late_same >= alpha:
        return "developmental_acceleration"
    if p_other_lineage >= alpha:
        return "fate_deviation"
    if p_any_normal >= alpha:
        return "fate_deviation"
    return "abnormal_off_normal"


def _best_group(candidates: list[tuple[str, float]]) -> tuple[str, float]:
    if not candidates:
        return "", 0.0
    group_id, value = max(candidates, key=lambda item: item[1])
    return group_id, float(value)


def classify_cells_against_reference(
    embeddings: np.ndarray,
    obs: pd.DataFrame,
    groups: dict[str, NormalityGroup],
    *,
    score_method: str = "knn_distance",
    alpha: float = 0.05,
    k: int = 15,
    regularization: float = 0.01,
) -> pd.DataFrame:
    rows = []
    for row_position, (_, cell_obs) in enumerate(obs.iterrows()):
        time_numeric = float(pd.to_numeric(cell_obs.get("time_numeric"), errors="coerce"))
        lineage = str(cell_obs.get("lineage"))
        z = embeddings[row_position : row_position + 1]
        p_by_group: dict[str, float] = {}
        for group_id, group in groups.items():
            score = score_cells(
                z,
                group.train_embeddings,
                method=score_method,
                k=k,
                regularization=regularization,
            )[0]
            p_by_group[group_id] = conformal_p_value(group.calibration_scores[score_method], float(score))

        current_same = []
        early_same = []
        late_same = []
        other_lineage = []
        for group_id, p_value in p_by_group.items():
            group = groups[group_id]
            same_lineage = group.lineage == lineage
            same_time = np.isclose(group.time_numeric, time_numeric, equal_nan=False)
            if same_lineage and same_time:
                current_same.append((group_id, p_value))
            elif same_lineage and group.time_numeric < time_numeric:
                early_same.append((group_id, p_value))
            elif same_lineage and group.time_numeric > time_numeric:
                late_same.append((group_id, p_value))
            elif not same_lineage:
                other_lineage.append((group_id, p_value))

        current_group, p_current = _best_group(current_same)
        early_group, p_early = _best_group(early_same)
        late_group, p_late = _best_group(late_same)
        other_group, p_other = _best_group(other_lineage)
        any_group, p_any = _best_group(list(p_by_group.items()))
        normality_class = classify_cell_from_pvalues(
            p_current_same=p_current,
            p_early_same=p_early,
            p_late_same=p_late,
            p_other_lineage=p_other,
            p_any_normal=p_any,
            alpha=alpha,
        )
        assigned_lookup = {
            "within_stage_normal": current_group,
            "developmental_delay": early_group,
            "developmental_acceleration": late_group,
            "fate_deviation": other_group or any_group,
            "abnormal_off_normal": any_group,
        }
        out = cell_obs.to_dict()
        out.update(
            {
                "normality_class": normality_class,
                "score_method": score_method,
                "alpha": alpha,
                "p_current_same": p_current,
                "p_early_same": p_early,
                "p_late_same": p_late,
                "p_other_lineage": p_other,
                "p_any_normal": p_any,
                "reference_current_same": current_group,
                "reference_early_same": early_group,
                "reference_late_same": late_group,
                "reference_other_lineage": other_group,
                "assigned_reference_group": assigned_lookup[normality_class],
            }
        )
        rows.append(out)
    return pd.DataFrame(rows)


def summarize_classes(frame: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    counts = frame.groupby(group_cols + ["normality_class"], dropna=False).size().reset_index(name="n_cells")
    totals = counts.groupby(group_cols, dropna=False)["n_cells"].transform("sum")
    counts["fraction"] = counts["n_cells"] / totals
    return counts
