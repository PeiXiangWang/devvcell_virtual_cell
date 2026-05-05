"""Conformal p-values and calibration summaries."""

from __future__ import annotations

import numpy as np


def conformal_p_value(calibration_scores: np.ndarray, score: float) -> float:
    calibration = np.asarray(calibration_scores, dtype=float)
    if calibration.size == 0:
        raise ValueError("At least one calibration score is required.")
    return float((1 + np.sum(calibration >= score)) / (1 + calibration.size))


def conformal_p_values(calibration_scores: np.ndarray, scores: np.ndarray) -> np.ndarray:
    calibration = np.asarray(calibration_scores, dtype=float)
    query = np.asarray(scores, dtype=float)
    if calibration.size == 0:
        raise ValueError("At least one calibration score is required.")
    return (1 + (calibration[None, :] >= query[:, None]).sum(axis=1)) / (1 + calibration.size)


def conformal_threshold(calibration_scores: np.ndarray, alpha: float) -> float:
    """Largest accepted score under the conformal p-value rule."""

    calibration = np.sort(np.asarray(calibration_scores, dtype=float))
    if calibration.size == 0:
        raise ValueError("At least one calibration score is required.")
    accepted = [score for score in calibration if conformal_p_value(calibration, float(score)) >= alpha]
    if not accepted:
        return float(calibration[0])
    return float(max(accepted))


def false_positive_rate(calibration_scores: np.ndarray, test_scores: np.ndarray, alpha: float) -> float:
    p_values = conformal_p_values(calibration_scores, test_scores)
    return float(np.mean(p_values < alpha)) if p_values.size else float("nan")
