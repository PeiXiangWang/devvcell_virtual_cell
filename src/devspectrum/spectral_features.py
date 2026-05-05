"""Spectral feature helpers."""

from __future__ import annotations

import numpy as np


def energy(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sum(values * values))


def spectral_entropy(coefficients: np.ndarray, *, skip_intercept: bool = True) -> float:
    coefficients = np.asarray(coefficients, dtype=float)
    if skip_intercept and coefficients.size > 1:
        coefficients = coefficients[1:]
    powers = coefficients * coefficients
    total = powers.sum()
    if total <= 0:
        return 0.0
    probabilities = powers / total
    entropy = -np.sum(probabilities * np.log2(np.maximum(probabilities, 1e-12)))
    max_entropy = np.log2(max(probabilities.size, 2))
    return float(entropy / max_entropy)


def dominant_basis_index(coefficients: np.ndarray, *, skip_intercept: bool = True) -> int:
    coefficients = np.asarray(coefficients, dtype=float)
    start = 1 if skip_intercept and coefficients.size > 1 else 0
    if coefficients.size <= start:
        return 0
    return int(start + np.argmax(np.abs(coefficients[start:])))


def phase_shift_proxy(times: np.ndarray, values: np.ndarray) -> float:
    """Return an expression-weighted time center as a simple phase proxy."""

    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    if values.size == 0 or np.all(~np.isfinite(values)):
        return float(np.nanmean(times))
    shifted = values - np.nanmin(values)
    weights = shifted + 1e-9
    if np.all(~np.isfinite(weights)) or weights.sum() <= 0:
        return float(np.nanmean(times))
    return float(np.sum(times * weights) / np.sum(weights))


def reconstruction_error(observed: np.ndarray, reconstructed: np.ndarray) -> float:
    observed = np.asarray(observed, dtype=float)
    reconstructed = np.asarray(reconstructed, dtype=float)
    return float(np.mean((observed - reconstructed) ** 2))


def amplitude(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.size == 0 or np.all(~np.isfinite(values)):
        return 0.0
    return float(np.nanmax(values) - np.nanmin(values)) if values.size else 0.0


def spectral_distance(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    n = min(left.size, right.size)
    if n == 0:
        return 0.0
    return float(np.linalg.norm(left[:n] - right[:n]))
