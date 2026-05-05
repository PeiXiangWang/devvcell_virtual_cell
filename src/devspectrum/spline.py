"""Spline and interpolation baselines for short time courses."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import UnivariateSpline, interp1d

from devspectrum.spectral_features import amplitude, phase_shift_proxy, reconstruction_error


@dataclass
class SplineFit:
    times: np.ndarray
    values: np.ndarray
    degree: int
    spline: object | None

    def predict(self, times: np.ndarray) -> np.ndarray:
        query = np.asarray(times, dtype=float)
        if self.spline is not None:
            return np.asarray(self.spline(query), dtype=float)
        interpolator = interp1d(self.times, self.values, bounds_error=False, fill_value="extrapolate")
        return np.asarray(interpolator(query), dtype=float)


def fit_spline(times: np.ndarray, values: np.ndarray, *, degree: int = 3, smoothing: float = 0.0) -> SplineFit:
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    order = np.argsort(times)
    times = times[order]
    values = values[order]
    unique_times, unique_indices = np.unique(times, return_index=True)
    unique_values = values[unique_indices]
    if unique_times.size > max(1, degree):
        spline = UnivariateSpline(unique_times, unique_values, k=min(degree, unique_times.size - 1), s=smoothing)
    else:
        spline = None
    return SplineFit(times=unique_times, values=unique_values, degree=degree, spline=spline)


def spline_feature_row(times: np.ndarray, values: np.ndarray, *, degree: int = 3) -> tuple[dict, SplineFit]:
    fit = fit_spline(times, values, degree=degree)
    reconstructed = fit.predict(times)
    row = {
        "basis_method": "spline",
        "low_frequency_energy": 0.0,
        "high_frequency_energy": 0.0,
        "spectral_entropy": 0.0,
        "dominant_basis_index": 0,
        "phase_shift_proxy": phase_shift_proxy(times, reconstructed),
        "amplitude": amplitude(reconstructed),
        "transient_burst_score": 0.0,
        "reconstruction_error": reconstruction_error(values, reconstructed),
        "spline_reconstruction_error": reconstruction_error(values, reconstructed),
    }
    return row, fit
