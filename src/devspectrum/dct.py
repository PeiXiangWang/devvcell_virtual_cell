"""Continuous cosine-basis regression for short developmental time courses."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from devspectrum.spectral_features import (
    amplitude,
    dominant_basis_index,
    energy,
    phase_shift_proxy,
    reconstruction_error,
    spectral_entropy,
)


def scale_times(times: np.ndarray, *, minimum: float | None = None, maximum: float | None = None) -> tuple[np.ndarray, float, float]:
    times = np.asarray(times, dtype=float)
    lo = float(np.nanmin(times) if minimum is None else minimum)
    hi = float(np.nanmax(times) if maximum is None else maximum)
    denom = max(hi - lo, 1e-12)
    return (times - lo) / denom, lo, hi


def cosine_design(scaled_times: np.ndarray, max_basis: int) -> np.ndarray:
    scaled_times = np.asarray(scaled_times, dtype=float)
    columns = [np.ones_like(scaled_times)]
    for k in range(1, max_basis + 1):
        columns.append(np.cos(np.pi * k * scaled_times))
    return np.column_stack(columns)


@dataclass
class DCTFit:
    coefficients: np.ndarray
    time_min: float
    time_max: float
    max_basis: int

    def predict(self, times: np.ndarray) -> np.ndarray:
        scaled, _, _ = scale_times(times, minimum=self.time_min, maximum=self.time_max)
        return cosine_design(scaled, self.max_basis) @ self.coefficients


def fit_dct(times: np.ndarray, values: np.ndarray, *, max_basis: int = 4) -> DCTFit:
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    order = np.argsort(times)
    times = times[order]
    values = values[order]
    basis = min(int(max_basis), max(0, times.size - 1))
    scaled, lo, hi = scale_times(times)
    design = cosine_design(scaled, basis)
    coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
    return DCTFit(coefficients=coefficients, time_min=lo, time_max=hi, max_basis=basis)


def dct_feature_row(times: np.ndarray, values: np.ndarray, *, max_basis: int = 4) -> tuple[dict, DCTFit]:
    fit = fit_dct(times, values, max_basis=max_basis)
    reconstructed = fit.predict(times)
    coeff = fit.coefficients
    low = coeff[1 : min(coeff.size, 3)]
    high = coeff[min(coeff.size, 3) :]
    row = {
        "basis_method": "dct",
        "low_frequency_energy": energy(low),
        "high_frequency_energy": energy(high),
        "spectral_entropy": spectral_entropy(coeff),
        "dominant_basis_index": dominant_basis_index(coeff),
        "phase_shift_proxy": phase_shift_proxy(times, reconstructed),
        "amplitude": amplitude(reconstructed),
        "transient_burst_score": 0.0,
        "reconstruction_error": reconstruction_error(values, reconstructed),
    }
    for index, value in enumerate(coeff):
        row[f"c{index}"] = float(value)
    return row, fit
