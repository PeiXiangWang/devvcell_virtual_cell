"""Small Haar-wavelet utilities without external wavelet dependencies."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from devspectrum.spectral_features import amplitude, energy, phase_shift_proxy, reconstruction_error, spectral_entropy


def _next_power_of_two(n: int) -> int:
    power = 1
    while power < n:
        power *= 2
    return power


def _haar_decompose(values: np.ndarray) -> tuple[float, list[np.ndarray]]:
    current = np.asarray(values, dtype=float)
    details: list[np.ndarray] = []
    while current.size > 1:
        avg = (current[0::2] + current[1::2]) / 2.0
        diff = (current[0::2] - current[1::2]) / 2.0
        details.append(diff)
        current = avg
    return float(current[0]), details


def _haar_reconstruct(approx: float, details: list[np.ndarray]) -> np.ndarray:
    current = np.asarray([approx], dtype=float)
    for detail in reversed(details):
        up = np.empty(detail.size * 2, dtype=float)
        up[0::2] = current + detail
        up[1::2] = current - detail
        current = up
    return current


@dataclass
class WaveletFit:
    sorted_times: np.ndarray
    padded_values: np.ndarray
    approx: float
    details: list[np.ndarray]
    original_n: int

    def reconstruct_observed_grid(self) -> np.ndarray:
        return _haar_reconstruct(self.approx, self.details)[: self.original_n]

    def predict(self, times: np.ndarray) -> np.ndarray:
        reconstructed = self.reconstruct_observed_grid()
        return np.interp(np.asarray(times, dtype=float), self.sorted_times, reconstructed)


def fit_wavelet(times: np.ndarray, values: np.ndarray, *, max_level: int | None = None) -> WaveletFit:
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    order = np.argsort(times)
    sorted_times = times[order]
    sorted_values = values[order]
    padded_n = _next_power_of_two(sorted_values.size)
    if padded_n > sorted_values.size:
        pad_values = np.repeat(sorted_values[-1], padded_n - sorted_values.size)
        padded = np.concatenate([sorted_values, pad_values])
    else:
        padded = sorted_values.copy()
    approx, details = _haar_decompose(padded)
    if max_level is not None:
        keep = max(0, int(max_level))
        details = [detail if i < keep else np.zeros_like(detail) for i, detail in enumerate(details)]
    return WaveletFit(sorted_times=sorted_times, padded_values=padded, approx=approx, details=details, original_n=sorted_values.size)


def wavelet_feature_row(times: np.ndarray, values: np.ndarray, *, max_level: int | None = 2) -> tuple[dict, WaveletFit]:
    fit = fit_wavelet(times, values, max_level=max_level)
    reconstructed = fit.predict(times)
    detail_energies = np.asarray([energy(detail) for detail in fit.details], dtype=float)
    all_coefficients = np.concatenate([[fit.approx], *fit.details]) if fit.details else np.asarray([fit.approx])
    high = float(detail_energies[0]) if detail_energies.size else 0.0
    low = float(detail_energies[1:].sum()) if detail_energies.size > 1 else 0.0
    burst = float(max(np.max(np.abs(detail)) for detail in fit.details)) if fit.details else 0.0
    local_index = int(np.argmax(np.abs(fit.details[0]))) if fit.details and fit.details[0].size else 0
    row = {
        "basis_method": "wavelet",
        "wavelet_low_energy": low,
        "wavelet_high_energy": high,
        "low_frequency_energy": low,
        "high_frequency_energy": high,
        "spectral_entropy": spectral_entropy(all_coefficients, skip_intercept=False),
        "dominant_basis_index": int(np.argmax(detail_energies)) if detail_energies.size else 0,
        "phase_shift_proxy": phase_shift_proxy(times, reconstructed),
        "amplitude": amplitude(reconstructed),
        "transient_burst_score": burst,
        "dominant_time_localization": local_index,
        "reconstruction_error": reconstruction_error(values, reconstructed),
    }
    return row, fit
