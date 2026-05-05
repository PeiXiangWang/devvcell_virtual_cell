"""Common basis fitting interface."""

from __future__ import annotations

import numpy as np

from devspectrum.dct import dct_feature_row, fit_dct
from devspectrum.spline import fit_spline, spline_feature_row
from devspectrum.wavelet import fit_wavelet, wavelet_feature_row


def fit_basis(times: np.ndarray, values: np.ndarray, method: str, config: dict | None = None):
    config = config or {}
    if method == "dct":
        return fit_dct(times, values, max_basis=int(config.get("max_basis", 4)))
    if method == "wavelet":
        return fit_wavelet(times, values, max_level=config.get("max_level", 2))
    if method == "spline":
        return fit_spline(times, values, degree=int(config.get("degree", 3)))
    raise ValueError(f"Unknown basis method: {method}")


def basis_feature_row(times: np.ndarray, values: np.ndarray, method: str, config: dict | None = None) -> tuple[dict, object]:
    config = config or {}
    if method == "dct":
        return dct_feature_row(times, values, max_basis=int(config.get("max_basis", 4)))
    if method == "wavelet":
        return wavelet_feature_row(times, values, max_level=config.get("max_level", 2))
    if method == "spline":
        return spline_feature_row(times, values, degree=int(config.get("degree", 3)))
    raise ValueError(f"Unknown basis method: {method}")
