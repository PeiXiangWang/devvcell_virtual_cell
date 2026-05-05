import numpy as np

from devspectrum.dct import dct_feature_row, fit_dct
from devspectrum.wavelet import wavelet_feature_row


def test_dct_basis_reconstructs_smooth_short_trajectory():
    times = np.asarray([3.0, 3.5, 4.0, 4.5, 5.0])
    values = 1.0 + 0.5 * np.cos(np.pi * (times - times.min()) / (times.max() - times.min()))
    row, fit = dct_feature_row(times, values, max_basis=2)
    predicted = fit.predict(times)
    assert np.mean((values - predicted) ** 2) < 1e-10
    assert row["low_frequency_energy"] > row["high_frequency_energy"]
    assert 0.0 <= row["spectral_entropy"] <= 1.0


def test_wavelet_feature_detects_local_burst():
    times = np.asarray([3.0, 3.5, 4.0, 4.5, 5.0])
    values = np.asarray([0.1, 0.1, 2.0, 0.1, 0.1])
    row, fit = wavelet_feature_row(times, values, max_level=2)
    assert row["transient_burst_score"] > 0.1
    assert fit.predict(np.asarray([4.0])).shape == (1,)
