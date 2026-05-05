import numpy as np

from devspectrum.spectral_features import dominant_basis_index, phase_shift_proxy, spectral_distance, spectral_entropy


def test_spectral_feature_helpers_are_bounded_and_directional():
    coeff = np.asarray([1.0, 2.0, 0.0, 0.5])
    assert dominant_basis_index(coeff) == 1
    assert 0.0 <= spectral_entropy(coeff) <= 1.0
    assert spectral_distance(np.asarray([1, 2]), np.asarray([1, 4])) == 2.0
    early = phase_shift_proxy(np.asarray([1, 2, 3]), np.asarray([3, 1, 0]))
    late = phase_shift_proxy(np.asarray([1, 2, 3]), np.asarray([0, 1, 3]))
    assert early < late
