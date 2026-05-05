import numpy as np

from devguard.conformal import conformal_p_value, conformal_p_values, false_positive_rate


def test_conformal_p_value_uses_upper_tail_nonconformity():
    calibration = np.array([1.0, 2.0, 3.0, 4.0])
    assert conformal_p_value(calibration, 2.0) == 0.8
    assert conformal_p_value(calibration, 5.0) == 0.2


def test_conformal_p_values_vectorized():
    calibration = np.array([1.0, 2.0, 3.0])
    values = conformal_p_values(calibration, np.array([1.0, 3.5]))
    np.testing.assert_allclose(values, np.array([1.0, 0.25]))


def test_false_positive_rate_flags_low_p_values():
    calibration = np.array([1.0, 2.0, 3.0, 4.0])
    test = np.array([1.0, 6.0])
    assert false_positive_rate(calibration, test, alpha=0.3) == 0.5
