from devspectrum.reconstruction import leave_one_timepoint_predictions, summarize_reconstruction
from devspectrum.timeseries import make_quick_timeseries


def test_leave_one_timepoint_reconstruction_outputs_metrics():
    timeseries = make_quick_timeseries()
    predictions = leave_one_timepoint_predictions(
        timeseries,
        methods=["mean", "linear", "dct", "wavelet", "dct_wavelet"],
        basis_config={"dct": {"max_basis": 3}, "wavelet": {"max_level": 2}},
    )
    metrics, summary = summarize_reconstruction(predictions)
    assert not predictions.empty
    assert {"mse", "mae", "pearson", "spearman"}.issubset(metrics.columns)
    assert set(summary["method"]).issuperset({"mean", "linear", "dct"})
