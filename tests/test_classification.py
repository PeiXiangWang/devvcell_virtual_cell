from devguard.classification import classify_cell_from_pvalues


def test_classification_priority_current_same_first():
    label = classify_cell_from_pvalues(
        p_current_same=0.2,
        p_early_same=0.9,
        p_late_same=0.9,
        p_other_lineage=0.9,
        p_any_normal=0.9,
        alpha=0.05,
    )
    assert label == "within_stage_normal"


def test_classification_delay_before_acceleration():
    label = classify_cell_from_pvalues(
        p_current_same=0.01,
        p_early_same=0.06,
        p_late_same=0.7,
        p_other_lineage=0.7,
        p_any_normal=0.7,
        alpha=0.05,
    )
    assert label == "developmental_delay"


def test_classification_abnormal_when_no_reference_passes():
    label = classify_cell_from_pvalues(
        p_current_same=0.01,
        p_early_same=0.02,
        p_late_same=0.03,
        p_other_lineage=0.04,
        p_any_normal=0.04,
        alpha=0.05,
    )
    assert label == "abnormal_off_normal"
