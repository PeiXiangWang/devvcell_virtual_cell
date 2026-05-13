from pathlib import Path

from src.discovery.common import configure_paths


def test_run_all_forwards_quick_flag_to_plotting_and_manifest():
    ps1 = Path("reproducibility/run_all.ps1").read_text(encoding="utf-8")
    sh = Path("reproducibility/run_all.sh").read_text(encoding="utf-8")
    assert "python -m src.plotting.figures_main @quick" in ps1
    assert "python -m src.plotting.figures_extended @quick" in ps1
    assert "python -m src.utils.reproducibility @quick" in ps1
    assert "python -m src.plotting.figures_main ${QUICK_FLAG}" in sh
    assert "python -m src.plotting.figures_extended ${QUICK_FLAG}" in sh
    assert "python -m src.utils.reproducibility ${QUICK_FLAG}" in sh


def test_quick_manifest_uses_separate_output_path():
    text = Path("src/utils/reproducibility.py").read_text(encoding="utf-8")
    assert "manifest.quick_fixture.json" in text
    assert "reports/quick_fixture/install_failures.md" in text


def test_quick_evaluate_reports_do_not_target_main_reports():
    train_cfg, model_cfg, _ = configure_paths("configs/train.yaml", quick_fixture=True)
    assert train_cfg["module_contribution_report"] == "reports/quick_fixture/module_contribution_audit.md"
    assert train_cfg["scientific_gap_report"] == "reports/quick_fixture/scientific_gap_audit.md"
    assert model_cfg["metrics_path"] == "tables/quick_fixture/final_metrics.csv"
