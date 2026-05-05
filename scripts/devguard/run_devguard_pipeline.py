from __future__ import annotations

import argparse

import _bootstrap  # noqa: F401

from build_control_reference import build_control_reference
from build_devguard_figures import build_figures
from calibrate_conformal_boundaries import calibrate_boundaries
from classify_perturbed_cells import classify_perturbed_cells
from compute_developmental_tolerance_index import compute_dti
from control_organoid_heterogeneity import control_organoid_heterogeneity
from prepare_mouse_datasets import prepare_registry
from validate_spatial_fate_deviation import validate_spatial_fate_deviation


def run_pipeline(mode: str) -> None:
    if mode not in {"quick", "stress", "main", "full"}:
        raise ValueError("mode must be one of: quick, stress, main, full")
    if mode == "stress":
        prepare_registry("config/devguard/datasets_mouse.json", stress_fixture=True)
        build_control_reference("config/devguard/normality_model_stress.json")
        classify_perturbed_cells("config/devguard/perturbation_tests_stress.json")
        compute_dti(
            "results/devguard_stress/perturbation_classification/cell_normality_classes.csv",
            "results/devguard_stress/tolerance_index",
        )
        build_figures("config/devguard/figure_plan_stress.json")
        return
    quick = mode == "quick"
    prepare_registry("config/devguard/datasets_mouse.json", quick_fixture=quick)
    build_control_reference("config/devguard/normality_model.json")
    calibrate_boundaries("config/devguard/conformal_thresholds.json")
    classify_perturbed_cells("config/devguard/perturbation_tests.json")
    compute_dti()
    build_figures("config/devguard/figure_plan.json")
    if mode == "full":
        control_organoid_heterogeneity()
        validate_spatial_fate_deviation("results/devguard/perturbation_classification/cell_normality_classes.csv")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the DevGuard pipeline.")
    parser.add_argument("--mode", choices=["quick", "stress", "main", "full"], default="quick")
    args = parser.parse_args()
    run_pipeline(args.mode)


if __name__ == "__main__":
    main()
