from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import joblib

from devguard.io import ensure_dir, load_json, write_dataframe, write_manifest
from devguard.normality import quality_frame, thresholds_frame


def calibrate_boundaries(config_path: str | Path) -> None:
    config = load_json(config_path)
    model_path = Path(config["model_path"])
    output_dir = ensure_dir(config.get("output_dir", "results/devguard/normality_reference"))
    model = joblib.load(model_path)
    groups = model["groups"]
    alpha = float(config.get("alpha", model.get("alpha", 0.05)))
    thresholds = thresholds_frame(groups, alpha)
    quality = quality_frame(groups, alpha)
    threshold_path = output_dir / "conformal_thresholds.csv"
    quality_path = output_dir / "normality_reference_quality.csv"
    write_dataframe(thresholds, threshold_path)
    write_dataframe(quality, quality_path)
    write_manifest(
        output_dir / "conformal_calibration_manifest.json",
        name="calibrate_conformal_boundaries",
        inputs=[str(model_path)],
        outputs=[str(threshold_path), str(quality_path)],
        parameters=config,
        metrics={"max_heldout_control_fpr": float(quality["heldout_control_fpr"].max())},
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Calibrate DevGuard conformal boundaries.")
    parser.add_argument("--config", default="config/devguard/conformal_thresholds.json")
    args = parser.parse_args()
    calibrate_boundaries(args.config)


if __name__ == "__main__":
    main()
