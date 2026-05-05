from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import pandas as pd

from devguard.io import ensure_dir, load_json, write_manifest
from devguard.plotting import plot_class_summary, plot_dti_heatmap, plot_method_schematic


def build_figures(config_path: str | Path = "config/devguard/figure_plan.json") -> list[Path]:
    config = load_json(config_path)
    output_dir = ensure_dir(config.get("output_dir", "results/devguard/figures"))
    outputs: list[Path] = []
    method_path = Path(config["figures"][0]["output"])
    outputs.append(plot_method_schematic(method_path))
    condition_summary = Path(config["condition_summary_csv"])
    if condition_summary.exists():
        outputs.append(plot_class_summary(pd.read_csv(condition_summary), config["figures"][1]["output"]))
    dti_csv = Path(config["dti_csv"])
    if dti_csv.exists():
        outputs.append(plot_dti_heatmap(pd.read_csv(dti_csv), config["figures"][2]["output"]))
    write_manifest(
        output_dir / "figure_manifest.json",
        name="build_devguard_figures",
        inputs=[config.get("condition_summary_csv", ""), config.get("dti_csv", "")],
        outputs=[str(path) for path in outputs],
        parameters=config,
    )
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description="Build DevGuard figure drafts.")
    parser.add_argument("--config", default="config/devguard/figure_plan.json")
    args = parser.parse_args()
    build_figures(args.config)


if __name__ == "__main__":
    main()
