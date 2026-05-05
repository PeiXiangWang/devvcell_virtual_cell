from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from devguard.io import ensure_dir, write_dataframe, write_manifest


def _bh_fdr(p_values: list[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    q = np.empty_like(ranked)
    m = len(ranked)
    running = 1.0
    for i in range(m - 1, -1, -1):
        running = min(running, ranked[i] * m / (i + 1))
        q[i] = running
    out = np.empty_like(q)
    out[order] = q
    return out.tolist()


def _class_set(frame: pd.DataFrame) -> list[str]:
    preferred = [
        "within_stage_normal",
        "developmental_delay",
        "developmental_acceleration",
        "fate_deviation",
        "abnormal_off_normal",
    ]
    observed = set(frame["normality_class"].astype(str))
    return [cls for cls in preferred if cls in observed] + sorted(observed - set(preferred))


def compare_perturbation_to_control(
    perturbation_csv: str | Path,
    control_csv: str | Path,
    output_dir: str | Path,
    *,
    label: str,
) -> Path:
    pert = pd.read_csv(perturbation_csv)
    ctrl = pd.read_csv(control_csv)
    rows = []
    classes = _class_set(pd.concat([pert[["normality_class"]], ctrl[["normality_class"]]], axis=0))
    for cls in classes:
        pert_in = int((pert["normality_class"] == cls).sum())
        ctrl_in = int((ctrl["normality_class"] == cls).sum())
        pert_total = int(pert.shape[0])
        ctrl_total = int(ctrl.shape[0])
        table = [[pert_in, pert_total - pert_in], [ctrl_in, ctrl_total - ctrl_in]]
        odds_ratio, p_value = fisher_exact(table)
        rows.append(
            {
                "comparison": label,
                "normality_class": cls,
                "perturbation_cells_in_class": pert_in,
                "perturbation_total_cells": pert_total,
                "perturbation_fraction": pert_in / pert_total if pert_total else 0.0,
                "control_cells_in_class": ctrl_in,
                "control_total_cells": ctrl_total,
                "control_fraction": ctrl_in / ctrl_total if ctrl_total else 0.0,
                "odds_ratio": float(odds_ratio),
                "fisher_p_value": float(p_value),
            }
        )
    frame = pd.DataFrame(rows)
    frame["fdr_q_value"] = _bh_fdr(frame["fisher_p_value"].tolist()) if not frame.empty else []
    output = ensure_dir(output_dir)
    out_path = output / "perturbation_vs_control_class_enrichment.csv"
    write_dataframe(frame, out_path)
    write_manifest(
        output / "perturbation_vs_control_class_enrichment_manifest.json",
        name="compare_perturbation_to_control_baseline",
        inputs=[str(perturbation_csv), str(control_csv)],
        outputs=[str(out_path)],
        parameters={"label": label},
        metrics={"n_tests": int(frame.shape[0])},
    )
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare DevGuard perturbation classes against heldout control baseline.")
    parser.add_argument("--perturbation-csv", required=True)
    parser.add_argument("--control-csv", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--label", required=True)
    args = parser.parse_args()
    compare_perturbation_to_control(
        args.perturbation_csv,
        args.control_csv,
        args.output_dir,
        label=args.label,
    )


if __name__ == "__main__":
    main()
