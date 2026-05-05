"""Compute TS14-TS19 enrichment statistics for DevVCell response-recovery."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, mannwhitneyu


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from devvcell.io import load_json, resolve_project_path  # noqa: E402
from devvcell.tables import write_table  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/response_recovery.json")
    return parser.parse_args()


def safe_mannwhitney(window: pd.Series, outside: pd.Series) -> tuple[float, float]:
    x = pd.to_numeric(window, errors="coerce").dropna()
    y = pd.to_numeric(outside, errors="coerce").dropna()
    if len(x) == 0 or len(y) == 0:
        return float("nan"), float("nan")
    stat, pvalue = mannwhitneyu(x, y, alternative="two-sided")
    return float(stat), float(pvalue)


def odds_ratio(a: int, b: int, c: int, d: int) -> float:
    return float(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))


def main() -> None:
    args = parse_args()
    cfg = load_json(args.config)
    classes_path = resolve_project_path(cfg["output"]["response_recovery_classes"])
    classes = pd.read_csv(classes_path)
    window_cfg = cfg["stage_window_of_interest"]
    in_window = classes["stage_num"].astype(int).between(int(window_cfg["min_stage"]), int(window_cfg["max_stage"]))
    classes = classes.copy()
    classes["window"] = np.where(in_window, window_cfg["name"], "outside")

    rows: list[dict[str, object]] = []
    for metric in ["response_amplitude", "recovery_cost", "developmental_delay_score", "fate_deflection_index", "off_manifold_score"]:
        stat, pvalue = safe_mannwhitney(classes.loc[in_window, metric], classes.loc[~in_window, metric])
        rows.append(
            {
                "test": "mannwhitneyu",
                "target": metric,
                "window": window_cfg["name"],
                "window_n": int(in_window.sum()),
                "outside_n": int((~in_window).sum()),
                "window_mean": float(pd.to_numeric(classes.loc[in_window, metric], errors="coerce").mean()),
                "outside_mean": float(pd.to_numeric(classes.loc[~in_window, metric], errors="coerce").mean()),
                "statistic": stat,
                "p_value": pvalue,
                "odds_ratio": np.nan,
            }
        )

    for klass in sorted(classes["response_recovery_class"].dropna().unique()):
        a = int(((classes["response_recovery_class"] == klass) & in_window).sum())
        b = int(((classes["response_recovery_class"] != klass) & in_window).sum())
        c = int(((classes["response_recovery_class"] == klass) & ~in_window).sum())
        d = int(((classes["response_recovery_class"] != klass) & ~in_window).sum())
        table = np.array([[a, b], [c, d]])
        chi2, pvalue, _, _ = chi2_contingency(table, correction=False)
        rows.append(
            {
                "test": "chi2_class_enrichment",
                "target": klass,
                "window": window_cfg["name"],
                "window_n": int(in_window.sum()),
                "outside_n": int((~in_window).sum()),
                "window_mean": a / max(a + b, 1),
                "outside_mean": c / max(c + d, 1),
                "statistic": float(chi2),
                "p_value": float(pvalue),
                "odds_ratio": odds_ratio(a, b, c, d),
            }
        )

    output = write_table(pd.DataFrame(rows), cfg["output"]["window_enrichment_statistics"])
    print(f"Wrote window enrichment statistics: {output}")


if __name__ == "__main__":
    main()
