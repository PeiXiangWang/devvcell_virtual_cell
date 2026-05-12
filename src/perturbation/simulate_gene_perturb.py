from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from src.utils.config import ensure_dir, load_config, write_text


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/model.yaml")
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    _ = load_config(args.config)
    out_dir = "tables/quick_fixture" if args.quick_fixture else "tables"
    report_dir = "reports/quick_fixture" if args.quick_fixture else "reports"
    ensure_dir(out_dir)
    ensure_dir(report_dir)
    candidates = pd.DataFrame(
        [
            {"gene_or_axis": "Mki67/Top2a proliferation module", "expected_effect": "growth hazard calibration", "evidence": "cell-cycle marker proxy in current AnnData"},
            {"gene_or_axis": "Fgf8-Fgfr1", "expected_effect": "mesoderm/neural fate-balance shift", "evidence": "curated LR axis if genes are present"},
            {"gene_or_axis": "Cxcl12-Cxcr4", "expected_effect": "dispersion/niche migration shift", "evidence": "curated LR axis if genes are present"},
        ]
    )
    out_path = f"{out_dir}/gene_perturbation_candidates.csv"
    candidates.to_csv(out_path, index=False)
    write_text(
        f"{report_dir}/gene_perturbation_report.md",
        "# Gene Perturbation Simulation\n\nNo matched perturbation time-series is used in the default lineage prototype. Candidate perturbations are therefore hypotheses for future validation, not claims of predictive accuracy.\n\n"
        + candidates.to_markdown(index=False)
        + "\n",
    )
    print({"output": out_path})


if __name__ == "__main__":
    main()
