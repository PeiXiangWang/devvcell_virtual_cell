from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from src.utils.config import ensure_dir, load_config, write_text


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/model.yaml")
    args = parser.parse_args()
    _ = load_config(args.config)
    ensure_dir("tables")
    ensure_dir("reports")
    candidates = pd.DataFrame(
        [
            {"gene_or_axis": "Mki67/Top2a proliferation module", "expected_effect": "growth hazard calibration", "evidence": "cell-cycle marker proxy in current AnnData"},
            {"gene_or_axis": "Fgf8-Fgfr1", "expected_effect": "mesoderm/neural fate-balance shift", "evidence": "curated LR axis if genes are present"},
            {"gene_or_axis": "Cxcl12-Cxcr4", "expected_effect": "dispersion/niche migration shift", "evidence": "curated LR axis if genes are present"},
        ]
    )
    candidates.to_csv("tables/gene_perturbation_candidates.csv", index=False)
    write_text(
        "reports/gene_perturbation_report.md",
        "# Gene Perturbation Simulation\n\nNo matched perturbation time-series is used in the default lineage prototype. Candidate perturbations are therefore hypotheses for future validation, not claims of predictive accuracy.\n\n"
        + candidates.to_markdown(index=False)
        + "\n",
    )
    print({"output": "tables/gene_perturbation_candidates.csv"})


if __name__ == "__main__":
    main()

