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
    ensure_dir("reports")
    lr = pd.read_csv("tables/lr_knockout_predictions.csv") if Path("tables/lr_knockout_predictions.csv").exists() else pd.DataFrame()
    gene = pd.read_csv("tables/gene_perturbation_candidates.csv") if Path("tables/gene_perturbation_candidates.csv").exists() else pd.DataFrame()
    text = [
        "# Perturbation Evaluation",
        "",
        "No matched perturb-seq or drug time-series validation is part of the default SwarmLineage-OT quick dataset. The available scPerturb files are catalogued in the data audit and can support future expression-response benchmarks, but they do not validate developmental lineage dynamics here.",
        "",
        "## LR Perturbation Candidates",
        "",
        lr.to_markdown(index=False) if not lr.empty else "No LR predictions generated.",
        "",
        "## Gene Perturbation Candidates",
        "",
        gene.to_markdown(index=False) if not gene.empty else "No gene perturbation candidates generated.",
        "",
    ]
    write_text("reports/perturbation_evaluation.md", "\n".join(text))
    print({"report": "reports/perturbation_evaluation.md"})


if __name__ == "__main__":
    main()
