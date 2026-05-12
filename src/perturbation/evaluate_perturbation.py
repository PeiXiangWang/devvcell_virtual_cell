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
    report_dir = "reports/quick_fixture" if args.quick_fixture else "reports"
    table_dir = "tables/quick_fixture" if args.quick_fixture else "tables"
    ensure_dir(report_dir)
    lr = pd.read_csv(f"{table_dir}/lr_knockout_predictions.csv") if Path(f"{table_dir}/lr_knockout_predictions.csv").exists() else pd.DataFrame()
    gene = pd.read_csv(f"{table_dir}/gene_perturbation_candidates.csv") if Path(f"{table_dir}/gene_perturbation_candidates.csv").exists() else pd.DataFrame()
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
    write_text(f"{report_dir}/perturbation_evaluation.md", "\n".join(text))
    print({"report": f"{report_dir}/perturbation_evaluation.md"})


if __name__ == "__main__":
    main()
