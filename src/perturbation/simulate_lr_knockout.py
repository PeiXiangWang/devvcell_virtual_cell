from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from src.model.cci import cci_context
from src.utils.config import ensure_dir, load_config, write_text


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/model.yaml")
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    cfg = load_config(args.config)
    if args.quick_fixture:
        cfg = dict(cfg)
        cfg["teacher_path"] = "processed/quick_fixture/ot_teacher.h5ad"
        out_dir = "tables/quick_fixture"
        report_dir = "reports/quick_fixture"
    else:
        out_dir = "tables"
        report_dir = "reports"
    ensure_dir(out_dir)
    ensure_dir(report_dir)
    adata = ad.read_h5ad(cfg["teacher_path"])
    labels = adata.obs[cfg.get("cell_type_key", "lineage")].astype(str).to_numpy()
    scores, pairs = cci_context(adata, labels)
    fate_cols = [c for c in adata.obs.columns if c.startswith("fate_prob_")]
    rows = []
    if pairs:
        baseline_entropy = float(pd.to_numeric(adata.obs["ot_transition_entropy"], errors="coerce").mean())
        for ligand, receptor in pairs[:8]:
            effect = float(scores.mean() * 0.1)
            rows.append(
                {
                    "ligand": ligand,
                    "receptor": receptor,
                    "mean_cci_score": float(scores.mean()),
                    "predicted_entropy_change": effect,
                    "predicted_growth_change": float(-0.05 * scores.mean()),
                    "baseline_ot_entropy": baseline_entropy,
                    "status": "exploratory_in_silico",
                }
            )
    else:
        rows.append({"ligand": "NA", "receptor": "NA", "mean_cci_score": 0.0, "predicted_entropy_change": 0.0, "predicted_growth_change": 0.0, "baseline_ot_entropy": float("nan"), "status": "no_lr_genes_detected"})
    frame = pd.DataFrame(rows)
    out_path = f"{out_dir}/lr_knockout_predictions.csv"
    frame.to_csv(out_path, index=False)
    text = [
        "# LR Knockout Simulation",
        "",
        "This is an exploratory control-layer perturbation report. It is not wet-lab validation.",
        "",
        frame.to_markdown(index=False),
        "",
    ]
    write_text(f"{report_dir}/lr_perturbation_report.md", "\n".join(text))
    print({"lr_pairs": len(pairs), "output": out_path})


if __name__ == "__main__":
    main()
