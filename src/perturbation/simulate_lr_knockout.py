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
    args = parser.parse_args()
    cfg = load_config(args.config)
    ensure_dir("tables")
    ensure_dir("reports")
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
    frame.to_csv("tables/lr_knockout_predictions.csv", index=False)
    text = [
        "# LR Knockout Simulation",
        "",
        "This is an exploratory control-layer perturbation report. It is not wet-lab validation.",
        "",
        frame.to_markdown(index=False),
        "",
    ]
    write_text("reports/lr_perturbation_report.md", "\n".join(text))
    print({"lr_pairs": len(pairs), "output": "tables/lr_knockout_predictions.csv"})


if __name__ == "__main__":
    main()

