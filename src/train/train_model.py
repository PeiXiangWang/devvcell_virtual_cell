from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch

from src.model.simulator import VARIANTS, simulate_variant, train_velocity_model
from src.utils.config import ensure_dir, load_config, write_json
from src.utils.metrics import composition_rmse, energy_distance, knn_two_sample_accuracy, mmd_rbf, sinkhorn_distance


def _evaluate_prediction(pred: np.ndarray, true: np.ndarray, pred_labels: np.ndarray, true_labels: np.ndarray, cfg: dict, seed: int) -> dict:
    sample = int(cfg.get("metric_sample_size", 700))
    return {
        "sinkhorn": sinkhorn_distance(pred, true, max_n=min(sample, 450), seed=seed),
        "mmd_rbf": mmd_rbf(pred, true, max_n=sample, seed=seed),
        "energy": energy_distance(pred, true, max_n=sample, seed=seed),
        "knn_two_sample_accuracy": knn_two_sample_accuracy(pred, true, max_n=sample, seed=seed),
        "celltype_composition_rmse": composition_rmse(pred_labels, true_labels),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/model.yaml")
    args = parser.parse_args()
    cfg = load_config(args.config)
    ensure_dir("results/swarmlineage")
    ensure_dir(cfg.get("model_dir", "results/swarmlineage/models"))
    ensure_dir(cfg.get("simulation_dir", "results/swarmlineage/simulations"))
    ensure_dir("tables")
    adata = ad.read_h5ad(cfg["teacher_path"])
    rows = []
    train_summaries = []
    for seed in cfg.get("seeds", [7, 17, 23, 42, 99]):
        model, summary = train_velocity_model(adata, cfg, int(seed))
        model_path = Path(cfg.get("model_dir", "results/swarmlineage/models")) / f"velocity_seed{seed}.pt"
        torch.save({"state_dict": model.state_dict(), "summary": summary, "config": cfg}, model_path)
        summary.update({"seed": int(seed), "model_path": str(model_path)})
        train_summaries.append(summary)
        for variant in VARIANTS:
            pred, true, meta = simulate_variant(adata, cfg, model, variant, int(seed))
            metrics = _evaluate_prediction(pred, true, meta["pred_labels"], meta["true_labels"], cfg, int(seed))
            out_npz = Path(cfg.get("simulation_dir", "results/swarmlineage/simulations")) / f"{variant.name}_seed{seed}.npz"
            np.savez_compressed(out_npz, pred=pred.astype(np.float32), true=true.astype(np.float32), pred_labels=meta["pred_labels"], true_labels=meta["true_labels"])
            rows.append(
                {
                    "seed": int(seed),
                    "model": variant.name,
                    "holdout_time": meta["holdout_time"],
                    "source_time": meta["source_time"],
                    "next_time": meta["next_time"],
                    "n_pred": meta["n_pred"],
                    "n_true": meta["n_true"],
                    "simulation_file": str(out_npz),
                    **metrics,
                }
            )
    metrics_frame = pd.DataFrame(rows)
    metrics_path = Path(cfg.get("metrics_path", "tables/final_metrics.csv"))
    metrics_frame.to_csv(metrics_path, index=False)
    write_json("results/swarmlineage/training_summary.json", {"seeds": train_summaries, "metrics_path": str(metrics_path)})
    print(json.dumps({"metrics": str(metrics_path), "rows": int(metrics_frame.shape[0])}, indent=2))


if __name__ == "__main__":
    main()

