from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch

from src.model.simulator import VARIANTS, simulate_variant, train_dynamics_model
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


def _baseline_status(name: str, module: str, timeout: int = 45) -> dict:
    code = f"import {module}; print(getattr({module}, '__version__', 'unknown'))"
    try:
        out = subprocess.check_output([sys.executable, "-c", code], text=True, stderr=subprocess.STDOUT, timeout=timeout)
        return {"baseline": name, "module": module, "available": True, "executed": False, "status": f"importable:{out.strip()}"}
    except subprocess.TimeoutExpired:
        return {"baseline": name, "module": module, "available": False, "executed": False, "status": f"import_timeout>{timeout}s"}
    except Exception as exc:
        return {"baseline": name, "module": module, "available": False, "executed": False, "status": f"{type(exc).__name__}:{exc}"}


def _execution_matrix(metric_rows: list[dict], cfg: dict) -> pd.DataFrame:
    compared = sorted({row["model"] for row in metric_rows})
    rows = []
    for name in compared:
        rows.append({"baseline": name, "available": True, "executed": True, "status": "evaluated_in_this_run"})
    rows += [
        _baseline_status("CellRank2", "cellrank"),
        _baseline_status("TrajectoryNet", "TrajectoryNet", timeout=15),
        _baseline_status("MIOFlow", "mioflow", timeout=15),
        _baseline_status("TIGON", "tigon", timeout=15),
    ]
    frame = pd.DataFrame(rows)
    path = cfg.get("baseline_execution_matrix_path", "reports/baseline_execution_matrix.csv")
    ensure_dir(Path(path).parent)
    frame.to_csv(path, index=False)
    return frame


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/model.yaml")
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    cfg = load_config(args.config)
    if args.quick_fixture:
        cfg = dict(cfg)
        cfg["teacher_path"] = "processed/quick_fixture/ot_teacher.h5ad"
        cfg["model_dir"] = "results/quick_fixture/models"
        cfg["simulation_dir"] = "results/quick_fixture/simulations"
        cfg["metrics_path"] = "tables/quick_fixture/final_metrics.csv"
        cfg["event_log_path"] = "tables/quick_fixture/birth_death_event_log.csv"
        cfg["order_log_path"] = "tables/quick_fixture/rollout_order_parameters.csv"
        cfg["baseline_execution_matrix_path"] = "reports/quick_fixture/baseline_execution_matrix.csv"
        cfg["holdout_time"] = 14.0
        cfg["seeds"] = [7, 17]
        cfg["epochs"] = int(cfg.get("quick_epochs", 8))
        cfg["simulation_cells_per_seed"] = 70
        cfg["metric_sample_size"] = 80
    ensure_dir("results/swarmlineage")
    ensure_dir(cfg.get("model_dir", "results/swarmlineage/models"))
    ensure_dir(cfg.get("simulation_dir", "results/swarmlineage/simulations"))
    ensure_dir(Path(cfg.get("metrics_path", "tables/final_metrics.csv")).parent)
    adata = ad.read_h5ad(cfg["teacher_path"])
    rows: list[dict] = []
    train_summaries = []
    all_events = []
    all_order_rows = []
    for seed in cfg.get("seeds", [7, 17, 23, 42, 99]):
        intrinsic, intrinsic_summary = train_dynamics_model(adata, cfg, int(seed), "intrinsic")
        teacher, teacher_summary = train_dynamics_model(adata, cfg, int(seed), "teacher")
        models = {"intrinsic": intrinsic, "teacher": teacher}
        model_path = Path(cfg.get("model_dir", "results/swarmlineage/models")) / f"dynamics_seed{seed}.pt"
        torch.save(
            {
                "intrinsic_state_dict": intrinsic.state_dict(),
                "teacher_state_dict": teacher.state_dict(),
                "summary": {"intrinsic": intrinsic_summary, "teacher": teacher_summary},
                "config": cfg,
            },
            model_path,
        )
        train_summaries.append({"seed": int(seed), "model_path": str(model_path), "intrinsic": intrinsic_summary, "teacher": teacher_summary})
        for variant in VARIANTS:
            pred, true, meta = simulate_variant(adata, cfg, models, variant, int(seed))
            metrics = _evaluate_prediction(pred, true, meta["pred_labels"], meta["true_labels"], cfg, int(seed))
            rows.append(
                {
                    "seed": int(seed),
                    "model": variant.name,
                    "holdout_time": meta["holdout_time"],
                    "source_time": meta["source_time"],
                    "next_time": meta["next_time"],
                    "n_pred": meta["n_pred"],
                    "n_true": meta["n_true"],
                    "uses_teacher": variant.flags.use_teacher,
                    "uses_birth_death": variant.flags.use_birth_death,
                    "uses_diffusion": variant.flags.use_diffusion,
                    "uses_swarm": variant.flags.use_swarm,
                    "uses_cci": variant.flags.use_cci,
                    "uses_memory": variant.flags.use_memory,
                    **metrics,
                }
            )
            for event in meta["event_rows"]:
                all_events.append({"seed": int(seed), **event})
            for order_row in meta.get("order_rows", []):
                all_order_rows.append(order_row)
    metrics_frame = pd.DataFrame(rows)
    metrics_path = Path(cfg.get("metrics_path", "tables/final_metrics.csv"))
    metrics_frame.to_csv(metrics_path, index=False)
    event_path = Path(cfg.get("event_log_path", "tables/birth_death_event_log.csv"))
    ensure_dir(event_path.parent)
    pd.DataFrame(all_events).to_csv(event_path, index=False)
    order_path = Path(cfg.get("order_log_path", "tables/rollout_order_parameters.csv"))
    ensure_dir(order_path.parent)
    pd.DataFrame(all_order_rows).to_csv(order_path, index=False)
    matrix = _execution_matrix(rows, cfg)
    write_json(
        Path(cfg.get("model_dir", "results/swarmlineage/models")) / "training_summary.json",
        {
            "seeds": train_summaries,
            "metrics_path": str(metrics_path),
            "event_log_path": str(event_path),
            "order_log_path": str(order_path),
            "baseline_execution_matrix_rows": int(matrix.shape[0]),
        },
    )
    print(json.dumps({"metrics": str(metrics_path), "rows": int(metrics_frame.shape[0]), "events": int(len(all_events)), "order_rows": int(len(all_order_rows))}, indent=2))


if __name__ == "__main__":
    main()
