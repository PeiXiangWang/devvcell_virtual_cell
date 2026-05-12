from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch

from src.model.dynamics import DynamicsFlags, SwarmLineageDynamics
from src.model.simulator import feature_matrix
from src.model.swarm_rules import local_density
from src.utils.config import ensure_dir, load_config


def configure_paths(config_path: str, quick_fixture: bool = False) -> tuple[dict, dict]:
    train_cfg = load_config(config_path)
    model_cfg = load_config(train_cfg.get("model_config", "configs/model.yaml"))
    if quick_fixture:
        model_cfg = dict(model_cfg)
        train_cfg = dict(train_cfg)
        model_cfg["teacher_path"] = "processed/quick_fixture/ot_teacher.h5ad"
        model_cfg["model_dir"] = "results/quick_fixture/models"
        model_cfg["event_log_path"] = "tables/quick_fixture/birth_death_event_log.csv"
        model_cfg["metrics_path"] = "tables/quick_fixture/final_metrics.csv"
        train_cfg["discovery_prefix"] = "quick_fixture"
    return train_cfg, model_cfg


def output_dirs(train_cfg: dict) -> tuple[Path, Path, Path]:
    prefix = train_cfg.get("discovery_prefix")
    table_dir = Path("tables") / prefix if prefix else Path("tables")
    report_dir = Path("reports") / prefix if prefix else Path("reports")
    fig_dir = Path("figures") / ("quick_fixture/discovery" if prefix else "discovery")
    ensure_dir(table_dir)
    ensure_dir(report_dir)
    ensure_dir(fig_dir)
    return table_dir, report_dir, fig_dir


def load_teacher(model_cfg: dict) -> ad.AnnData:
    return ad.read_h5ad(model_cfg["teacher_path"])


def load_teacher_model(adata: ad.AnnData, model_cfg: dict, seed: int = 7) -> SwarmLineageDynamics | None:
    model_path = Path(model_cfg.get("model_dir", "results/swarmlineage/models")) / f"dynamics_seed{seed}.pt"
    if not model_path.exists():
        return None
    x = feature_matrix(adata, np.arange(min(adata.n_obs, 8)), model_cfg)
    latent_dim = adata.obsm[model_cfg.get("latent_key", "X_pca")].shape[1]
    model = SwarmLineageDynamics(x.shape[1], latent_dim, int(model_cfg.get("hidden_dim", 96)))
    try:
        state = torch.load(model_path, map_location="cpu", weights_only=True)
    except TypeError:
        state = torch.load(model_path, map_location="cpu")
    model.load_state_dict(state["teacher_state_dict"])
    model.eval()
    return model


def cell_feature_frame(adata: ad.AnnData, model_cfg: dict, seed: int = 7) -> pd.DataFrame:
    idx = np.arange(adata.n_obs)
    z = np.asarray(adata.obsm[model_cfg.get("latent_key", "X_pca")], dtype=float)
    obs = adata.obs.copy()
    entropy = pd.to_numeric(obs.get("ot_transition_entropy", 0.5), errors="coerce").fillna(0.5).to_numpy(dtype=float)
    growth = pd.to_numeric(obs.get("ot_growth", 1.0), errors="coerce").fillna(1.0).to_numpy(dtype=float)
    cycle = pd.to_numeric(obs.get("cell_cycle_score", 0.0), errors="coerce").fillna(0.0).to_numpy(dtype=float)
    cci = pd.to_numeric(obs.get("cci_signal", 0.0), errors="coerce").fillna(0.0).to_numpy(dtype=float)
    fate_cols = [c for c in obs.columns if c.startswith("fate_prob_")]
    fate = obs[fate_cols].to_numpy(dtype=float) if fate_cols else np.zeros((adata.n_obs, 1))
    fate_max = fate.max(axis=1) if fate.size else np.zeros(adata.n_obs)
    fate_entropy = -np.sum(np.clip(fate, 1e-12, 1.0) * np.log(np.clip(fate, 1e-12, 1.0)), axis=1) / max(np.log(fate.shape[1]), 1.0)
    density = local_density(z)
    velocity = np.asarray(adata.obsm.get("X_ot_velocity", np.zeros_like(z)), dtype=float)
    displacement = np.linalg.norm(velocity, axis=1)
    frame = pd.DataFrame(
        {
            "cell_id": obs.index.astype(str),
            "time_numeric": pd.to_numeric(obs.get(model_cfg.get("time_key", "time_numeric")), errors="coerce").to_numpy(dtype=float),
            "lineage": obs.get(model_cfg.get("cell_type_key", "lineage"), "unknown").astype(str).to_numpy(),
            "ot_transition_entropy": entropy,
            "local_density": density,
            "fate_probability_max": fate_max,
            "fate_entropy": fate_entropy,
            "cell_cycle_score": cycle,
            "cci_signal": cci,
            "ot_growth": growth,
            "ot_displacement": displacement,
        }
    )
    model = load_teacher_model(adata, model_cfg, seed=seed)
    if model is not None:
        x = torch.as_tensor(feature_matrix(adata, idx, model_cfg), dtype=torch.float32)
        flags = DynamicsFlags(use_teacher=True, use_birth_death=True, use_diffusion=True, use_cci=True, use_swarm=True, use_memory=True)
        with torch.no_grad():
            frame["learned_sigma"] = model.sigma(x, flags).cpu().numpy()
            frame["birth_hazard"] = model.birth_hazard(x, flags).cpu().numpy()
            frame["death_hazard"] = model.death_hazard(x, flags).cpu().numpy()
            frame["net_growth_hazard"] = frame["birth_hazard"] - frame["death_hazard"]
            frame["learned_displacement"] = np.linalg.norm(model.vector_field(x, flags=flags).cpu().numpy(), axis=1)
    else:
        frame["learned_sigma"] = 0.015 + 0.12 * frame["ot_transition_entropy"]
        frame["birth_hazard"] = np.maximum(np.log(np.maximum(frame["ot_growth"], 1e-3)), 0.0)
        frame["death_hazard"] = np.maximum(-np.log(np.maximum(frame["ot_growth"], 1e-3)), 0.0)
        frame["net_growth_hazard"] = frame["birth_hazard"] - frame["death_hazard"]
        frame["learned_displacement"] = frame["ot_displacement"]
    return frame


def linear_effects(frame: pd.DataFrame, response: str, predictors: list[str]) -> pd.DataFrame:
    rows = []
    data = frame[[response] + predictors].replace([np.inf, -np.inf], np.nan).dropna()
    if data.shape[0] < len(predictors) + 5:
        return pd.DataFrame(columns=["response", "predictor", "coef", "abs_coef", "r2", "n"])
    y = data[response].to_numpy(dtype=float)
    x = data[predictors].to_numpy(dtype=float)
    x = (x - x.mean(axis=0)) / np.maximum(x.std(axis=0), 1e-8)
    y_center = y - y.mean()
    beta, *_ = np.linalg.lstsq(np.c_[np.ones(x.shape[0]), x], y_center, rcond=None)
    pred = np.c_[np.ones(x.shape[0]), x] @ beta + y.mean()
    r2 = 1.0 - float(np.sum((y - pred) ** 2) / max(np.sum((y - y.mean()) ** 2), 1e-12))
    for p, c in zip(predictors, beta[1:]):
        rows.append({"response": response, "predictor": p, "coef": float(c), "abs_coef": float(abs(c)), "r2": r2, "n": int(data.shape[0])})
    return pd.DataFrame(rows)


def write_report(path: Path, title: str, body: list[str]) -> None:
    from src.utils.config import write_text

    write_text(path, "\n".join([f"# {title}", "", *body, ""]))
