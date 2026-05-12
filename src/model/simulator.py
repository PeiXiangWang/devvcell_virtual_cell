from __future__ import annotations

from dataclasses import dataclass

import anndata as ad
import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import NearestNeighbors

from src.model.birth_death import growth_resample_indices
from src.model.cci import cci_context, cci_delta
from src.model.diffusion import adaptive_sigma
from src.model.pheromone import pheromone_delta
from src.model.swarm_rules import local_density, swarm_delta


class VelocityMLP(torch.nn.Module):
    def __init__(self, input_dim: int, hidden_dim: int, output_dim: int):
        super().__init__()
        self.net = torch.nn.Sequential(
            torch.nn.Linear(input_dim, hidden_dim),
            torch.nn.SiLU(),
            torch.nn.Linear(hidden_dim, hidden_dim),
            torch.nn.SiLU(),
            torch.nn.Linear(hidden_dim, output_dim),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


@dataclass(frozen=True)
class Variant:
    name: str
    intrinsic: bool
    birth_death: bool
    adaptive_diffusion: bool
    swarm: bool
    cci: bool
    pheromone: bool
    negative_control: str | None = None


VARIANTS = [
    Variant("M0_ot_interpolation", False, False, False, False, False, False),
    Variant("M1_intrinsic_neural", True, False, False, False, False, False),
    Variant("M2_ot_teacher_force", True, False, False, False, False, False),
    Variant("M3_ot_birth_death", True, True, False, False, False, False),
    Variant("M4_ot_adaptive_diffusion", True, False, True, False, False, False),
    Variant("M5_ot_swarm", True, False, False, True, False, False),
    Variant("M6_ot_swarm_birth_death", True, True, False, True, False, False),
    Variant("M7_ot_swarm_birth_death_diffusion", True, True, True, True, False, False),
    Variant("M8_ot_swarm_birth_death_diffusion_cci", True, True, True, True, True, False),
    Variant("M9_full_pheromone", True, True, True, True, True, True),
    Variant("M10_shuffled_time_ot", True, True, True, True, False, False, "time"),
    Variant("M11_random_lr_labels", True, True, True, True, True, False, "labels"),
]


def feature_matrix(adata: ad.AnnData, idx: np.ndarray, cfg: dict, context: dict | None = None) -> np.ndarray:
    z = np.asarray(adata.obsm[cfg.get("latent_key", "X_pca")][idx], dtype=np.float32)
    obs = adata.obs.iloc[idx]
    time = pd.to_numeric(obs[cfg.get("time_key", "time_numeric")], errors="coerce").to_numpy(dtype=float)
    all_time = pd.to_numeric(adata.obs[cfg.get("time_key", "time_numeric")], errors="coerce").to_numpy(dtype=float)
    time_norm = (time - np.nanmin(all_time)) / max(np.nanmax(all_time) - np.nanmin(all_time), 1e-8)
    entropy = pd.to_numeric(obs.get("ot_transition_entropy", 0.5), errors="coerce").fillna(0.5).to_numpy(dtype=float)
    growth = pd.to_numeric(obs.get("ot_growth", 1.0), errors="coerce").fillna(1.0).to_numpy(dtype=float)
    cycle = pd.to_numeric(obs.get("cell_cycle_score", 0.0), errors="coerce").fillna(0.0).to_numpy(dtype=float)
    density = local_density(z)
    extra = np.c_[time_norm, entropy, growth, cycle, density].astype(np.float32)
    return np.hstack([z, extra])


def train_velocity_model(adata: ad.AnnData, cfg: dict, seed: int) -> tuple[VelocityMLP, dict]:
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    holdout = float(cfg.get("holdout_time"))
    target_time = pd.to_numeric(adata.obs["ot_target_time"], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(target_time) & (target_time != holdout)
    idx = np.where(valid)[0]
    if idx.size > 12000:
        idx = rng.choice(idx, size=12000, replace=False)
    x = feature_matrix(adata, idx, cfg)
    y = np.asarray(adata.obsm[cfg.get("teacher_velocity_key", "X_ot_velocity")][idx], dtype=np.float32)
    model = VelocityMLP(x.shape[1], int(cfg.get("hidden_dim", 96)), y.shape[1])
    opt = torch.optim.AdamW(model.parameters(), lr=float(cfg.get("learning_rate", 1e-3)), weight_decay=1e-4)
    batch_size = int(cfg.get("batch_size", 512))
    xt = torch.as_tensor(x, dtype=torch.float32)
    yt = torch.as_tensor(y, dtype=torch.float32)
    losses = []
    for _ in range(int(cfg.get("epochs", 45))):
        order = torch.randperm(xt.shape[0])
        epoch_losses = []
        for start in range(0, xt.shape[0], batch_size):
            batch = order[start : start + batch_size]
            pred = model(xt[batch])
            loss = torch.mean((pred - yt[batch]) ** 2) + 1e-4 * torch.mean(pred**2)
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 2.0)
            opt.step()
            epoch_losses.append(float(loss.detach().cpu()))
        losses.append(float(np.mean(epoch_losses)))
    return model, {"train_cells": int(idx.size), "final_loss": losses[-1], "loss_curve": losses}


def _classify_labels(train_z: np.ndarray, train_labels: np.ndarray, pred_z: np.ndarray) -> np.ndarray:
    nn = NearestNeighbors(n_neighbors=min(7, train_z.shape[0])).fit(train_z)
    idx = nn.kneighbors(pred_z, return_distance=False)
    out = []
    for row in idx:
        vals, counts = np.unique(train_labels[row], return_counts=True)
        out.append(vals[np.argmax(counts)])
    return np.asarray(out, dtype=object)


def _linear_interpolation(adata: ad.AnnData, cfg: dict, source_idx: np.ndarray, next_idx: np.ndarray, alpha: float, seed: int) -> np.ndarray:
    z = np.asarray(adata.obsm[cfg.get("latent_key", "X_pca")], dtype=float)
    labels = adata.obs[cfg.get("cell_type_key", "lineage")].astype(str).to_numpy()
    rng = np.random.default_rng(seed)
    pred = np.empty_like(z[source_idx])
    for label in np.unique(labels[source_idx]):
        s = source_idx[labels[source_idx] == label]
        t = next_idx[labels[next_idx] == label]
        if t.size == 0:
            t = next_idx
        choose = rng.choice(t, size=s.size, replace=True)
        pred[labels[source_idx] == label] = (1 - alpha) * z[s] + alpha * z[choose]
    return pred


def simulate_variant(adata: ad.AnnData, cfg: dict, model: VelocityMLP, variant: Variant, seed: int) -> tuple[np.ndarray, np.ndarray, dict]:
    rng = np.random.default_rng(seed)
    time = pd.to_numeric(adata.obs[cfg.get("time_key", "time_numeric")], errors="coerce").to_numpy(dtype=float)
    holdout = float(cfg.get("holdout_time"))
    times = sorted(np.unique(time[np.isfinite(time)]))
    prev_time = max(t for t in times if t < holdout)
    next_time = min(t for t in times if t > holdout)
    alpha = (holdout - prev_time) / max(next_time - prev_time, 1e-8)
    source_idx_all = np.where(time == prev_time)[0]
    true_idx = np.where(time == holdout)[0]
    next_idx = np.where(time == next_time)[0]
    n_sim = min(int(cfg.get("simulation_cells_per_seed", 900)), source_idx_all.size)
    source_idx = rng.choice(source_idx_all, size=n_sim, replace=False)
    labels_all = adata.obs[cfg.get("cell_type_key", "lineage")].astype(str).to_numpy()
    z_all = np.asarray(adata.obsm[cfg.get("latent_key", "X_pca")], dtype=float)
    if variant.name == "M0_ot_interpolation":
        pred = _linear_interpolation(adata, cfg, source_idx, next_idx, alpha, seed)
    else:
        x = torch.as_tensor(feature_matrix(adata, source_idx, cfg), dtype=torch.float32)
        with torch.no_grad():
            delta = model(x).cpu().numpy()
        if variant.negative_control == "time":
            delta = rng.permutation(delta)
        pred = z_all[source_idx] + delta
    src_labels = labels_all[source_idx].copy()
    if variant.negative_control == "labels":
        src_labels = rng.permutation(src_labels)
    density = local_density(z_all[source_idx])
    entropy = pd.to_numeric(adata.obs.iloc[source_idx].get("ot_transition_entropy", 0.5), errors="coerce").fillna(0.5).to_numpy(dtype=float)
    growth = pd.to_numeric(adata.obs.iloc[source_idx].get("ot_growth", 1.0), errors="coerce").fillna(1.0).to_numpy(dtype=float)
    if variant.swarm:
        teacher_vel = np.asarray(adata.obsm[cfg.get("teacher_velocity_key", "X_ot_velocity")][source_idx], dtype=float)
        pred += swarm_delta(pred, teacher_vel, src_labels)
    if variant.cci:
        cci_scores_all, pairs = cci_context(adata, labels_all)
        pred += cci_delta(pred, src_labels, cci_scores_all[source_idx], strength=0.08)
    else:
        pairs = []
    if variant.pheromone:
        fate_cols = [c for c in adata.obs.columns if c.startswith("fate_prob_")]
        fate = adata.obs.iloc[source_idx][fate_cols].to_numpy(dtype=float) if fate_cols else np.empty((source_idx.size, 0))
        centers = []
        for col in fate_cols:
            label = col.replace("fate_prob_", "")
            mask = labels_all == label
            if mask.any():
                centers.append(z_all[mask].mean(axis=0))
        if centers:
            pred += pheromone_delta(pred, fate, np.vstack(centers))
    if variant.adaptive_diffusion:
        sigma = adaptive_sigma(entropy, density, growth)
        pred += rng.normal(scale=sigma[:, None], size=pred.shape)
    else:
        pred += rng.normal(scale=0.015, size=pred.shape)
    if variant.birth_death:
        true_n = min(int(cfg.get("simulation_cells_per_seed", 900)), true_idx.size)
        resample = growth_resample_indices(pd.Series(src_labels), true_n, growth, seed)
        pred = pred[resample]
        src_labels = src_labels[resample]
    true_z = z_all[true_idx]
    true_labels = labels_all[true_idx]
    pred_labels = _classify_labels(z_all, labels_all, pred)
    meta = {
        "source_time": float(prev_time),
        "holdout_time": float(holdout),
        "next_time": float(next_time),
        "n_pred": int(pred.shape[0]),
        "n_true": int(true_z.shape[0]),
        "cci_pairs_used": [f"{a}-{b}" for a, b in pairs],
        "pred_labels": pred_labels,
        "true_labels": true_labels,
    }
    return pred, true_z, meta
