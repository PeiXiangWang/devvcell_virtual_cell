from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch
from sklearn.neighbors import NearestNeighbors

from src.model.birth_death import calibrate_count_hazard, stochastic_birth_death
from src.model.cci import cci_branch_delta, sender_receiver_graph
from src.model.dynamics import DynamicsFlags, SwarmLineageDynamics
from src.model.pheromone import MemoryField
from src.model.swarm_rules import local_density, swarm_delta
from src.ot_teacher.run_moscot import _sinkhorn_numpy, _scaled_cost


@dataclass(frozen=True)
class Variant:
    name: str
    baseline: str | None
    flags: DynamicsFlags
    model_role: str
    negative_control: str | None = None


VARIANTS = [
    Variant("M0_linear_label_interpolation", "linear", DynamicsFlags(), "none"),
    Variant("M0b_ot_interpolation", "ot", DynamicsFlags(), "none"),
    Variant("M1_intrinsic_neural", None, DynamicsFlags(), "intrinsic"),
    Variant("M2_ot_teacher_force", None, DynamicsFlags(use_teacher=True), "teacher"),
    Variant("M3_ot_birth_death", None, DynamicsFlags(use_teacher=True, use_birth_death=True), "teacher"),
    Variant("M4_ot_adaptive_diffusion", None, DynamicsFlags(use_teacher=True, use_diffusion=True), "teacher"),
    Variant("M5_ot_swarm", None, DynamicsFlags(use_teacher=True, use_swarm=True), "teacher"),
    Variant("M6_ot_swarm_birth_death", None, DynamicsFlags(use_teacher=True, use_swarm=True, use_birth_death=True), "teacher"),
    Variant("M7_ot_swarm_birth_death_diffusion", None, DynamicsFlags(use_teacher=True, use_swarm=True, use_birth_death=True, use_diffusion=True), "teacher"),
    Variant("M8_ot_swarm_birth_death_diffusion_cci", None, DynamicsFlags(use_teacher=True, use_swarm=True, use_birth_death=True, use_diffusion=True, use_cci=True), "teacher"),
    Variant("M9_full_memory", None, DynamicsFlags(use_teacher=True, use_swarm=True, use_birth_death=True, use_diffusion=True, use_cci=True, use_memory=True), "teacher"),
    Variant("M10_shuffled_time_ot", None, DynamicsFlags(use_teacher=True, use_swarm=True, use_birth_death=True, use_diffusion=True), "teacher", "time"),
    Variant("M11_random_lr_labels", None, DynamicsFlags(use_teacher=True, use_swarm=True, use_birth_death=True, use_diffusion=True, use_cci=True), "teacher", "labels"),
]


def _time_values(adata: ad.AnnData, cfg: dict) -> np.ndarray:
    return pd.to_numeric(adata.obs[cfg.get("time_key", "time_numeric")], errors="coerce").to_numpy(dtype=float)


def _fate_matrix(adata: ad.AnnData, idx: np.ndarray) -> np.ndarray:
    cols = [c for c in adata.obs.columns if c.startswith("fate_prob_")]
    if not cols:
        return np.empty((idx.size, 0), dtype=float)
    return adata.obs.iloc[idx][cols].to_numpy(dtype=float)


def _fate_entropy(fate: np.ndarray) -> np.ndarray:
    if fate.size == 0:
        return np.zeros(fate.shape[0], dtype=float)
    p = np.clip(fate, 1e-12, 1.0)
    return -np.sum(p * np.log(p), axis=1) / max(np.log(p.shape[1]), 1.0)


def feature_matrix_from_arrays(
    z: np.ndarray,
    time_norm: np.ndarray,
    entropy: np.ndarray,
    growth: np.ndarray,
    cycle: np.ndarray,
    density: np.ndarray,
    cci_signal: np.ndarray,
    fate_entropy: np.ndarray,
) -> np.ndarray:
    return np.hstack(
        [
            z.astype(np.float32),
            np.c_[
                time_norm,
                np.nan_to_num(entropy, nan=0.5),
                np.nan_to_num(growth, nan=1.0),
                np.nan_to_num(cycle, nan=0.0),
                np.nan_to_num(density, nan=0.5),
                np.nan_to_num(cci_signal, nan=0.0),
                np.nan_to_num(fate_entropy, nan=0.0),
            ].astype(np.float32),
        ]
    )


def feature_matrix(adata: ad.AnnData, idx: np.ndarray, cfg: dict) -> np.ndarray:
    z = np.asarray(adata.obsm[cfg.get("latent_key", "X_pca")][idx], dtype=np.float32)
    obs = adata.obs.iloc[idx]
    all_time = _time_values(adata, cfg)
    time = pd.to_numeric(obs[cfg.get("time_key", "time_numeric")], errors="coerce").to_numpy(dtype=float)
    time_norm = (time - np.nanmin(all_time)) / max(np.nanmax(all_time) - np.nanmin(all_time), 1e-8)
    entropy = pd.to_numeric(obs.get("ot_transition_entropy", 0.5), errors="coerce").fillna(0.5).to_numpy(dtype=float)
    growth = pd.to_numeric(obs.get("ot_growth", 1.0), errors="coerce").fillna(1.0).to_numpy(dtype=float)
    cycle = pd.to_numeric(obs.get("cell_cycle_score", 0.0), errors="coerce").fillna(0.0).to_numpy(dtype=float)
    density = local_density(z)
    cci = pd.to_numeric(obs.get("cci_signal", 0.0), errors="coerce").fillna(0.0).to_numpy(dtype=float)
    fate = _fate_matrix(adata, idx)
    return feature_matrix_from_arrays(z, time_norm, entropy, growth, cycle, density, cci, _fate_entropy(fate))


def _intrinsic_targets(adata: ad.AnnData, cfg: dict, idx: np.ndarray) -> np.ndarray:
    z = np.asarray(adata.obsm[cfg.get("latent_key", "X_pca")], dtype=float)
    time = _time_values(adata, cfg)
    labels = adata.obs[cfg.get("cell_type_key", "lineage")].astype(str).to_numpy()
    times = sorted(np.unique(time[np.isfinite(time)]))
    target = np.zeros((idx.size, z.shape[1]), dtype=np.float32)
    global_centers = {t: z[time == t].mean(axis=0) for t in times}
    for row, cell_idx in enumerate(idx):
        t = time[cell_idx]
        future = [x for x in times if x > t]
        if not future:
            continue
        nt = future[0]
        lab = labels[cell_idx]
        mask = (time == nt) & (labels == lab)
        center = z[mask].mean(axis=0) if mask.any() else global_centers[nt]
        target[row] = center - z[cell_idx]
    return target


def train_dynamics_model(adata: ad.AnnData, cfg: dict, seed: int, role: str) -> tuple[SwarmLineageDynamics, dict]:
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    holdout = float(cfg.get("holdout_time"))
    time = _time_values(adata, cfg)
    split = adata.obs.get("split_role", pd.Series("train", index=adata.obs_names)).astype(str).to_numpy()
    valid = (split != "eval_holdout") & np.isfinite(time) & (~np.isclose(time, holdout))
    if role == "teacher":
        target_time = pd.to_numeric(adata.obs.get("ot_target_time", np.nan), errors="coerce").to_numpy(dtype=float)
        valid &= np.isfinite(target_time)
    idx = np.where(valid)[0]
    if idx.size > 12000:
        idx = rng.choice(idx, size=12000, replace=False)
    x = feature_matrix(adata, idx, cfg)
    if role == "teacher":
        y = np.asarray(adata.obsm[cfg.get("teacher_velocity_key", "X_ot_velocity")][idx], dtype=np.float32)
    else:
        y = _intrinsic_targets(adata, cfg, idx)
    model = SwarmLineageDynamics(x.shape[1], y.shape[1], int(cfg.get("hidden_dim", 96)))
    opt = torch.optim.AdamW(model.parameters(), lr=float(cfg.get("learning_rate", 1e-3)), weight_decay=1e-4)
    batch_size = int(cfg.get("batch_size", 512))
    epochs = int(cfg.get("epochs", 45))
    xt = torch.as_tensor(x, dtype=torch.float32)
    yt = torch.as_tensor(y, dtype=torch.float32)
    entropy = torch.as_tensor(pd.to_numeric(adata.obs.iloc[idx].get("ot_transition_entropy", 0.5), errors="coerce").fillna(0.5).to_numpy(dtype=float), dtype=torch.float32)
    growth = torch.as_tensor(pd.to_numeric(adata.obs.iloc[idx].get("ot_growth", 1.0), errors="coerce").fillna(1.0).to_numpy(dtype=float), dtype=torch.float32)
    losses = []
    flags = DynamicsFlags(use_teacher=(role == "teacher"), use_birth_death=True, use_diffusion=True)
    for _ in range(epochs):
        order = torch.randperm(xt.shape[0])
        epoch_losses = []
        for start in range(0, xt.shape[0], batch_size):
            batch = order[start : start + batch_size]
            pred = model.vector_field(xt[batch], flags=flags)
            bary = torch.mean((pred - yt[batch]) ** 2)
            birth = model.birth_hazard(xt[batch], flags)
            death = model.death_hazard(xt[batch], flags)
            sigma = model.sigma(xt[batch], flags)
            growth_loss = torch.mean((birth - death - torch.log(torch.clamp(growth[batch], min=1e-3))) ** 2)
            diffusion_loss = torch.mean((sigma - (0.015 + 0.12 * entropy[batch])) ** 2)
            manifold = 1e-4 * torch.mean(pred**2)
            loss = bary + 0.2 * growth_loss + 0.2 * diffusion_loss + manifold
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 2.0)
            opt.step()
            epoch_losses.append(float(loss.detach().cpu()))
        losses.append(float(np.mean(epoch_losses)))
    return model, {"role": role, "train_cells": int(idx.size), "final_loss": losses[-1], "loss_curve": losses}


def _classify_labels(train_z: np.ndarray, train_labels: np.ndarray, pred_z: np.ndarray) -> np.ndarray:
    nn = NearestNeighbors(n_neighbors=min(7, train_z.shape[0])).fit(train_z)
    idx = nn.kneighbors(pred_z, return_distance=False)
    out = []
    for row in idx:
        vals, counts = np.unique(train_labels[row], return_counts=True)
        out.append(vals[np.argmax(counts)])
    return np.asarray(out, dtype=object)


def _split_times(adata: ad.AnnData, cfg: dict) -> tuple[float, float, float]:
    time = _time_values(adata, cfg)
    holdout = float(cfg.get("holdout_time"))
    times = sorted(np.unique(time[np.isfinite(time)]))
    return max(t for t in times if t < holdout), holdout, min(t for t in times if t > holdout)


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
        pred[labels[source_idx] == label] = (1.0 - alpha) * z[s] + alpha * z[choose]
    return pred


def _ot_interpolation(adata: ad.AnnData, cfg: dict, source_idx: np.ndarray, next_idx: np.ndarray, alpha: float) -> np.ndarray:
    z = np.asarray(adata.obsm[cfg.get("latent_key", "X_pca")], dtype=float)
    x, y = z[source_idx], z[next_idx]
    a = np.full(x.shape[0], 1.0 / x.shape[0])
    b = np.full(y.shape[0], 1.0 / y.shape[0])
    plan = _sinkhorn_numpy(a, b, _scaled_cost(x, y), reg=float(cfg.get("epsilon", 0.08)))
    transition = plan / np.maximum(plan.sum(axis=1, keepdims=True), 1e-12)
    bary = transition @ y
    return (1.0 - alpha) * x + alpha * bary


def _agent_state(adata: ad.AnnData, idx: np.ndarray, cfg: dict) -> dict:
    obs = adata.obs.iloc[idx]
    return {
        "entropy": pd.to_numeric(obs.get("ot_transition_entropy", 0.5), errors="coerce").fillna(0.5).to_numpy(dtype=float),
        "growth": pd.to_numeric(obs.get("ot_growth", 1.0), errors="coerce").fillna(1.0).to_numpy(dtype=float),
        "cycle": pd.to_numeric(obs.get("cell_cycle_score", 0.0), errors="coerce").fillna(0.0).to_numpy(dtype=float),
        "cci": pd.to_numeric(obs.get("cci_signal", 0.0), errors="coerce").fillna(0.0).to_numpy(dtype=float),
        "fate": _fate_matrix(adata, idx),
    }


def _rollout(
    adata: ad.AnnData,
    cfg: dict,
    model: SwarmLineageDynamics,
    variant: Variant,
    source_idx: np.ndarray,
    seed: int,
    prev_time: float,
    holdout: float,
    next_time: float,
) -> tuple[np.ndarray, np.ndarray, list[dict]]:
    rng = np.random.default_rng(seed)
    z_all = np.asarray(adata.obsm[cfg.get("latent_key", "X_pca")], dtype=float)
    labels_all = adata.obs[cfg.get("cell_type_key", "lineage")].astype(str).to_numpy()
    z = z_all[source_idx].copy()
    labels = labels_all[source_idx].copy()
    state = _agent_state(adata, source_idx, cfg)
    if variant.negative_control == "labels":
        labels = rng.permutation(labels)
    all_z = z_all
    all_labels = labels_all
    cci_graph, cci_all, lr_pairs = sender_receiver_graph(adata, all_z, all_labels)
    state["cci"] = cci_all[source_idx] if cci_all.size == adata.n_obs else state["cci"]
    memory = MemoryField(max(state["fate"].shape[1], 1))
    event_rows: list[dict] = []
    steps = int(cfg.get("rollout_steps", 4))
    dt = float(cfg.get("dt", 0.25))
    for step in range(steps):
        frac = (step + 1) / max(steps, 1)
        current_time = prev_time + frac * (holdout - prev_time)
        density = local_density(z)
        time_norm = np.full(z.shape[0], (current_time - prev_time) / max(next_time - prev_time, 1e-8))
        fate_entropy = _fate_entropy(state["fate"])
        x = feature_matrix_from_arrays(z, time_norm, state["entropy"], state["growth"], state["cycle"], density, state["cci"], fate_entropy)
        sw = swarm_delta(z, np.zeros_like(z), labels) if variant.flags.use_swarm else np.zeros_like(z)
        cci_delta = cci_branch_delta(z, labels, state["cci"]) if variant.flags.use_cci else np.zeros_like(z)
        if variant.flags.use_memory:
            memory.step(z, state["fate"], dt)
            mem = memory.gradient_delta(z, state["fate"])
        else:
            mem = np.zeros_like(z)
        with torch.no_grad():
            xt = torch.as_tensor(x, dtype=torch.float32)
            v = model.vector_field(
                xt,
                swarm_delta=torch.as_tensor(sw, dtype=torch.float32),
                cci_delta=torch.as_tensor(cci_delta, dtype=torch.float32),
                memory_delta=torch.as_tensor(mem, dtype=torch.float32),
                flags=variant.flags,
            ).cpu().numpy()
            birth = model.birth_hazard(xt, variant.flags).cpu().numpy()
            death = model.death_hazard(xt, variant.flags).cpu().numpy()
            sigma = model.sigma(xt, variant.flags).cpu().numpy()
        if variant.negative_control == "time":
            v = rng.permutation(v)
        z = z + dt * v
        if variant.flags.use_diffusion:
            z = z + rng.normal(scale=np.maximum(sigma, 1e-5)[:, None] * np.sqrt(dt), size=z.shape)
        if variant.flags.use_birth_death:
            bd = stochastic_birth_death(birth, death, dt, seed + 101 * step, current_time, labels)
            daughter_idx = bd.daughter_parent_indices
            keep = bd.keep_mask
            daughters = z[daughter_idx] + rng.normal(scale=0.02, size=(daughter_idx.size, z.shape[1]))
            daughter_labels = labels[daughter_idx].copy()
            z = np.vstack([z[keep], daughters]) if daughters.size else z[keep]
            labels = np.r_[labels[keep], daughter_labels] if daughters.size else labels[keep]
            for key in ("entropy", "growth", "cycle", "cci"):
                vals = state[key]
                state[key] = np.r_[vals[keep], vals[daughter_idx]] if daughter_idx.size else vals[keep]
            fate = state["fate"]
            if fate.size:
                daughter_fate = fate[daughter_idx].copy()
                if daughter_fate.size:
                    noise = rng.normal(scale=0.03, size=daughter_fate.shape)
                    daughter_fate = np.clip(daughter_fate + noise, 1e-6, None)
                    daughter_fate /= daughter_fate.sum(axis=1, keepdims=True)
                state["fate"] = np.vstack([fate[keep], daughter_fate]) if daughter_idx.size else fate[keep]
            event_rows.extend([{**row, "variant": variant.name, "step": step} for row in bd.event_rows])
            if z.shape[0] == 0:
                z = z_all[source_idx[:1]].copy()
                labels = labels_all[source_idx[:1]].copy()
                state = _agent_state(adata, source_idx[:1], cfg)
    return z, labels, event_rows


def simulate_variant(adata: ad.AnnData, cfg: dict, models: dict[str, SwarmLineageDynamics], variant: Variant, seed: int) -> tuple[np.ndarray, np.ndarray, dict]:
    rng = np.random.default_rng(seed)
    time = _time_values(adata, cfg)
    prev_time, holdout, next_time = _split_times(adata, cfg)
    alpha = (holdout - prev_time) / max(next_time - prev_time, 1e-8)
    source_idx_all = np.where(time == prev_time)[0]
    true_idx = np.where(time == holdout)[0]
    next_idx = np.where(time == next_time)[0]
    n_sim = min(int(cfg.get("simulation_cells_per_seed", 900)), source_idx_all.size)
    source_idx = rng.choice(source_idx_all, size=n_sim, replace=False)
    labels_all = adata.obs[cfg.get("cell_type_key", "lineage")].astype(str).to_numpy()
    if variant.baseline == "linear":
        pred = _linear_interpolation(adata, cfg, source_idx, next_idx, alpha, seed)
        pred_labels = labels_all[source_idx]
        event_rows: list[dict] = []
    elif variant.baseline == "ot":
        pred = _ot_interpolation(adata, cfg, source_idx, next_idx, alpha)
        pred_labels = labels_all[source_idx]
        event_rows = []
    else:
        model = models[variant.model_role]
        pred, pred_labels, event_rows = _rollout(adata, cfg, model, variant, source_idx, seed, prev_time, holdout, next_time)
    z_all = np.asarray(adata.obsm[cfg.get("latent_key", "X_pca")], dtype=float)
    true_z = z_all[true_idx]
    true_labels = labels_all[true_idx]
    pred_labels = _classify_labels(z_all, labels_all, pred) if pred.shape[0] else np.array([], dtype=object)
    meta = {
        "source_time": float(prev_time),
        "holdout_time": float(holdout),
        "next_time": float(next_time),
        "n_pred": int(pred.shape[0]),
        "n_true": int(true_z.shape[0]),
        "pred_labels": pred_labels,
        "true_labels": true_labels,
        "event_rows": event_rows,
    }
    return pred, true_z, meta

