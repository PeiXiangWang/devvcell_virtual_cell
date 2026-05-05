"""Train and evaluate cell-level developmental transition baselines.

This is the first cell-level model layer for DevVCell. It uses cross-sectional
single-cell states, constructs adjacent-stage pseudo-targets within each broad
system using either nearest neighbors or entropic OT/Sinkhorn barycenters, and
compares identity, mean-shift, ridge and optional MLP transition operators in a
low-dimensional expression space.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import joblib
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from devvcell.io import load_json, resolve_project_path, write_json  # noqa: E402


try:
    import torch
    import torch.nn as nn
    from torch.utils.data import DataLoader, TensorDataset
except ImportError:  # pragma: no cover - torch is optional.
    torch = None
    nn = None
    DataLoader = None
    TensorDataset = None


@dataclass(frozen=True)
class PairSet:
    x: np.ndarray
    y: np.ndarray
    system: np.ndarray
    src_stage: np.ndarray
    tgt_stage: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        default="config/cell_level_baseline.json",
        help="Path to the cell-level experiment config.",
    )
    return parser.parse_args()


def as_csr_float32(matrix) -> sparse.csr_matrix:
    if sparse.issparse(matrix):
        return matrix.tocsr().astype(np.float32)
    return sparse.csr_matrix(np.asarray(matrix, dtype=np.float32))


def fit_latent(X: sparse.csr_matrix, train_cell_mask: np.ndarray, latent_dim: int, seed: int) -> tuple[np.ndarray, TruncatedSVD, StandardScaler]:
    n_components = min(latent_dim, X.shape[1] - 1, max(2, int(train_cell_mask.sum()) - 1))
    svd = TruncatedSVD(n_components=n_components, random_state=seed)
    train_latent = svd.fit_transform(X[train_cell_mask])
    scaler = StandardScaler()
    scaler.fit(train_latent)
    latent = scaler.transform(svd.transform(X)).astype(np.float32)
    return latent, svd, scaler


def squared_euclidean_cost(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    x_norm = np.sum(x * x, axis=1, keepdims=True)
    y_norm = np.sum(y * y, axis=1, keepdims=True).T
    cost = x_norm + y_norm - 2.0 * np.matmul(x, y.T)
    return np.maximum(cost, 0.0).astype(np.float32)


def sinkhorn_barycentric_targets(
    source: np.ndarray,
    target_pool: np.ndarray,
    epsilon: float,
    iterations: int,
) -> tuple[np.ndarray, dict[str, float]]:
    """Return OT barycentric target states for each source cell.

    The cost matrix is divided by its positive median before Sinkhorn scaling,
    so the same epsilon value is usable across systems and stage transitions.
    """

    cost = squared_euclidean_cost(source, target_pool)
    positive = cost[cost > 0]
    median_cost = float(np.median(positive)) if len(positive) else 1.0
    if not np.isfinite(median_cost) or median_cost <= 0:
        median_cost = 1.0
    cost = cost / median_cost

    eps = max(float(epsilon), 1e-4)
    kernel = np.exp(-cost / eps).astype(np.float64)
    kernel = np.maximum(kernel, 1e-300)

    n_src, n_tgt = kernel.shape
    a = np.full(n_src, 1.0 / n_src, dtype=np.float64)
    b = np.full(n_tgt, 1.0 / n_tgt, dtype=np.float64)
    u = np.ones(n_src, dtype=np.float64)
    v = np.ones(n_tgt, dtype=np.float64)
    tiny = 1e-300

    for _ in range(int(iterations)):
        u = a / np.maximum(kernel @ v, tiny)
        v = b / np.maximum(kernel.T @ u, tiny)

    plan = (u[:, None] * kernel) * v[None, :]
    row_mass = np.maximum(plan.sum(axis=1, keepdims=True), tiny)
    barycentric = (plan @ target_pool.astype(np.float64)) / row_mass
    diagnostics = {
        "cost_median": median_cost,
        "transport_entropy": float(-np.sum(plan * np.log(np.maximum(plan, tiny)))),
        "row_mass_abs_error": float(np.mean(np.abs(plan.sum(axis=1) - a))),
        "col_mass_abs_error": float(np.mean(np.abs(plan.sum(axis=0) - b))),
    }
    return barycentric.astype(np.float32), diagnostics


def build_transition_targets(
    source: np.ndarray,
    target_pool: np.ndarray,
    config: dict,
) -> tuple[np.ndarray, dict[str, object]]:
    pairing_cfg = config["model"].get("pairing", {})
    method = str(pairing_cfg.get("method", "nearest")).lower()
    if method in {"nearest", "nearest_neighbor", "nn"}:
        nn_model = NearestNeighbors(n_neighbors=1, metric="euclidean")
        nn_model.fit(target_pool)
        nn_pos = nn_model.kneighbors(source, return_distance=False).ravel()
        return target_pool[nn_pos].astype(np.float32), {"pairing_method": "nearest"}

    if method in {"sinkhorn", "ot", "entropic_ot"}:
        target, diagnostics = sinkhorn_barycentric_targets(
            source,
            target_pool,
            epsilon=float(pairing_cfg.get("sinkhorn_epsilon", 0.12)),
            iterations=int(pairing_cfg.get("sinkhorn_iterations", 120)),
        )
        return target, {"pairing_method": "sinkhorn", **diagnostics}

    raise ValueError(f"Unknown pairing method: {method!r}")


def build_pair_set(
    latent: np.ndarray,
    obs: pd.DataFrame,
    config: dict,
    train: bool,
    rng: np.random.Generator,
) -> tuple[PairSet, pd.DataFrame]:
    heldout_targets = set(int(s) for s in config["heldout_target_stages"])
    max_pairs = int(config["model"]["max_pairs_per_transition"])
    pairing_cfg = config["model"].get("pairing", {})
    target_pool_multiplier = max(1, int(pairing_cfg.get("target_pool_multiplier", 2)))
    rows: list[dict[str, object]] = []
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    systems: list[str] = []
    src_stages: list[int] = []
    tgt_stages: list[int] = []

    for system_name in sorted(obs["devvcell_system"].astype(str).unique()):
        if system_name == "unassigned":
            continue
        system_obs = obs[obs["devvcell_system"].astype(str) == system_name]
        available_stages = sorted(system_obs["stage_num"].astype(int).unique())
        for src_stage in available_stages:
            tgt_stage = int(src_stage) + 1
            if tgt_stage not in available_stages:
                continue
            is_heldout = tgt_stage in heldout_targets
            if train and is_heldout:
                continue
            if (not train) and (not is_heldout):
                continue

            src_idx = system_obs.index[system_obs["stage_num"].astype(int) == src_stage].to_numpy()
            tgt_idx = system_obs.index[system_obs["stage_num"].astype(int) == tgt_stage].to_numpy()
            pair_count = min(max_pairs, len(src_idx), len(tgt_idx))
            if pair_count < 10:
                continue

            src_sample = rng.choice(src_idx, size=pair_count, replace=False)
            tgt_pool = rng.choice(
                tgt_idx,
                size=min(len(tgt_idx), max_pairs * target_pool_multiplier),
                replace=False,
            )
            source_states = latent[src_sample]
            target_states, pairing_info = build_transition_targets(source_states, latent[tgt_pool], config)

            xs.append(source_states)
            ys.append(target_states)
            systems.extend([system_name] * pair_count)
            src_stages.extend([int(src_stage)] * pair_count)
            tgt_stages.extend([int(tgt_stage)] * pair_count)
            rows.append(
                {
                    "split": "train" if train else "heldout",
                    "system": system_name,
                    "src_stage": int(src_stage),
                    "tgt_stage": int(tgt_stage),
                    "pairs": int(pair_count),
                    "src_cells_available": int(len(src_idx)),
                    "tgt_cells_available": int(len(tgt_idx)),
                    "target_pool_cells": int(len(tgt_pool)),
                    **pairing_info,
                }
            )

    if not xs:
        raise RuntimeError("No pseudo-pairs were built. Check stages and heldout_target_stages.")

    pair_set = PairSet(
        x=np.vstack(xs).astype(np.float32),
        y=np.vstack(ys).astype(np.float32),
        system=np.asarray(systems, dtype=object),
        src_stage=np.asarray(src_stages, dtype=np.int32),
        tgt_stage=np.asarray(tgt_stages, dtype=np.int32),
    )
    return pair_set, pd.DataFrame(rows)


def fit_mean_shift(train_pairs: PairSet) -> tuple[dict[str, np.ndarray], np.ndarray]:
    deltas = train_pairs.y - train_pairs.x
    global_delta = deltas.mean(axis=0)
    by_system: dict[str, np.ndarray] = {}
    for system_name in sorted(set(train_pairs.system)):
        mask = train_pairs.system == system_name
        by_system[str(system_name)] = deltas[mask].mean(axis=0)
    return by_system, global_delta


def predict_mean_shift(pair_set: PairSet, by_system: dict[str, np.ndarray], global_delta: np.ndarray) -> np.ndarray:
    pred = np.empty_like(pair_set.x)
    for system_name in sorted(set(pair_set.system)):
        mask = pair_set.system == system_name
        delta = by_system.get(str(system_name), global_delta)
        pred[mask] = pair_set.x[mask] + delta
    return pred


def save_model_artifacts(
    model_dir: Path,
    svd: TruncatedSVD,
    scaler: StandardScaler,
    ridge: Ridge,
    shift_by_system: dict[str, np.ndarray],
    global_shift: np.ndarray,
    config: dict,
) -> None:
    joblib.dump(svd, model_dir / "state_svd.joblib")
    joblib.dump(scaler, model_dir / "state_scaler.joblib")
    joblib.dump(ridge, model_dir / "transition_ridge.joblib")
    np.savez(
        model_dir / "mean_shift_vectors.npz",
        global_shift=global_shift.astype(np.float32),
        **{f"system__{name}": vec.astype(np.float32) for name, vec in shift_by_system.items()},
    )
    with (model_dir / "model_config_snapshot.json").open("w", encoding="utf-8") as handle:
        json.dump(config, handle, ensure_ascii=False, indent=2)


class TransitionMLP(nn.Module):
    def __init__(self, dim: int, hidden_dim: int) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(dim, hidden_dim),
            nn.GELU(),
            nn.LayerNorm(hidden_dim),
            nn.Linear(hidden_dim, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, dim),
        )

    def forward(self, x):  # noqa: D401
        return self.net(x)


class ContextResidualMLP(nn.Module):
    """Stage/system-conditioned residual transition operator."""

    def __init__(self, state_dim: int, context_dim: int, hidden_dim: int, dropout: float) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(state_dim + context_dim, hidden_dim),
            nn.GELU(),
            nn.LayerNorm(hidden_dim),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.GELU(),
            nn.LayerNorm(hidden_dim),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, state_dim),
        )

    def forward(self, x, context):  # noqa: D401
        residual = self.net(torch.cat([x, context], dim=1))
        return x + residual


def train_mlp(train_pairs: PairSet, config: dict, seed: int):
    if torch is None:
        raise RuntimeError("Torch is not installed, but MLP training is enabled.")

    mlp_cfg = config["model"]["mlp"]
    torch.manual_seed(seed)
    model = TransitionMLP(train_pairs.x.shape[1], int(mlp_cfg["hidden_dim"]))
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=float(mlp_cfg["learning_rate"]),
        weight_decay=float(mlp_cfg["weight_decay"]),
    )
    dataset = TensorDataset(torch.from_numpy(train_pairs.x), torch.from_numpy(train_pairs.y))
    loader = DataLoader(dataset, batch_size=int(mlp_cfg["batch_size"]), shuffle=True)
    loss_fn = nn.MSELoss()
    history: list[float] = []

    model.train()
    for _ in range(int(mlp_cfg["epochs"])):
        epoch_loss = 0.0
        n_seen = 0
        for xb, yb in loader:
            optimizer.zero_grad(set_to_none=True)
            loss = loss_fn(model(xb), yb)
            loss.backward()
            optimizer.step()
            batch_n = xb.shape[0]
            epoch_loss += float(loss.detach()) * batch_n
            n_seen += batch_n
        history.append(epoch_loss / max(1, n_seen))
    return model, history


def predict_mlp(model, x: np.ndarray) -> np.ndarray:
    model.eval()
    with torch.no_grad():
        return model(torch.from_numpy(x.astype(np.float32))).cpu().numpy().astype(np.float32)


def build_context_features(pair_set: PairSet, config: dict, system_order: list[str] | None = None) -> tuple[np.ndarray, list[str]]:
    if system_order is None:
        configured = [str(name) for name in config.get("systems", {}).keys()]
        observed = sorted(str(name) for name in set(pair_set.system))
        system_order = [name for name in configured if name in observed] + [name for name in observed if name not in configured]

    stage_values = [int(s) for s in config.get("stages", [])]
    stage_min = min(stage_values) if stage_values else int(np.min(pair_set.src_stage))
    stage_max = max(stage_values) if stage_values else int(np.max(pair_set.tgt_stage))
    stage_span = max(1, stage_max - stage_min)

    src_norm = ((pair_set.src_stage.astype(np.float32) - stage_min) / stage_span)[:, None]
    tgt_norm = ((pair_set.tgt_stage.astype(np.float32) - stage_min) / stage_span)[:, None]
    delta_norm = ((pair_set.tgt_stage.astype(np.float32) - pair_set.src_stage.astype(np.float32)) / stage_span)[:, None]

    onehot = np.zeros((len(pair_set.x), len(system_order)), dtype=np.float32)
    system_to_idx = {name: idx for idx, name in enumerate(system_order)}
    for idx, system_name in enumerate(pair_set.system.astype(str)):
        pos = system_to_idx.get(system_name)
        if pos is not None:
            onehot[idx, pos] = 1.0

    context = np.hstack([src_norm, tgt_norm, delta_norm, onehot]).astype(np.float32)
    return context, list(system_order)


def train_context_residual_mlp(train_pairs: PairSet, config: dict, seed: int):
    if torch is None:
        raise RuntimeError("Torch is not installed, but context residual MLP training is enabled.")

    model_cfg = config["model"]["context_residual_mlp"]
    context, system_order = build_context_features(train_pairs, config)
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed + 1729)
    model = ContextResidualMLP(
        state_dim=train_pairs.x.shape[1],
        context_dim=context.shape[1],
        hidden_dim=int(model_cfg["hidden_dim"]),
        dropout=float(model_cfg.get("dropout", 0.0)),
    )
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=float(model_cfg["learning_rate"]),
        weight_decay=float(model_cfg["weight_decay"]),
    )
    n_pairs = len(train_pairs.x)
    validation_fraction = float(model_cfg.get("validation_fraction", 0.0))
    validation_size = int(round(n_pairs * validation_fraction))
    validation_size = min(max(validation_size, 0), max(0, n_pairs - 1))
    order = rng.permutation(n_pairs)
    val_idx = order[:validation_size]
    train_idx = order[validation_size:]
    if len(train_idx) == 0:
        train_idx = order
        val_idx = np.array([], dtype=np.int64)

    dataset = TensorDataset(
        torch.from_numpy(train_pairs.x[train_idx]),
        torch.from_numpy(context[train_idx]),
        torch.from_numpy(train_pairs.y[train_idx]),
    )
    loader = DataLoader(dataset, batch_size=int(model_cfg["batch_size"]), shuffle=True)
    val_tensors = None
    if len(val_idx):
        val_tensors = (
            torch.from_numpy(train_pairs.x[val_idx]),
            torch.from_numpy(context[val_idx]),
            torch.from_numpy(train_pairs.y[val_idx]),
        )
    loss_fn = nn.MSELoss()
    history: list[dict[str, float | int | None]] = []
    best_state = deepcopy(model.state_dict())
    best_val_loss = float("inf")
    best_epoch = 0
    patience = int(model_cfg.get("patience", 0))
    min_delta = float(model_cfg.get("min_delta", 0.0))
    stale_epochs = 0

    for epoch in range(int(model_cfg["epochs"])):
        model.train()
        epoch_loss = 0.0
        n_seen = 0
        for xb, cb, yb in loader:
            optimizer.zero_grad(set_to_none=True)
            pred = model(xb, cb)
            loss = loss_fn(pred, yb)
            loss.backward()
            optimizer.step()
            batch_n = xb.shape[0]
            epoch_loss += float(loss.detach()) * batch_n
            n_seen += batch_n
        train_loss = epoch_loss / max(1, n_seen)
        val_loss = None
        if val_tensors is not None:
            model.eval()
            with torch.no_grad():
                xb, cb, yb = val_tensors
                val_loss = float(loss_fn(model(xb, cb), yb).detach())
            if val_loss + min_delta < best_val_loss:
                best_val_loss = val_loss
                best_state = deepcopy(model.state_dict())
                best_epoch = epoch + 1
                stale_epochs = 0
            else:
                stale_epochs += 1
        else:
            best_state = deepcopy(model.state_dict())
            best_epoch = epoch + 1
        history.append({"epoch": epoch + 1, "train_loss": train_loss, "validation_loss": val_loss})
        if patience > 0 and val_tensors is not None and stale_epochs >= patience:
            break

    model.load_state_dict(best_state)
    context_info = {
        "system_order": system_order,
        "context_columns": ["src_stage_norm", "tgt_stage_norm", "stage_delta_norm"]
        + [f"system__{name}" for name in system_order],
        "validation_fraction": validation_fraction,
        "n_train_internal_pairs": int(len(train_idx)),
        "n_validation_internal_pairs": int(len(val_idx)),
        "best_epoch": int(best_epoch),
        "best_validation_loss": best_val_loss if np.isfinite(best_val_loss) else None,
    }
    return model, context_info, history


def predict_context_residual_mlp(model, pair_set: PairSet, config: dict, system_order: list[str]) -> np.ndarray:
    context, _ = build_context_features(pair_set, config, system_order=system_order)
    model.eval()
    with torch.no_grad():
        return model(
            torch.from_numpy(pair_set.x.astype(np.float32)),
            torch.from_numpy(context.astype(np.float32)),
        ).cpu().numpy().astype(np.float32)


def rbf_mmd(x: np.ndarray, y: np.ndarray, max_cells: int, rng: np.random.Generator) -> float:
    if len(x) > max_cells:
        x = x[rng.choice(np.arange(len(x)), size=max_cells, replace=False)]
    if len(y) > max_cells:
        y = y[rng.choice(np.arange(len(y)), size=max_cells, replace=False)]
    if len(x) < 2 or len(y) < 2:
        return float("nan")

    pooled = np.vstack([x, y])
    diffs = pooled[:, None, :] - pooled[None, :, :]
    sq_dists = np.sum(diffs * diffs, axis=2)
    median_sq = np.median(sq_dists[sq_dists > 0])
    gamma = 1.0 / (2.0 * median_sq) if np.isfinite(median_sq) and median_sq > 0 else 1.0

    def kernel(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        d = a[:, None, :] - b[None, :, :]
        return np.exp(-gamma * np.sum(d * d, axis=2))

    kxx = kernel(x, x)
    kyy = kernel(y, y)
    kxy = kernel(x, y)
    return float(kxx.mean() + kyy.mean() - 2.0 * kxy.mean())


def evaluate_predictions(
    pair_set: PairSet,
    predictions: dict[str, np.ndarray],
    config: dict,
    rng: np.random.Generator,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    max_cells = int(config["model"]["mmd_max_cells"])
    for model_name, pred in predictions.items():
        for system_name in sorted(set(pair_set.system)):
            for tgt_stage in sorted(set(pair_set.tgt_stage[pair_set.system == system_name])):
                mask = (pair_set.system == system_name) & (pair_set.tgt_stage == tgt_stage)
                if mask.sum() == 0:
                    continue
                y_true = pair_set.y[mask]
                y_pred = pred[mask]
                centroid_mse = float(mean_squared_error(y_true.mean(axis=0), y_pred.mean(axis=0)))
                rows.append(
                    {
                        "model": model_name,
                        "split": "heldout",
                        "system": str(system_name),
                        "src_stage": int(pair_set.src_stage[mask][0]),
                        "tgt_stage": int(tgt_stage),
                        "n_pairs": int(mask.sum()),
                        "pair_latent_mse": float(mean_squared_error(y_true, y_pred)),
                        "centroid_latent_mse": centroid_mse,
                        "rbf_mmd": rbf_mmd(y_pred, y_true, max_cells=max_cells, rng=rng),
                    }
                )
    return pd.DataFrame(rows)


def plot_metrics(metrics: pd.DataFrame, figures_dir: Path) -> None:
    summary = (
        metrics.groupby("model", as_index=False)["pair_latent_mse"]
        .mean()
        .sort_values("pair_latent_mse", ascending=True)
    )
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.bar(summary["model"], summary["pair_latent_mse"], color=["#6c757d", "#4e79a7", "#59a14f", "#b07aa1"][: len(summary)])
    ax.set_ylabel("Heldout pseudo-pair latent MSE")
    ax.set_xlabel("Transition model")
    ax.set_title("Cell-level adjacent-stage transition baselines")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(figures_dir / "cell_level_transition_baseline_mse.png", dpi=220)
    plt.close(fig)


def train_cell_mask(obs: pd.DataFrame, heldout_targets: set[int]) -> np.ndarray:
    stages = obs["stage_num"].astype(int).to_numpy()
    return ~np.isin(stages, list(heldout_targets))


def main() -> None:
    args = parse_args()
    config = load_json(args.config)
    seed = int(config["seed"])
    rng = np.random.default_rng(seed)

    subset_path = resolve_project_path(config["processed_subset_path"])
    results_dir = resolve_project_path(config["results_dir"])
    tables_dir = results_dir / "tables"
    figures_dir = results_dir / "figures"
    model_dir = results_dir / "models"
    for path in [tables_dir, figures_dir, model_dir]:
        path.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(subset_path)
    X = as_csr_float32(adata.X)
    obs = adata.obs.reset_index(drop=True).copy()
    if "devvcell_system" not in obs.columns or "stage_num" not in obs.columns:
        raise ValueError("Subset is missing devvcell_system or stage_num. Run export_cell_level_subset.py first.")

    heldout_targets = set(int(s) for s in config["heldout_target_stages"])
    latent, svd, scaler = fit_latent(
        X,
        train_cell_mask(obs, heldout_targets),
        int(config["model"]["latent_dim"]),
        seed,
    )
    train_pairs, train_manifest = build_pair_set(latent, obs, config, train=True, rng=rng)
    heldout_pairs, heldout_manifest = build_pair_set(latent, obs, config, train=False, rng=rng)

    ridge = Ridge(alpha=float(config["model"]["ridge_alpha"]))
    ridge.fit(train_pairs.x, train_pairs.y)
    shift_by_system, global_shift = fit_mean_shift(train_pairs)
    save_model_artifacts(model_dir, svd, scaler, ridge, shift_by_system, global_shift, config)

    predictions = {
        "identity": heldout_pairs.x,
        "mean_shift": predict_mean_shift(heldout_pairs, shift_by_system, global_shift),
        "ridge": ridge.predict(heldout_pairs.x).astype(np.float32),
    }
    training_history = None
    if bool(config["model"]["mlp"]["enabled"]):
        mlp, training_history = train_mlp(train_pairs, config, seed)
        predictions["mlp"] = predict_mlp(mlp, heldout_pairs.x)
        torch.save(mlp.state_dict(), model_dir / "transition_mlp.pt")

    context_residual_history = None
    context_residual_info = None
    context_cfg = config["model"].get("context_residual_mlp", {})
    if bool(context_cfg.get("enabled", False)):
        context_model, context_residual_info, context_residual_history = train_context_residual_mlp(train_pairs, config, seed)
        predictions["context_residual_mlp"] = predict_context_residual_mlp(
            context_model,
            heldout_pairs,
            config,
            context_residual_info["system_order"],
        )
        torch.save(
            {
                "state_dict": context_model.state_dict(),
                "context_info": context_residual_info,
                "config": context_cfg,
            },
            model_dir / "transition_context_residual_mlp.pt",
        )

    metrics = evaluate_predictions(heldout_pairs, predictions, config, rng)
    pair_manifest = pd.concat([train_manifest, heldout_manifest], ignore_index=True)

    metrics.to_csv(tables_dir / "cell_level_transition_metrics.csv", index=False)
    pair_manifest.to_csv(tables_dir / "cell_level_pair_manifest.csv", index=False)
    plot_metrics(metrics, figures_dir)

    def rel(path: Path) -> str:
        try:
            return str(path.relative_to(PROJECT_ROOT)).replace("\\", "/")
        except ValueError:
            return str(path)

    summary = {
        "project": config["project"],
        "subset_path": str(subset_path),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "latent_dim": int(latent.shape[1]),
        "svd_explained_variance_ratio_sum": float(np.sum(svd.explained_variance_ratio_)),
        "heldout_target_stages": sorted(heldout_targets),
        "n_train_pairs": int(len(train_pairs.x)),
        "n_heldout_pairs": int(len(heldout_pairs.x)),
        "model_mean_metrics": metrics.groupby("model").mean(numeric_only=True).reset_index().to_dict(orient="records"),
        "mlp_training_loss": training_history,
        "context_residual_mlp_training_loss": context_residual_history,
        "context_residual_mlp_context": context_residual_info,
        "outputs": {
            "metrics": rel(tables_dir / "cell_level_transition_metrics.csv"),
            "pair_manifest": rel(tables_dir / "cell_level_pair_manifest.csv"),
            "figure": rel(figures_dir / "cell_level_transition_baseline_mse.png"),
            "model_artifacts": rel(model_dir),
        },
    }
    write_json(results_dir / "training_summary.json", summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
