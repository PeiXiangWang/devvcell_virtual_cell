"""Run a small external perturbation-response benchmark on scPerturb.

The default benchmark uses Datlinger/Bock 2021 from scPerturb. It trains on
guide suffix ``_1`` and evaluates guide suffix ``_2`` for the same perturbed
genes and stimulation contexts, so the test split measures guide transfer
rather than random cell-level interpolation.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import StandardScaler


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from devvcell.io import load_json, resolve_project_path, write_json  # noqa: E402
from train_cell_transition_baseline import as_csr_float32, sinkhorn_barycentric_targets  # noqa: E402


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
class PerturbationPairSet:
    x: np.ndarray
    y: np.ndarray
    gene: np.ndarray
    context: np.ndarray
    perturbation_label: np.ndarray
    split: np.ndarray
    condition: np.ndarray


class GeneContextResidualMLP(nn.Module):
    def __init__(self, input_dim: int, latent_dim: int, hidden_dim: int, dropout: float):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, latent_dim),
        )

    def forward(self, features):
        return self.net(features)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/external_perturbation_benchmark.json")
    return parser.parse_args()


def perturbation_gene(label: str, suffixes: tuple[str, str]) -> str | None:
    if label.lower() in {"control", "nan", "*", "none"}:
        return None
    for suffix in suffixes:
        if label.endswith(suffix):
            return label[: -len(suffix)]
    match = re.match(r"^(.+)_\d+$", label)
    return match.group(1) if match else None


def choose_feature_genes(var: pd.DataFrame, var_names: pd.Index, perturb_genes: list[str], max_genes: int) -> list[str]:
    score_col = "ncounts" if "ncounts" in var.columns else "ncells" if "ncells" in var.columns else None
    if score_col is not None:
        scores = pd.to_numeric(var[score_col], errors="coerce").fillna(0.0)
        ranked = scores.sort_values(ascending=False).index.astype(str).tolist()
    else:
        ranked = var_names.astype(str).tolist()

    selected: list[str] = []
    seen: set[str] = set()
    for gene in perturb_genes + ranked:
        if gene in var_names and gene not in seen:
            selected.append(gene)
            seen.add(gene)
        if len(selected) >= max_genes:
            break
    return selected


def fit_latent(X, latent_dim: int, seed: int) -> tuple[np.ndarray, TruncatedSVD, StandardScaler]:
    X = as_csr_float32(X)
    n_components = min(latent_dim, X.shape[1] - 1, X.shape[0] - 1)
    svd = TruncatedSVD(n_components=n_components, random_state=seed)
    raw = svd.fit_transform(X)
    scaler = StandardScaler()
    latent = scaler.fit_transform(raw).astype(np.float32)
    return latent, svd, scaler


def build_pair_targets(source: np.ndarray, target_pool: np.ndarray, cfg: dict) -> tuple[np.ndarray, dict[str, object]]:
    method = str(cfg["model"].get("pairing_method", "sinkhorn")).lower()
    if method == "sinkhorn":
        target, info = sinkhorn_barycentric_targets(
            source,
            target_pool,
            epsilon=float(cfg["model"].get("sinkhorn_epsilon", 0.12)),
            iterations=int(cfg["model"].get("sinkhorn_iterations", 100)),
        )
        return target, {"pairing_method": "sinkhorn", **info}
    if method in {"random", "sample"}:
        return target_pool[: len(source)].astype(np.float32), {"pairing_method": "random"}
    raise ValueError(f"Unsupported pairing method: {method!r}")


def build_pairs(
    latent: np.ndarray,
    obs: pd.DataFrame,
    genes: list[str],
    suffix: str,
    split: str,
    cfg: dict,
    rng: np.random.Generator,
) -> tuple[PerturbationPairSet, pd.DataFrame]:
    data_cfg = cfg["data"]
    model_cfg = cfg["model"]
    perturb_col = data_cfg["perturbation_column"]
    context_col = data_cfg["context_column"]
    control_label = str(data_cfg["control_label"])
    min_cells = int(data_cfg["min_cells_per_condition"])
    max_pairs = int(model_cfg["max_pairs_per_condition"])
    target_pool_multiplier = max(1, int(model_cfg.get("target_pool_multiplier", 2)))

    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    gene_rows: list[str] = []
    context_rows: list[str] = []
    label_rows: list[str] = []
    split_rows: list[str] = []
    condition_rows: list[str] = []
    manifest_rows: list[dict[str, object]] = []

    contexts = sorted(obs[context_col].astype(str).unique())
    for gene in genes:
        label = f"{gene}{suffix}"
        for context in contexts:
            ctrl_idx = obs.index[
                (obs[perturb_col].astype(str) == control_label) & (obs[context_col].astype(str) == context)
            ].to_numpy()
            tgt_idx = obs.index[
                (obs[perturb_col].astype(str) == label) & (obs[context_col].astype(str) == context)
            ].to_numpy()
            pair_count = min(max_pairs, len(ctrl_idx), len(tgt_idx))
            if pair_count < min_cells:
                continue

            src_sample = rng.choice(ctrl_idx, size=pair_count, replace=False)
            pool_size = min(len(tgt_idx), pair_count * target_pool_multiplier)
            tgt_pool = rng.choice(tgt_idx, size=pool_size, replace=False)
            source = latent[src_sample]
            target, target_info = build_pair_targets(source, latent[tgt_pool], cfg)
            condition = f"{gene}|{context}|{label}"

            xs.append(source)
            ys.append(target)
            gene_rows.extend([gene] * pair_count)
            context_rows.extend([context] * pair_count)
            label_rows.extend([label] * pair_count)
            split_rows.extend([split] * pair_count)
            condition_rows.extend([condition] * pair_count)
            manifest_rows.append(
                {
                    "split": split,
                    "gene": gene,
                    "context": context,
                    "perturbation_label": label,
                    "pairs": int(pair_count),
                    "control_cells_available": int(len(ctrl_idx)),
                    "target_cells_available": int(len(tgt_idx)),
                    "target_pool_cells": int(pool_size),
                    **target_info,
                }
            )

    if not xs:
        raise RuntimeError(f"No {split} pairs were built. Check labels and min_cells_per_condition.")
    pair_set = PerturbationPairSet(
        x=np.vstack(xs).astype(np.float32),
        y=np.vstack(ys).astype(np.float32),
        gene=np.asarray(gene_rows, dtype=object),
        context=np.asarray(context_rows, dtype=object),
        perturbation_label=np.asarray(label_rows, dtype=object),
        split=np.asarray(split_rows, dtype=object),
        condition=np.asarray(condition_rows, dtype=object),
    )
    return pair_set, pd.DataFrame(manifest_rows)


def context_features(pair_set: PerturbationPairSet, gene_order: list[str], context_order: list[str]) -> np.ndarray:
    gene_pos = {name: idx for idx, name in enumerate(gene_order)}
    context_pos = {name: idx for idx, name in enumerate(context_order)}
    gene_onehot = np.zeros((len(pair_set.x), len(gene_order)), dtype=np.float32)
    context_onehot = np.zeros((len(pair_set.x), len(context_order)), dtype=np.float32)
    for idx, value in enumerate(pair_set.gene):
        gene_onehot[idx, gene_pos[str(value)]] = 1.0
    for idx, value in enumerate(pair_set.context):
        context_onehot[idx, context_pos[str(value)]] = 1.0
    return np.hstack([pair_set.x, gene_onehot, context_onehot]).astype(np.float32)


def train_mlp(features: np.ndarray, deltas: np.ndarray, cfg: dict, seed: int):
    if torch is None:
        return None, {"enabled": False, "reason": "torch_not_available"}

    mlp_cfg = cfg["model"]["mlp"]
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(features))
    val_size = max(1, int(len(order) * float(mlp_cfg.get("validation_fraction", 0.15))))
    val_idx = order[:val_size]
    train_idx = order[val_size:]

    torch.manual_seed(seed)
    model = GeneContextResidualMLP(
        input_dim=features.shape[1],
        latent_dim=deltas.shape[1],
        hidden_dim=int(mlp_cfg["hidden_dim"]),
        dropout=float(mlp_cfg.get("dropout", 0.0)),
    )
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=float(mlp_cfg["learning_rate"]),
        weight_decay=float(mlp_cfg.get("weight_decay", 0.0)),
    )
    loss_fn = nn.MSELoss()
    train_ds = TensorDataset(torch.from_numpy(features[train_idx]), torch.from_numpy(deltas[train_idx]))
    val_x = torch.from_numpy(features[val_idx])
    val_y = torch.from_numpy(deltas[val_idx])
    loader = DataLoader(train_ds, batch_size=int(mlp_cfg["batch_size"]), shuffle=True)

    best_loss = float("inf")
    best_state = None
    stale = 0
    history: list[dict[str, float]] = []
    patience = int(mlp_cfg.get("patience", 20))
    min_delta = float(mlp_cfg.get("min_delta", 1e-6))

    for epoch in range(1, int(mlp_cfg["epochs"]) + 1):
        model.train()
        batch_losses = []
        for xb, yb in loader:
            optimizer.zero_grad()
            loss = loss_fn(model(xb), yb)
            loss.backward()
            optimizer.step()
            batch_losses.append(float(loss.detach().cpu()))

        model.eval()
        with torch.no_grad():
            val_loss = float(loss_fn(model(val_x), val_y).cpu())
        history.append({"epoch": int(epoch), "train_loss": float(np.mean(batch_losses)), "validation_loss": val_loss})
        if val_loss < best_loss - min_delta:
            best_loss = val_loss
            best_state = {name: tensor.detach().cpu().clone() for name, tensor in model.state_dict().items()}
            stale = 0
        else:
            stale += 1
        if patience > 0 and stale >= patience:
            break

    if best_state is not None:
        model.load_state_dict(best_state)
    return model, {
        "enabled": True,
        "epochs_ran": int(len(history)),
        "best_validation_loss": float(best_loss),
        "history": history,
    }


def predict_mlp(model, features: np.ndarray, x: np.ndarray) -> np.ndarray:
    model.eval()
    with torch.no_grad():
        delta = model(torch.from_numpy(features)).cpu().numpy().astype(np.float32)
    return x + delta


def fit_gene_context_shift(train_pairs: PerturbationPairSet) -> tuple[dict[tuple[str, str], np.ndarray], np.ndarray]:
    deltas = train_pairs.y - train_pairs.x
    global_delta = deltas.mean(axis=0)
    by_key: dict[tuple[str, str], np.ndarray] = {}
    for gene in sorted(set(train_pairs.gene)):
        for context in sorted(set(train_pairs.context)):
            mask = (train_pairs.gene == gene) & (train_pairs.context == context)
            if np.any(mask):
                by_key[(str(gene), str(context))] = deltas[mask].mean(axis=0)
    return by_key, global_delta


def predict_gene_context_shift(pair_set: PerturbationPairSet, by_key: dict[tuple[str, str], np.ndarray], global_delta: np.ndarray) -> np.ndarray:
    pred = np.empty_like(pair_set.x)
    for gene in sorted(set(pair_set.gene)):
        for context in sorted(set(pair_set.context)):
            mask = (pair_set.gene == gene) & (pair_set.context == context)
            if np.any(mask):
                pred[mask] = pair_set.x[mask] + by_key.get((str(gene), str(context)), global_delta)
    return pred


def rbf_mmd(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) == 0 or len(y) == 0:
        return float("nan")
    joined = np.vstack([x, y])
    diffs = joined[:, None, :] - joined[None, :, :]
    dist2 = np.sum(diffs * diffs, axis=2)
    positive = dist2[dist2 > 0]
    gamma = 1.0 / max(float(np.median(positive)) if len(positive) else 1.0, 1e-6)
    kxx = np.exp(-gamma * dist2[: len(x), : len(x)]).mean()
    kyy = np.exp(-gamma * dist2[len(x) :, len(x) :]).mean()
    kxy = np.exp(-gamma * dist2[: len(x), len(x) :]).mean()
    return float(kxx + kyy - 2.0 * kxy)


def cosine_or_corr(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 0 or not np.isfinite(denom):
        return float("nan")
    return float(np.dot(a, b) / denom)


def evaluate_predictions(pair_set: PerturbationPairSet, predictions: dict[str, np.ndarray]) -> tuple[pd.DataFrame, pd.DataFrame]:
    condition_rows = []
    for model_name, pred in predictions.items():
        for condition in sorted(set(pair_set.condition)):
            mask = pair_set.condition == condition
            x = pair_set.x[mask]
            y = pair_set.y[mask]
            p = pred[mask]
            condition_rows.append(
                {
                    "model": model_name,
                    "condition": condition,
                    "gene": str(pair_set.gene[mask][0]),
                    "context": str(pair_set.context[mask][0]),
                    "perturbation_label": str(pair_set.perturbation_label[mask][0]),
                    "n_pairs": int(mask.sum()),
                    "pair_latent_mse": float(mean_squared_error(y, p)),
                    "centroid_latent_mse": float(np.mean((y.mean(axis=0) - p.mean(axis=0)) ** 2)),
                    "effect_cosine": cosine_or_corr(p.mean(axis=0) - x.mean(axis=0), y.mean(axis=0) - x.mean(axis=0)),
                    "rbf_mmd": rbf_mmd(p, y),
                }
            )
    condition_table = pd.DataFrame(condition_rows)
    summary = (
        condition_table.groupby("model", as_index=False)
        .agg(
            n_conditions=("condition", "nunique"),
            mean_pair_latent_mse=("pair_latent_mse", "mean"),
            mean_centroid_latent_mse=("centroid_latent_mse", "mean"),
            mean_effect_cosine=("effect_cosine", "mean"),
            mean_rbf_mmd=("rbf_mmd", "mean"),
        )
        .sort_values("mean_pair_latent_mse")
    )
    return condition_table, summary


def plot_summary(summary: pd.DataFrame, figures_dir: Path) -> None:
    labels = {
        "identity": "Identity",
        "global_mean_shift": "Global shift",
        "gene_context_shift": "Guide-transfer shift",
        "gene_context_ridge_residual": "Context ridge residual",
        "gene_context_mlp_residual": "Context neural residual",
    }
    ordered = summary.sort_values("mean_pair_latent_mse").copy()
    ordered["label"] = ordered["model"].map(labels).fillna(ordered["model"])
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    colors = ["#4e79a7", "#59a14f", "#f28e2b", "#b07aa1", "#9c755f"][: len(ordered)]
    ax.bar(ordered["label"], ordered["mean_pair_latent_mse"], color=colors)
    ax.set_ylabel("Heldout guide pseudo-target latent MSE")
    ax.set_xlabel("")
    ax.set_title("External scPerturb guide-transfer benchmark")
    ax.tick_params(axis="x", rotation=18)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(figures_dir / "external_perturbation_benchmark_mse.png", dpi=220)
    plt.close(fig)


def json_safe(value):
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    if isinstance(value, tuple):
        return [json_safe(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        numeric = float(value)
        return numeric if np.isfinite(numeric) else None
    return value


def main() -> None:
    args = parse_args()
    cfg = load_json(resolve_project_path(args.config))
    seed = int(cfg["model"]["seed"])
    rng = np.random.default_rng(seed)

    h5ad_path = resolve_project_path(cfg["data"]["h5ad_path"])
    out_dir = resolve_project_path(cfg["outputs"]["output_dir"])
    tables_dir = out_dir / "tables"
    figures_dir = out_dir / "figures"
    models_dir = out_dir / "models"
    for directory in [tables_dir, figures_dir, models_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(h5ad_path, backed="r")
    obs = adata.obs.copy().reset_index(drop=True)
    var = adata.var.copy()
    var_names = pd.Index(adata.var_names.astype(str))
    perturb_col = cfg["data"]["perturbation_column"]
    suffixes = (str(cfg["data"]["train_suffix"]), str(cfg["data"]["test_suffix"]))

    label_counts = obs[perturb_col].astype(str).value_counts()
    genes = sorted(
        {
            gene
            for label in label_counts.index
            for gene in [perturbation_gene(str(label), suffixes)]
            if gene is not None
        }
    )
    train_suffix, test_suffix = suffixes
    eligible_genes = [
        gene
        for gene in genes
        if f"{gene}{train_suffix}" in label_counts.index and f"{gene}{test_suffix}" in label_counts.index
    ]
    selected_genes = choose_feature_genes(var, var_names, eligible_genes, int(cfg["model"]["max_genes"]))
    selected_positions = [int(var_names.get_loc(gene)) for gene in selected_genes]
    X = adata[:, selected_positions].X
    X = as_csr_float32(X)
    adata.file.close()

    latent, svd, scaler = fit_latent(X, int(cfg["model"]["latent_dim"]), seed)
    joblib.dump(svd, models_dir / "external_perturbation_svd.joblib")
    joblib.dump(scaler, models_dir / "external_perturbation_scaler.joblib")

    train_pairs, train_manifest = build_pairs(latent, obs, eligible_genes, train_suffix, "train", cfg, rng)
    test_pairs, test_manifest = build_pairs(latent, obs, eligible_genes, test_suffix, "heldout_guide", cfg, rng)
    pair_manifest = pd.concat([train_manifest, test_manifest], ignore_index=True)
    pair_manifest.to_csv(tables_dir / "external_perturbation_pair_manifest.csv", index=False)

    gene_order = eligible_genes
    context_order = sorted(obs[cfg["data"]["context_column"]].astype(str).unique())
    train_features = context_features(train_pairs, gene_order, context_order)
    test_features = context_features(test_pairs, gene_order, context_order)
    train_delta = train_pairs.y - train_pairs.x

    by_key, global_delta = fit_gene_context_shift(train_pairs)
    ridge = Ridge(alpha=float(cfg["model"].get("ridge_alpha", 1.0)))
    ridge.fit(train_features, train_delta)
    joblib.dump(ridge, models_dir / "external_perturbation_ridge_residual.joblib")

    predictions = {
        "identity": test_pairs.x,
        "global_mean_shift": test_pairs.x + global_delta,
        "gene_context_shift": predict_gene_context_shift(test_pairs, by_key, global_delta),
        "gene_context_ridge_residual": test_pairs.x + ridge.predict(test_features).astype(np.float32),
    }

    mlp_info = {"enabled": False}
    if bool(cfg["model"].get("mlp", {}).get("enabled", False)):
        mlp, mlp_info = train_mlp(train_features, train_delta.astype(np.float32), cfg, seed)
        if mlp is not None:
            torch.save(
                {
                    "model_state": mlp.state_dict(),
                    "gene_order": gene_order,
                    "context_order": context_order,
                    "config": cfg["model"]["mlp"],
                },
                models_dir / "external_perturbation_gene_context_mlp.pt",
            )
            predictions["gene_context_mlp_residual"] = predict_mlp(mlp, test_features, test_pairs.x)

    condition_metrics, summary = evaluate_predictions(test_pairs, predictions)
    condition_metrics.to_csv(tables_dir / "external_perturbation_condition_metrics.csv", index=False)
    summary.to_csv(tables_dir / "external_perturbation_metrics.csv", index=False)
    plot_summary(summary, figures_dir)

    result = {
        "analysis": "external_scperturb_guide_transfer_benchmark",
        "source_h5ad": str(h5ad_path.relative_to(PROJECT_ROOT)).replace("\\", "/"),
        "split_design": {
            "train_suffix": train_suffix,
            "heldout_test_suffix": test_suffix,
            "eligible_genes": eligible_genes,
            "contexts": context_order,
        },
        "latent": {
            "n_cells": int(latent.shape[0]),
            "n_feature_genes": int(len(selected_genes)),
            "latent_dim": int(latent.shape[1]),
            "explained_variance_ratio_sum": float(np.sum(svd.explained_variance_ratio_)),
        },
        "pairing": {
            "train_pairs": int(len(train_pairs.x)),
            "heldout_pairs": int(len(test_pairs.x)),
            "conditions": int(len(set(test_pairs.condition))),
            "method": cfg["model"].get("pairing_method", "sinkhorn"),
        },
        "best_model_by_pair_latent_mse": str(summary.iloc[0]["model"]),
        "model_summary": summary.to_dict(orient="records"),
        "mlp_training": mlp_info,
        "outputs": {
            "metrics": "results/external_perturbation_v1/tables/external_perturbation_metrics.csv",
            "condition_metrics": "results/external_perturbation_v1/tables/external_perturbation_condition_metrics.csv",
            "pair_manifest": "results/external_perturbation_v1/tables/external_perturbation_pair_manifest.csv",
            "figure": "results/external_perturbation_v1/figures/external_perturbation_benchmark_mse.png",
        },
    }
    safe_result = json_safe(result)
    write_json(out_dir / "external_perturbation_summary.json", safe_result)
    print(json.dumps(safe_result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
