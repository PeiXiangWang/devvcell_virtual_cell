from __future__ import annotations

import math
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import KNeighborsClassifier


def _sample_rows(x: np.ndarray, max_n: int, seed: int) -> np.ndarray:
    if x.shape[0] <= max_n:
        return np.asarray(x, dtype=float)
    rng = np.random.default_rng(seed)
    idx = rng.choice(x.shape[0], size=max_n, replace=False)
    return np.asarray(x[idx], dtype=float)


def mmd_rbf(x: np.ndarray, y: np.ndarray, gamma: float | None = None, max_n: int = 700, seed: int = 0) -> float:
    x = _sample_rows(np.asarray(x), max_n, seed)
    y = _sample_rows(np.asarray(y), max_n, seed + 1)
    xy = np.vstack([x, y])
    if gamma is None:
        d = pairwise_distances(_sample_rows(xy, min(500, xy.shape[0]), seed + 2), metric="sqeuclidean")
        med = np.median(d[d > 0]) if np.any(d > 0) else 1.0
        gamma = 1.0 / max(med, 1e-8)
    kxx = np.exp(-gamma * pairwise_distances(x, x, metric="sqeuclidean"))
    kyy = np.exp(-gamma * pairwise_distances(y, y, metric="sqeuclidean"))
    kxy = np.exp(-gamma * pairwise_distances(x, y, metric="sqeuclidean"))
    return float(kxx.mean() + kyy.mean() - 2.0 * kxy.mean())


def energy_distance(x: np.ndarray, y: np.ndarray, max_n: int = 700, seed: int = 0) -> float:
    x = _sample_rows(np.asarray(x), max_n, seed)
    y = _sample_rows(np.asarray(y), max_n, seed + 1)
    dxy = cdist(x, y).mean()
    dxx = cdist(x, x).mean()
    dyy = cdist(y, y).mean()
    return float(2.0 * dxy - dxx - dyy)


def sinkhorn_distance(x: np.ndarray, y: np.ndarray, reg: float = 0.08, max_n: int = 450, seed: int = 0) -> float:
    x = _sample_rows(np.asarray(x), max_n, seed)
    y = _sample_rows(np.asarray(y), max_n, seed + 1)
    cost = pairwise_distances(x, y, metric="sqeuclidean")
    med = float(np.median(cost[cost > 0])) if np.any(cost > 0) else 1.0
    cost = cost / max(med, 1e-8)
    a = np.full(x.shape[0], 1.0 / x.shape[0])
    b = np.full(y.shape[0], 1.0 / y.shape[0])
    kernel = np.exp(-cost / max(reg, 1e-8))
    kernel = np.maximum(kernel, 1e-300)
    u = np.ones_like(a)
    v = np.ones_like(b)
    for _ in range(600):
        u_prev = u.copy()
        u = a / np.maximum(kernel @ v, 1e-300)
        v = b / np.maximum(kernel.T @ u, 1e-300)
        if np.max(np.abs(u - u_prev)) < 1e-7:
            break
    plan = (u[:, None] * kernel) * v[None, :]
    plan = plan / np.maximum(plan.sum(), 1e-300)
    return float(np.sum(plan * cost))


def composition_rmse(pred_labels: Iterable[object], true_labels: Iterable[object]) -> float:
    pred = pd.Series(list(pred_labels), dtype="object")
    true = pd.Series(list(true_labels), dtype="object")
    labels = sorted(set(pred.dropna().astype(str)) | set(true.dropna().astype(str)))
    if not labels:
        return float("nan")
    pred_freq = pred.astype(str).value_counts(normalize=True).reindex(labels, fill_value=0.0)
    true_freq = true.astype(str).value_counts(normalize=True).reindex(labels, fill_value=0.0)
    return float(np.sqrt(np.mean((pred_freq.to_numpy() - true_freq.to_numpy()) ** 2)))


def knn_two_sample_accuracy(x: np.ndarray, y: np.ndarray, k: int = 15, max_n: int = 700, seed: int = 0) -> float:
    x = _sample_rows(np.asarray(x), max_n, seed)
    y = _sample_rows(np.asarray(y), max_n, seed + 1)
    z = np.vstack([x, y])
    labels = np.r_[np.zeros(x.shape[0], dtype=int), np.ones(y.shape[0], dtype=int)]
    rng = np.random.default_rng(seed)
    order = rng.permutation(z.shape[0])
    split = int(0.7 * z.shape[0])
    train_idx, test_idx = order[:split], order[split:]
    clf = KNeighborsClassifier(n_neighbors=min(k, max(1, split - 1)))
    clf.fit(z[train_idx], labels[train_idx])
    return float(clf.score(z[test_idx], labels[test_idx]))


def fate_entropy(prob: np.ndarray) -> np.ndarray:
    p = np.clip(np.asarray(prob, dtype=float), 1e-12, 1.0)
    ent = -np.sum(p * np.log(p), axis=1)
    denom = math.log(p.shape[1]) if p.ndim == 2 and p.shape[1] > 1 else 1.0
    return ent / denom
