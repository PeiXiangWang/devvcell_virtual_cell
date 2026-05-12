from __future__ import annotations

import numpy as np
from sklearn.neighbors import NearestNeighbors


def local_density(z: np.ndarray, k: int = 20) -> np.ndarray:
    if z.shape[0] <= 2:
        return np.zeros(z.shape[0])
    k = min(k, z.shape[0] - 1)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(z)
    dist, _ = nn.kneighbors(z)
    radius = np.mean(dist[:, 1:], axis=1)
    dens = 1.0 / np.maximum(radius, 1e-6)
    lo, hi = np.quantile(dens, [0.05, 0.95])
    return np.clip((dens - lo) / max(hi - lo, 1e-8), 0.0, 1.0)


def swarm_delta(z: np.ndarray, velocity: np.ndarray, labels: np.ndarray, strength: float = 0.18, k: int = 15) -> np.ndarray:
    if z.shape[0] <= 2:
        return np.zeros_like(z)
    k = min(k, z.shape[0] - 1)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(z)
    dist, idx = nn.kneighbors(z)
    neigh = idx[:, 1:]
    delta = np.zeros_like(z)
    for i in range(z.shape[0]):
        nidx = neigh[i]
        same = labels[nidx] == labels[i]
        if np.any(same):
            center = z[nidx[same]].mean(axis=0)
            align = velocity[nidx[same]].mean(axis=0) if velocity is not None else 0
            delta[i] += 0.55 * (center - z[i]) + 0.35 * align
        close = nidx[dist[i, 1:] < np.quantile(dist[:, 1:].ravel(), 0.15)]
        if close.size:
            repel = (z[i] - z[close]).mean(axis=0)
            norm = np.linalg.norm(repel)
            if norm > 0:
                delta[i] += 0.25 * repel / norm
    return strength * delta

