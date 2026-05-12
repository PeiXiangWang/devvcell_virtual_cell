from __future__ import annotations

import numpy as np
from sklearn.neighbors import NearestNeighbors


class MemoryField:
    """Fate-specific phenomenological field on the current latent kNN graph."""

    def __init__(self, n_fates: int, decay: float = 0.08, diffusion: float = 0.25):
        self.n_fates = int(max(n_fates, 1))
        self.decay = float(decay)
        self.diffusion = float(diffusion)
        self.values: np.ndarray | None = None

    def ensure(self, n_agents: int) -> None:
        if self.values is None or self.values.shape[0] != n_agents:
            self.values = np.zeros((n_agents, self.n_fates), dtype=float)

    def step(self, z: np.ndarray, fate_probs: np.ndarray, dt: float, k: int = 12) -> None:
        self.ensure(z.shape[0])
        fate = fate_probs if fate_probs.size else np.full((z.shape[0], self.n_fates), 1.0 / self.n_fates)
        if fate.shape[1] != self.n_fates:
            fate = np.resize(fate, (z.shape[0], self.n_fates))
        self.values += dt * fate
        self.values *= max(0.0, 1.0 - self.decay * dt)
        if z.shape[0] > 2:
            k = min(k, z.shape[0] - 1)
            nn = NearestNeighbors(n_neighbors=k + 1).fit(z)
            neigh = nn.kneighbors(z, return_distance=False)[:, 1:]
            smooth = self.values[neigh].mean(axis=1)
            self.values += self.diffusion * dt * (smooth - self.values)

    def gradient_delta(self, z: np.ndarray, fate_probs: np.ndarray, strength: float = 0.05, k: int = 12) -> np.ndarray:
        self.ensure(z.shape[0])
        if z.shape[0] <= 2:
            return np.zeros_like(z)
        fate = fate_probs if fate_probs.size else np.full((z.shape[0], self.n_fates), 1.0 / self.n_fates)
        k = min(k, z.shape[0] - 1)
        nn = NearestNeighbors(n_neighbors=k + 1).fit(z)
        neigh = nn.kneighbors(z, return_distance=False)[:, 1:]
        agent_field = np.sum(self.values * fate[:, : self.n_fates], axis=1)
        delta = np.zeros_like(z)
        for i in range(z.shape[0]):
            weights = np.maximum(agent_field[neigh[i]] - agent_field[i], 0.0)
            if weights.sum() > 0:
                target = np.average(z[neigh[i]], axis=0, weights=weights)
                delta[i] = target - z[i]
        return strength * delta
