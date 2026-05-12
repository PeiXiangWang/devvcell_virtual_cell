from __future__ import annotations

import numpy as np


def pheromone_delta(z: np.ndarray, fate_probs: np.ndarray, fate_centers: np.ndarray, strength: float = 0.06) -> np.ndarray:
    if fate_probs.size == 0 or fate_centers.size == 0:
        return np.zeros_like(z)
    attractor = fate_probs @ fate_centers
    return strength * (attractor - z)

