from __future__ import annotations

import numpy as np


def adaptive_sigma(entropy: np.ndarray, density: np.ndarray, growth: np.ndarray, base: float = 0.04) -> np.ndarray:
    entropy = np.nan_to_num(np.asarray(entropy, dtype=float), nan=0.5)
    density = np.nan_to_num(np.asarray(density, dtype=float), nan=0.5)
    growth = np.nan_to_num(np.asarray(growth, dtype=float), nan=1.0)
    sigma = base * (1.0 + 0.9 * entropy + 0.35 * density + 0.15 * np.maximum(growth - 1.0, 0.0))
    return np.clip(sigma, base * 0.4, base * 5.0)

