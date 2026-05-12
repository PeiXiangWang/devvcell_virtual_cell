from __future__ import annotations

import numpy as np
import pandas as pd


def growth_resample_indices(labels: pd.Series, target_n: int, growth: np.ndarray, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    growth = np.nan_to_num(np.asarray(growth, dtype=float), nan=1.0, posinf=1.0, neginf=1.0)
    weights = np.clip(growth, 0.05, 20.0)
    if labels is not None:
        label_counts = labels.astype(str).value_counts(normalize=True)
        weights = weights * np.array([1.0 / max(label_counts.get(str(v), 1e-6), 1e-6) for v in labels.astype(str)])
    weights = weights / np.maximum(weights.sum(), 1e-12)
    return rng.choice(np.arange(weights.size), size=target_n, replace=True, p=weights)


def hazard_scores(growth: np.ndarray, cell_cycle: np.ndarray, density: np.ndarray, entropy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    growth = np.nan_to_num(growth, nan=1.0)
    cell_cycle = np.nan_to_num(cell_cycle, nan=0.0)
    density = np.nan_to_num(density, nan=0.5)
    entropy = np.nan_to_num(entropy, nan=0.5)
    birth = np.log1p(np.maximum(growth, 0)) + 0.15 * cell_cycle + 0.2 * entropy - 0.2 * density
    death = 0.1 + 0.25 * density - 0.1 * np.log1p(np.maximum(growth, 0))
    return birth, death

