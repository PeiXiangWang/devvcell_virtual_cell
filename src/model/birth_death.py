from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class BirthDeathResult:
    keep_mask: np.ndarray
    daughter_parent_indices: np.ndarray
    event_rows: list[dict]


def stochastic_birth_death(
    birth_hazard: np.ndarray,
    death_hazard: np.ndarray,
    dt: float,
    seed: int,
    time_value: float,
    labels: np.ndarray,
) -> BirthDeathResult:
    """Sample continuous-time birth/death events for a finite agent population."""
    rng = np.random.default_rng(seed)
    birth_hazard = np.clip(np.nan_to_num(birth_hazard, nan=0.0), 0.0, 20.0)
    death_hazard = np.clip(np.nan_to_num(death_hazard, nan=0.0), 0.0, 20.0)
    p_birth = 1.0 - np.exp(-birth_hazard * dt)
    p_death = 1.0 - np.exp(-death_hazard * dt)
    death = rng.random(death_hazard.size) < p_death
    birth = (rng.random(birth_hazard.size) < p_birth) & (~death)
    daughters = np.where(birth)[0]
    rows = []
    for idx in np.where(death)[0]:
        rows.append({"time": float(time_value), "event": "death", "parent_index": int(idx), "lineage": str(labels[idx])})
    for idx in daughters:
        rows.append({"time": float(time_value), "event": "birth", "parent_index": int(idx), "lineage": str(labels[idx])})
    return BirthDeathResult(keep_mask=~death, daughter_parent_indices=daughters.astype(int), event_rows=rows)


def calibrate_count_hazard(growth: np.ndarray, cell_cycle: np.ndarray, density: np.ndarray, entropy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    growth = np.nan_to_num(growth, nan=1.0)
    cell_cycle = np.nan_to_num(cell_cycle, nan=0.0)
    density = np.nan_to_num(density, nan=0.5)
    entropy = np.nan_to_num(entropy, nan=0.5)
    birth = np.maximum(np.log(np.maximum(growth, 1e-3)), 0.0) + 0.08 * cell_cycle + 0.05 * entropy
    death = np.maximum(-np.log(np.maximum(growth, 1e-3)), 0.0) + 0.04 * density
    return birth, death
