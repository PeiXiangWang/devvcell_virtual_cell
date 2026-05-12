from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def bootstrap_ci(values: np.ndarray, n_boot: int = 2000, seed: int = 0, alpha: float = 0.05) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    means = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        means[i] = rng.choice(values, size=values.size, replace=True).mean()
    return float(np.quantile(means, alpha / 2)), float(np.quantile(means, 1 - alpha / 2))


def paired_test(frame: pd.DataFrame, metric: str, baseline: str, challenger: str) -> dict[str, float | str]:
    sub = frame[frame["model"].isin([baseline, challenger])]
    piv = sub.pivot_table(index="seed", columns="model", values=metric, aggfunc="mean").dropna()
    if piv.empty or baseline not in piv or challenger not in piv:
        return {"metric": metric, "baseline": baseline, "challenger": challenger, "n": 0, "p_value": np.nan, "effect": np.nan}
    diff = piv[baseline].to_numpy() - piv[challenger].to_numpy()
    try:
        stat = stats.wilcoxon(piv[baseline], piv[challenger], zero_method="wilcox", alternative="greater")
        p_value = float(stat.pvalue)
    except Exception:
        p_value = float("nan")
    lo, hi = bootstrap_ci(diff, seed=17)
    return {
        "metric": metric,
        "baseline": baseline,
        "challenger": challenger,
        "n": int(piv.shape[0]),
        "baseline_mean": float(piv[baseline].mean()),
        "challenger_mean": float(piv[challenger].mean()),
        "effect_baseline_minus_challenger": float(np.mean(diff)),
        "effect_ci_low": lo,
        "effect_ci_high": hi,
        "p_value": p_value,
    }


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    p = np.asarray([np.nan if v is None else v for v in p_values], dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    valid = np.where(np.isfinite(p))[0]
    if valid.size == 0:
        return out.tolist()
    order = valid[np.argsort(p[valid])]
    ranked = p[order]
    q = ranked * valid.size / (np.arange(valid.size) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out[order] = np.clip(q, 0, 1)
    return out.tolist()

