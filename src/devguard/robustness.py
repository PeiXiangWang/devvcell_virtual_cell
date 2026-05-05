"""Robustness summaries for DevGuard cell-classification outputs."""

from __future__ import annotations

import numpy as np
import pandas as pd

CLASS_ORDER = [
    "within_stage_normal",
    "developmental_delay",
    "developmental_acceleration",
    "fate_deviation",
    "abnormal_off_normal",
]


def class_order(frame: pd.DataFrame, class_col: str = "normality_class") -> list[str]:
    observed = set(frame[class_col].dropna().astype(str))
    return [item for item in CLASS_ORDER if item in observed] + sorted(observed - set(CLASS_ORDER))


def class_fraction_table(
    frame: pd.DataFrame,
    group_cols: list[str],
    *,
    class_col: str = "normality_class",
    weight_col: str | None = None,
) -> pd.DataFrame:
    """Summarize class fractions, optionally using cell weights."""

    data = frame.copy()
    if weight_col is None:
        data["_devguard_weight"] = 1.0
    else:
        data["_devguard_weight"] = pd.to_numeric(data[weight_col], errors="coerce").fillna(0.0)
    counts = (
        data.groupby(group_cols + [class_col], dropna=False, observed=True)["_devguard_weight"]
        .sum()
        .reset_index(name="weighted_cells")
    )
    totals = counts.groupby(group_cols, dropna=False, observed=True)["weighted_cells"].transform("sum")
    counts["fraction"] = counts["weighted_cells"] / totals.where(totals > 0, np.nan)
    counts = counts.rename(columns={class_col: "normality_class"})
    return counts.sort_values(group_cols + ["normality_class"]).reset_index(drop=True)


def sample_class_summary(
    frame: pd.DataFrame,
    *,
    condition_col: str = "perturbation_name",
    sample_col: str = "sample_id",
    class_col: str = "normality_class",
) -> pd.DataFrame:
    return class_fraction_table(frame, [condition_col, sample_col], class_col=class_col).rename(
        columns={"weighted_cells": "n_cells"}
    )


def bootstrap_sample_class_fractions(
    frame: pd.DataFrame,
    *,
    condition_col: str = "perturbation_name",
    sample_col: str = "sample_id",
    class_col: str = "normality_class",
    n_bootstrap: int = 1000,
    seed: int = 42,
) -> pd.DataFrame:
    """Resample sample units and estimate class-fraction intervals."""

    rng = np.random.default_rng(seed)
    rows = []
    classes = class_order(frame, class_col)
    for condition, group in frame.groupby(condition_col, dropna=False, observed=True):
        sample_labels = group[sample_col].astype("string").fillna("NA").astype(str)
        units = np.asarray(sorted(sample_labels.unique()))
        observed_counts = group[class_col].value_counts()
        observed_total = max(int(group.shape[0]), 1)
        for cls in classes:
            boot_values = []
            if units.size > 1:
                for _ in range(n_bootstrap):
                    sampled_units = rng.choice(units, size=units.size, replace=True)
                    sample = pd.concat([group.loc[sample_labels.eq(unit)] for unit in sampled_units], axis=0)
                    boot_values.append(float((sample[class_col] == cls).mean()) if sample.shape[0] else np.nan)
                bootstrap_level = "sample"
            else:
                positions = np.arange(group.shape[0])
                for _ in range(n_bootstrap):
                    sampled = group.iloc[rng.choice(positions, size=positions.size, replace=True)]
                    boot_values.append(float((sampled[class_col] == cls).mean()) if sampled.shape[0] else np.nan)
                bootstrap_level = "cell"
            sample_fractions = (
                group.assign(_sample=sample_labels)
                .groupby("_sample", dropna=False, observed=True)[class_col]
                .apply(lambda values: float((values == cls).mean()))
            )
            values = np.asarray(boot_values, dtype=float)
            rows.append(
                {
                    condition_col: condition,
                    "normality_class": cls,
                    "n_cells": int(group.shape[0]),
                    "n_sample_units": int(units.size),
                    "observed_cell_fraction": float(observed_counts.get(cls, 0) / observed_total),
                    "sample_fraction_mean": float(sample_fractions.mean()) if not sample_fractions.empty else np.nan,
                    "sample_fraction_median": float(sample_fractions.median()) if not sample_fractions.empty else np.nan,
                    "sample_fraction_min": float(sample_fractions.min()) if not sample_fractions.empty else np.nan,
                    "sample_fraction_max": float(sample_fractions.max()) if not sample_fractions.empty else np.nan,
                    "bootstrap_fraction_mean": float(np.nanmean(values)),
                    "bootstrap_ci_lower": float(np.nanquantile(values, 0.025)),
                    "bootstrap_ci_upper": float(np.nanquantile(values, 0.975)),
                    "bootstrap_level": bootstrap_level,
                    "n_bootstrap": int(n_bootstrap),
                }
            )
    return pd.DataFrame(rows)


def balanced_downsample_class_fractions(
    frame: pd.DataFrame,
    *,
    condition_col: str = "perturbation_name",
    strata_cols: list[str] | None = None,
    class_col: str = "normality_class",
    n_iterations: int = 200,
    max_cells_per_stratum: int | None = 200,
    seed: int = 42,
) -> pd.DataFrame:
    """Repeatedly downsample conditions to matched cell counts within strata."""

    strata_cols = strata_cols or []
    rng = np.random.default_rng(seed)
    classes = class_order(frame, class_col)
    conditions = sorted(frame[condition_col].dropna().astype(str).unique())
    data = frame.copy()
    data[condition_col] = data[condition_col].astype(str)

    if strata_cols:
        count_table = data.groupby([condition_col] + strata_cols, dropna=False, observed=True).size().reset_index(name="n")
        stratum_counts = count_table.pivot_table(index=strata_cols, columns=condition_col, values="n", fill_value=0, aggfunc="sum")
        for condition in conditions:
            if condition not in stratum_counts.columns:
                stratum_counts[condition] = 0
        common = stratum_counts[(stratum_counts[conditions] > 0).all(axis=1)].copy()
        common["draw_n"] = common[conditions].min(axis=1)
        if max_cells_per_stratum is not None:
            common["draw_n"] = np.minimum(common["draw_n"], int(max_cells_per_stratum))
        common = common[common["draw_n"] > 0]
        strata_records = common.reset_index().to_dict("records")
    else:
        draw_n = data.groupby(condition_col).size().min()
        if max_cells_per_stratum is not None:
            draw_n = min(int(draw_n), int(max_cells_per_stratum))
        strata_records = [{"draw_n": int(draw_n)}] if draw_n > 0 else []

    rows = []
    for iteration in range(n_iterations):
        pieces = []
        for condition in conditions:
            condition_data = data[data[condition_col] == condition]
            for record in strata_records:
                draw_n = int(record["draw_n"])
                if strata_cols:
                    mask = np.ones(condition_data.shape[0], dtype=bool)
                    for column in strata_cols:
                        mask &= condition_data[column].astype(str).to_numpy() == str(record[column])
                    pool = condition_data.loc[mask]
                else:
                    pool = condition_data
                if pool.shape[0] < draw_n or draw_n <= 0:
                    continue
                pieces.append(pool.iloc[rng.choice(pool.shape[0], size=draw_n, replace=False)])
        sampled = pd.concat(pieces, axis=0) if pieces else data.iloc[[]]
        for condition, condition_data in sampled.groupby(condition_col, dropna=False, observed=True):
            total = max(int(condition_data.shape[0]), 1)
            counts = condition_data[class_col].value_counts()
            for cls in classes:
                rows.append(
                    {
                        "iteration": iteration,
                        condition_col: condition,
                        "normality_class": cls,
                        "fraction": float(counts.get(cls, 0) / total),
                        "n_cells": int(condition_data.shape[0]),
                        "n_strata": int(len(strata_records)),
                    }
                )
    return pd.DataFrame(rows)


def summarize_downsample_iterations(frame: pd.DataFrame, *, condition_col: str = "perturbation_name") -> pd.DataFrame:
    group_cols = [condition_col, "normality_class"]
    summary = (
        frame.groupby(group_cols, dropna=False, observed=True)
        .agg(
            fraction_mean=("fraction", "mean"),
            fraction_median=("fraction", "median"),
            fraction_ci_lower=("fraction", lambda values: float(np.quantile(values, 0.025))),
            fraction_ci_upper=("fraction", lambda values: float(np.quantile(values, 0.975))),
            n_iterations=("iteration", "nunique"),
            mean_cells=("n_cells", "mean"),
            n_strata=("n_strata", "max"),
        )
        .reset_index()
    )
    return summary
