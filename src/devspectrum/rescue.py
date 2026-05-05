"""In silico spectral rescue candidate ranking."""

from __future__ import annotations

import pandas as pd


def rank_rescue_candidates(
    residuals: pd.DataFrame,
    *,
    cohort: str = "Tal1_chimera",
    target_reduction: float = 0.5,
) -> tuple[pd.DataFrame, str]:
    data = residuals[residuals["cohort"].astype(str).eq(cohort)].copy()
    if data.empty:
        return pd.DataFrame(), "# Spectral Rescue Candidates\n\nNo residuals available.\n"
    module = (
        data.groupby("module_name", dropna=False, observed=True)
        .agg(
            mean_abs_residual=("absolute_spectral_residual", "mean"),
            total_abs_residual=("absolute_spectral_residual", "sum"),
            max_abs_residual=("absolute_spectral_residual", "max"),
            n_lineage_class_records=("module_name", "size"),
            affected_lineages=("lineage", lambda values: ";".join(sorted(set(values.astype(str))))),
        )
        .reset_index()
        .sort_values("total_abs_residual", ascending=False)
    )
    total = module["total_abs_residual"].sum()
    module["cumulative_abs_residual"] = module["total_abs_residual"].cumsum()
    module["cumulative_fraction_of_residual"] = module["cumulative_abs_residual"] / max(total, 1e-12)
    module["selected_for_50pct_correction"] = module["cumulative_fraction_of_residual"] <= target_reduction
    if not module.empty:
        first_above = module.index[module["cumulative_fraction_of_residual"] >= target_reduction]
        if len(first_above):
            module.loc[first_above[0], "selected_for_50pct_correction"] = True
    lines = [
        "# Tal1 Spectral Rescue Candidate Report",
        "",
        "These are in silico spectral rescue hypotheses, not validated rescue mechanisms.",
        "",
        "## Minimal module set for 50% residual coverage",
        "",
        module[module["selected_for_50pct_correction"]][
            ["module_name", "total_abs_residual", "cumulative_fraction_of_residual", "affected_lineages"]
        ].to_markdown(index=False),
        "",
        "## Top residual modules",
        "",
        module.head(10).to_markdown(index=False),
        "",
    ]
    return module, "\n".join(lines)
