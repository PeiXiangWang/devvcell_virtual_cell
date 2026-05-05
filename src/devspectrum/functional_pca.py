"""Optional placeholder for future functional PCA extensions."""

from __future__ import annotations

import pandas as pd


def functional_pca_summary(timeseries: pd.DataFrame) -> pd.DataFrame:
    """Return an empty table in the MVP; reserved for future data-driven modes."""

    return pd.DataFrame(columns=["component", "explained_variance_ratio"])
