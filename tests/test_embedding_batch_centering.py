import numpy as np
import pandas as pd

from devguard.embedding import apply_batch_centering, fit_batch_centering


def test_batch_centering_uses_fallback_batch_column():
    embeddings = np.array([[1.0, 2.0], [3.0, 4.0], [11.0, 12.0], [13.0, 14.0]])
    obs = pd.DataFrame({"source_dataset_id": ["a", "a", "b", "b"]})
    state = fit_batch_centering(embeddings, obs, column="source_dataset_id")

    query = np.array([[2.0, 3.0], [12.0, 13.0], [100.0, 100.0]])
    query_obs = pd.DataFrame({"dataset_id": ["a", "b", "unknown"]})
    centered = apply_batch_centering(query, query_obs, state, fallback_columns=["dataset_id"])

    np.testing.assert_allclose(centered[0], [0.0, 0.0])
    np.testing.assert_allclose(centered[1], [0.0, 0.0])
    np.testing.assert_allclose(centered[2], [93.0, 92.0])
