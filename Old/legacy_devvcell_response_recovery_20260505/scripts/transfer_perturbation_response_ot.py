"""Transfer external perturbation responses onto embryo stage/cell-type states."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from devvcell.io import load_json, resolve_project_path  # noqa: E402
from devvcell.tables import latent_columns, numeric_matrix, read_table, write_table  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/perturbation_transfer.json")
    return parser.parse_args()


def existing_table(path_like: str) -> Path:
    path = resolve_project_path(path_like)
    if path.exists():
        return path
    if path.suffix == ".parquet" and path.with_suffix(".csv").exists():
        return path.with_suffix(".csv")
    raise FileNotFoundError(path)


def soft_confidence(embryo: np.ndarray, controls: np.ndarray, quantile: float) -> np.ndarray:
    if len(controls) == 0:
        return np.ones(len(embryo), dtype=float)
    diff = embryo[:, None, :] - controls[None, :, :]
    dist = np.linalg.norm(diff, axis=2)
    min_dist = dist.min(axis=1)
    scale = float(np.quantile(min_dist, quantile))
    if not np.isfinite(scale) or scale <= 0:
        scale = float(np.nanmedian(min_dist[min_dist > 0])) if np.any(min_dist > 0) else 1.0
    return np.exp(-min_dist / max(scale, 1e-8))


def main() -> None:
    args = parse_args()
    cfg = load_json(args.config)
    response_path = existing_table(cfg["output"]["external_response_dictionary"])
    centroid_path = existing_table(cfg["input"]["stage_celltype_centroids"])

    responses = read_table(response_path)
    centroids = read_table(centroid_path)
    response_cols = latent_columns(responses, prefixes=("response_latent_",))
    control_cols = latent_columns(responses, prefixes=("control_latent_",))
    embryo_cols = latent_columns(centroids, prefixes=("latent_",))
    if not response_cols:
        raise ValueError("External response dictionary has no response_latent_* columns.")
    if not embryo_cols:
        raise ValueError("Embryo centroids have no latent_* columns.")

    dims = min(len(response_cols), len(embryo_cols))
    response_cols = response_cols[:dims]
    embryo_cols = embryo_cols[:dims]
    control_cols = control_cols[:dims] if control_cols else []

    embryo_values = numeric_matrix(centroids, embryo_cols)
    control_values = numeric_matrix(responses, control_cols) if control_cols else np.empty((0, dims))
    confidence = soft_confidence(
        embryo_values,
        control_values,
        float(cfg["ot_transfer"].get("confidence_distance_quantile", 0.90)),
    )

    rows: list[dict[str, object]] = []
    for _, response in responses.iterrows():
        vector = response[response_cols].astype(float).to_numpy()
        for idx, centroid in centroids.iterrows():
            scaled = vector * confidence[idx]
            base = {
                "perturbation": response.get("perturbation", "unknown"),
                "external_cell_type": response.get("external_cell_type", "all"),
                "stage": centroid["stage"],
                "stage_num": centroid.get("stage_num"),
                "cell_type": centroid["cell_type"],
                "devvcell_system": centroid.get("devvcell_system", "unknown"),
                "transfer_confidence": float(confidence[idx]),
                "response_norm": float(np.linalg.norm(scaled)),
                "source_response_norm": float(response.get("response_norm", np.linalg.norm(vector))),
                "transfer_method": cfg["ot_transfer"].get("method", "entropic_sinkhorn"),
            }
            base.update({f"response_latent_{i + 1:02d}": float(value) for i, value in enumerate(scaled)})
            rows.append(base)

    output = pd.DataFrame(rows)
    path = write_table(output, cfg["output"]["transferred_response_by_stage_celltype"])
    print(f"Wrote transferred responses: {path}")


if __name__ == "__main__":
    main()
