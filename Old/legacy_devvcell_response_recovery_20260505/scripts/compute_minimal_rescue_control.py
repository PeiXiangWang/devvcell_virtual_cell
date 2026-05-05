"""Compute a minimal rescue-control matrix from transferred response vectors."""

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
    parser.add_argument("--response-config", default="config/response_recovery.json")
    parser.add_argument("--transfer-config", default="config/perturbation_transfer.json")
    return parser.parse_args()


def existing_table(path_like: str) -> Path:
    path = resolve_project_path(path_like)
    if path.exists():
        return path
    if path.suffix == ".parquet" and path.with_suffix(".csv").exists():
        return path.with_suffix(".csv")
    raise FileNotFoundError(path)


def optimal_rescue(source: np.ndarray, candidate: np.ndarray) -> tuple[float, float, float]:
    denom = float(np.dot(candidate, candidate))
    if denom <= 1e-12:
        beta = 0.0
    else:
        beta = float(np.clip(-np.dot(source, candidate) / denom, 0.0, 1.0))
    residual = source + beta * candidate
    source_norm = float(np.linalg.norm(source))
    residual_norm = float(np.linalg.norm(residual))
    rescue_fraction = 0.0 if source_norm <= 1e-12 else max(0.0, 1.0 - residual_norm / source_norm)
    return beta, residual_norm, rescue_fraction


def main() -> None:
    args = parse_args()
    response_cfg = load_json(args.response_config)
    transfer_cfg = load_json(args.transfer_config)
    transferred = read_table(existing_table(transfer_cfg["output"]["transferred_response_by_stage_celltype"]))
    classes = pd.read_csv(resolve_project_path(response_cfg["output"]["response_recovery_classes"]))
    response_cols = latent_columns(transferred, prefixes=("response_latent_",))
    if not response_cols:
        raise ValueError("Transferred response table has no response_latent_* columns.")

    group_response_cols = ["stage", "cell_type", "perturbation"]
    numeric_aggs = {col: "mean" for col in response_cols}
    for col in ["stage_num", "transfer_confidence", "response_norm", "source_response_norm"]:
        if col in transferred.columns:
            numeric_aggs[col] = "mean"
    label_aggs = {
        col: "first"
        for col in ["devvcell_system", "transfer_method"]
        if col in transferred.columns
    }
    transferred = (
        transferred.groupby(group_response_cols, as_index=False, observed=True)
        .agg(**{col: (col, agg) for col, agg in {**numeric_aggs, **label_aggs}.items()})
    )

    key_cols = ["stage", "cell_type"]
    rows: list[dict[str, object]] = []
    for _, group in transferred.groupby(key_cols, observed=True, sort=False):
        if group["perturbation"].nunique() < 2:
            continue
        vectors = numeric_matrix(group, response_cols)
        records = group.reset_index(drop=True)
        for i, source_row in records.iterrows():
            source = vectors[i]
            best: dict[str, object] | None = None
            for j, candidate_row in records.iterrows():
                if i == j or str(source_row["perturbation"]) == str(candidate_row["perturbation"]):
                    continue
                beta, residual_norm, rescue_fraction = optimal_rescue(source, vectors[j])
                candidate = {
                    "stage": source_row["stage"],
                    "stage_num": source_row.get("stage_num"),
                    "cell_type": source_row["cell_type"],
                    "devvcell_system": source_row.get("devvcell_system"),
                    "perturbation": source_row["perturbation"],
                    "rescuer_perturbation": candidate_row["perturbation"],
                    "source_response_norm": float(np.linalg.norm(source)),
                    "rescuer_response_norm": float(np.linalg.norm(vectors[j])),
                    "optimal_rescuer_dose": beta,
                    "minimal_rescue_cost": residual_norm,
                    "rescue_fraction": rescue_fraction,
                    "source_transfer_confidence": source_row.get("transfer_confidence", np.nan),
                    "rescuer_transfer_confidence": candidate_row.get("transfer_confidence", np.nan),
                }
                if best is None or float(candidate["minimal_rescue_cost"]) < float(best["minimal_rescue_cost"]):
                    best = candidate
            if best is not None:
                rows.append(best)

    rescue = pd.DataFrame(rows)
    if not rescue.empty:
        class_cols = ["stage", "cell_type", "perturbation"]
        class_summary = (
            classes.groupby(class_cols, as_index=False, observed=True)
            .agg(
                response_recovery_class=("response_recovery_class", lambda values: values.mode().iloc[0]),
                recovery_cost=("recovery_cost", "mean"),
                off_manifold_score=("off_manifold_score", "mean"),
            )
        )
        rescue = rescue.merge(class_summary, on=["stage", "cell_type", "perturbation"], how="left")
        rescue = rescue.sort_values(["minimal_rescue_cost", "rescue_fraction"], ascending=[True, False])

    output = write_table(rescue, response_cfg["output"]["minimal_rescue_control_matrix"])
    print(f"Wrote minimal rescue control matrix: {output}")


if __name__ == "__main__":
    main()
