"""Run core DevVCell response-recovery ablations without overwriting main outputs."""

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
from devvcell.response_recovery import classify_from_latent_tables, summarize_response_recovery  # noqa: E402
from devvcell.tables import latent_columns, read_table, write_table  # noqa: E402
from run_response_recovery_pipeline import ensure_lite_outputs, proxy_centroids, proxy_transferred  # noqa: E402


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


def classify_ablation(name: str, transferred: pd.DataFrame, centroids: pd.DataFrame, cfg: dict, out_dir: Path) -> pd.DataFrame:
    classes = classify_from_latent_tables(transferred, centroids, cfg)
    summary = summarize_response_recovery(classes)
    write_table(classes, out_dir / f"{name}_classes.csv")
    write_table(summary, out_dir / f"{name}_stage_summary.csv")
    return classes


def rows_for_mean_response(response_dictionary: pd.DataFrame, centroids: pd.DataFrame, mode: str) -> pd.DataFrame:
    response_cols = latent_columns(response_dictionary, prefixes=("response_latent_",))
    centroid_cols = latent_columns(centroids, prefixes=("latent_",))
    dims = min(len(response_cols), len(centroid_cols))
    response_cols = response_cols[:dims]

    if mode == "perturbation_mean":
        grouped = response_dictionary.groupby("perturbation", as_index=False, observed=True)[response_cols].mean()
    elif mode == "global_mean":
        perturbations = sorted(response_dictionary["perturbation"].astype(str).unique())
        mean_values = response_dictionary[response_cols].mean(axis=0)
        grouped = pd.DataFrame(
            [{"perturbation": perturbation, **{col: float(mean_values[col]) for col in response_cols}} for perturbation in perturbations]
        )
    else:
        raise ValueError(mode)

    rows: list[dict[str, object]] = []
    for _, response in grouped.iterrows():
        for _, centroid in centroids.iterrows():
            base = {
                "perturbation": response["perturbation"],
                "external_cell_type": mode,
                "stage": centroid["stage"],
                "stage_num": centroid.get("stage_num"),
                "cell_type": centroid["cell_type"],
                "devvcell_system": centroid.get("devvcell_system", "unknown"),
                "transfer_confidence": 1.0,
                "response_norm": float(np.linalg.norm(response[response_cols].astype(float).to_numpy())),
                "source_response_norm": float(np.linalg.norm(response[response_cols].astype(float).to_numpy())),
                "transfer_method": mode,
            }
            for idx, col in enumerate(response_cols):
                base[f"response_latent_{idx + 1:02d}"] = float(response[col])
            rows.append(base)
    return pd.DataFrame(rows)


def aggregate_classes(name: str, classes: pd.DataFrame, cfg: dict) -> dict[str, object]:
    window = cfg["stage_window_of_interest"]
    in_window = classes["stage_num"].astype(int).between(int(window["min_stage"]), int(window["max_stage"]))
    row: dict[str, object] = {
        "ablation": name,
        "n_cases": int(len(classes)),
        "window_cases": int(in_window.sum()),
        "outside_cases": int((~in_window).sum()),
        "mean_response_amplitude": float(classes["response_amplitude"].mean()),
        "mean_recovery_cost": float(classes["recovery_cost"].mean()),
        "mean_off_manifold_score": float(classes["off_manifold_score"].mean()),
    }
    for klass, count in classes["response_recovery_class"].value_counts().items():
        row[f"class_count_{klass}"] = int(count)
        row[f"class_fraction_{klass}"] = float(count / len(classes))
        row[f"window_fraction_{klass}"] = float(((classes["response_recovery_class"] == klass) & in_window).sum() / max(in_window.sum(), 1))
        row[f"outside_fraction_{klass}"] = float(((classes["response_recovery_class"] == klass) & ~in_window).sum() / max((~in_window).sum(), 1))
    return row


def main() -> None:
    args = parse_args()
    response_cfg = load_json(args.response_config)
    transfer_cfg = load_json(args.transfer_config)
    rng = np.random.default_rng(int(response_cfg["seed"]))
    out_dir = resolve_project_path("results/response_recovery/ablations")
    out_dir.mkdir(parents=True, exist_ok=True)

    centroids = read_table(existing_table(response_cfg["output"]["stage_celltype_centroids"]))
    transferred = read_table(existing_table(transfer_cfg["output"]["transferred_response_by_stage_celltype"]))
    dictionary = read_table(existing_table(transfer_cfg["output"]["external_response_dictionary"]))
    response_cols = latent_columns(transferred, prefixes=("response_latent_",))

    ablation_rows: list[dict[str, object]] = []

    main_classes = pd.read_csv(resolve_project_path(response_cfg["output"]["response_recovery_classes"]))
    write_table(main_classes, out_dir / "main_external_ot_classes.csv")
    ablation_rows.append(aggregate_classes("main_external_ot", main_classes, response_cfg))

    shuffled = transferred.copy()
    permutation = rng.permutation(len(shuffled))
    shuffled.loc[:, response_cols] = shuffled.loc[permutation, response_cols].to_numpy()
    shuffled["transfer_method"] = "shuffled_response_vectors"
    classes = classify_ablation("shuffled_response_vectors", shuffled, centroids, response_cfg, out_dir)
    ablation_rows.append(aggregate_classes("shuffled_response_vectors", classes, response_cfg))

    no_stage_transfer = rows_for_mean_response(dictionary, centroids, "perturbation_mean")
    classes = classify_ablation("no_stage_celltype_transfer", no_stage_transfer, centroids, response_cfg, out_dir)
    ablation_rows.append(aggregate_classes("no_stage_celltype_transfer", classes, response_cfg))

    global_mean = rows_for_mean_response(dictionary, centroids, "global_mean")
    classes = classify_ablation("global_mean_response", global_mean, centroids, response_cfg, out_dir)
    ablation_rows.append(aggregate_classes("global_mean_response", classes, response_cfg))

    stage_df, priority = ensure_lite_outputs(response_cfg)
    quick_centroids = proxy_centroids(stage_df)
    quick_transferred = proxy_transferred(stage_df, priority, response_cfg)
    classes = classify_ablation("no_external_rdeg_proxy", quick_transferred, quick_centroids, response_cfg, out_dir)
    ablation_rows.append(aggregate_classes("no_external_rdeg_proxy", classes, response_cfg))

    summary = pd.DataFrame(ablation_rows).fillna(0)
    output = write_table(summary, out_dir / "ablation_summary.csv")
    print(f"Wrote ablation summary: {output}")


if __name__ == "__main__":
    main()
