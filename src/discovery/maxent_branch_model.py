from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

from src.utils.config import ensure_dir, write_text


ROOT = Path(__file__).resolve().parents[2]


@dataclass
class MaxEntResult:
    fit: pd.DataFrame
    prediction: pd.DataFrame
    tier: str


def _unit(x: np.ndarray) -> np.ndarray:
    return x / np.maximum(np.linalg.norm(x, axis=1, keepdims=True), 1e-8)


def _neighbor_pairs(z: np.ndarray, k: int) -> np.ndarray:
    if z.shape[0] <= 2:
        return np.zeros((0, 2), dtype=int)
    k = min(k, z.shape[0] - 1)
    idx = NearestNeighbors(n_neighbors=k + 1).fit(z).kneighbors(z, return_distance=False)[:, 1:]
    src = np.repeat(np.arange(z.shape[0]), idx.shape[1])
    return np.column_stack([src, idx.ravel()])


def _pair_stats(z: np.ndarray, velocity: np.ndarray, labels: np.ndarray, k: int) -> dict:
    pairs = _neighbor_pairs(z, k)
    if pairs.size == 0:
        return {"pair_count": 0, "same_label_rate": 0.0, "velocity_alignment": 0.0, "random_same_label_rate": 0.0}
    vel = _unit(velocity)
    same = labels[pairs[:, 0]] == labels[pairs[:, 1]]
    align = np.sum(vel[pairs[:, 0]] * vel[pairs[:, 1]], axis=1)
    rng = np.random.default_rng(17)
    random_targets = rng.integers(0, z.shape[0], size=pairs.shape[0])
    random_same = labels[pairs[:, 0]] == labels[random_targets]
    return {
        "pair_count": int(pairs.shape[0]),
        "same_label_rate": float(np.mean(same)),
        "velocity_alignment": float(np.mean(align)),
        "random_same_label_rate": float(np.mean(random_same)),
    }


def fit_predict(
    dataset_frames: dict[str, dict],
    best_k: int,
    out_table_dir: str | Path = "tables",
    out_report_dir: str | Path = "reports",
    out_figure_dir: str | Path = "figures/discovery",
) -> MaxEntResult:
    """Fit a minimal pairwise topological model and predict branch order signs.

    The model is deliberately small:

    E = -J sum_(i,j in kNN) 1[state_i = state_j]
        -h sum_(i,j in kNN) cos(v_i, v_j)

    J and h are estimated from local pairwise statistics only. The model is not
    trained on branch event labels.
    """
    fit_rows: list[dict] = []
    pred_rows: list[dict] = []
    for dataset, payload in dataset_frames.items():
        z = payload["z"]
        v = payload["velocity"]
        labels = payload["lineage"].astype(str).to_numpy()
        order_effect = payload.get("lineage_separation_effect", np.nan)
        align_effect = payload.get("alignment_effect", np.nan)
        stats = _pair_stats(z, v, labels, best_k)
        same = min(max(stats["same_label_rate"], 1e-5), 1 - 1e-5)
        random_same = min(max(stats["random_same_label_rate"], 1e-5), 1 - 1e-5)
        j_hat = float(np.log(same / (1 - same)) - np.log(random_same / (1 - random_same)))
        h_hat = float(stats["velocity_alignment"])
        predicted_alignment_effect = float(np.tanh(h_hat + 0.25 * j_hat) * 0.05)
        predicted_condensation = bool(j_hat > 0 and h_hat > 0)
        predicted_separation_effect = float(-abs(j_hat) / (1.0 + abs(j_hat)) if predicted_condensation else abs(j_hat) / (1.0 + abs(j_hat)))
        fit_rows.append(
            {
                "dataset": dataset,
                "k": best_k,
                "J_pairwise_label": j_hat,
                "h_velocity_alignment": h_hat,
                "pair_count": stats["pair_count"],
                "same_label_rate": stats["same_label_rate"],
                "random_same_label_rate": stats["random_same_label_rate"],
            }
        )
        pred_rows.append(
            {
                "dataset": dataset,
                "k": best_k,
                "observed_lineage_separation_effect": order_effect,
                "predicted_lineage_separation_effect": predicted_separation_effect,
                "observed_alignment_effect": align_effect,
                "predicted_alignment_effect": predicted_alignment_effect,
                "condensation_direction_match": bool(np.sign(order_effect) == np.sign(predicted_separation_effect)) if np.isfinite(order_effect) else False,
                "alignment_direction_match": bool(np.sign(align_effect) == np.sign(predicted_alignment_effect)) if np.isfinite(align_effect) else False,
            }
        )
    fit = pd.DataFrame(fit_rows)
    pred = pd.DataFrame(pred_rows)
    internal_match = pred[pred["dataset"].eq("internal")]["condensation_direction_match"].any()
    e1_match = pred[pred["dataset"].eq("E1")]["condensation_direction_match"].any()
    if internal_match and e1_match:
        tier = "acceptable"
    elif internal_match or e1_match:
        tier = "weak"
    else:
        tier = "fail"
    out_table_dir = ensure_dir(out_table_dir)
    out_report_dir = ensure_dir(out_report_dir)
    out_figure_dir = ensure_dir(out_figure_dir)
    fit.to_csv(Path(out_table_dir) / "maxent_model_fit.csv", index=False)
    pred.to_csv(Path(out_table_dir) / "maxent_branch_prediction.csv", index=False)
    fig, ax = plt.subplots(figsize=(6, 4))
    if not pred.empty:
        ax.scatter(pred["observed_lineage_separation_effect"], pred["predicted_lineage_separation_effect"], s=70)
        for row in pred.itertuples(index=False):
            ax.annotate(row.dataset, (row.observed_lineage_separation_effect, row.predicted_lineage_separation_effect), fontsize=8)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("observed separation effect")
    ax.set_ylabel("predicted separation effect")
    ax.set_title("Pairwise topological minimal model")
    fig.tight_layout()
    fig.savefig(Path(out_figure_dir) / "maxent_predicted_vs_observed_order_parameters.png", dpi=180)
    main_dir = ensure_dir("figures/main")
    fig.savefig(Path(main_dir) / "figure9_maxent_minimal_model.png", dpi=180)
    plt.close(fig)
    write_text(
        Path(out_report_dir) / "maxent_minimal_model_report.md",
        "# Maximum-Entropy Minimal Model\n\n"
        f"- maxent_model_tier: {tier}\n"
        f"- selected_k: {best_k}\n\n"
        "The model uses only pairwise topological label similarity and velocity alignment. It is not trained on branch-event labels. Because it is a prototype, a weak or failed result is interpreted as limited support for the minimal-model explanation rather than as a biological negative result.\n\n"
        "## Fit\n\n"
        + fit.to_markdown(index=False)
        + "\n\n## Prediction\n\n"
        + pred.to_markdown(index=False)
        + "\n",
    )
    return MaxEntResult(fit=fit, prediction=pred, tier=tier)
