from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch

from src.model.dynamics import DynamicsFlags, SwarmLineageDynamics
from src.model.simulator import feature_matrix
from src.model.swarm_rules import local_density
from src.utils.config import ensure_dir, load_config


TIER_ORDER = {"fail": 0, "weak": 1, "acceptable": 2, "strong": 3}
INTERPRETATION_LEVELS = {
    "retained_computational_hypothesis",
    "encoded_control_law_recovery",
    "rollout_supported_mechanistic_probe",
    "exploratory_sensitivity",
    "demonstration_only",
    "unsupported",
}


def tier_at_least(tier: str, minimum: str) -> bool:
    return TIER_ORDER.get(str(tier), 0) >= TIER_ORDER.get(str(minimum), 0)


def configure_paths(config_path: str, quick_fixture: bool = False) -> tuple[dict, dict, dict]:
    train_cfg = load_config(config_path)
    model_cfg = load_config(train_cfg.get("model_config", "configs/model.yaml"))
    discovery_cfg = load_config(train_cfg.get("discovery_config", "configs/discovery.yaml"))
    if quick_fixture:
        model_cfg = dict(model_cfg)
        train_cfg = dict(train_cfg)
        model_cfg["teacher_path"] = "processed/quick_fixture/ot_teacher.h5ad"
        model_cfg["model_dir"] = "results/quick_fixture/models"
        model_cfg["event_log_path"] = "tables/quick_fixture/birth_death_event_log.csv"
        model_cfg["order_log_path"] = "tables/quick_fixture/rollout_order_parameters.csv"
        model_cfg["metrics_path"] = "tables/quick_fixture/final_metrics.csv"
        train_cfg["discovery_prefix"] = "quick_fixture"
    return train_cfg, model_cfg, discovery_cfg


def output_dirs(train_cfg: dict, discovery_cfg: dict | None = None) -> tuple[Path, Path, Path]:
    discovery_cfg = discovery_cfg or {}
    out_cfg = discovery_cfg.get("output", {})
    prefix = train_cfg.get("discovery_prefix")
    table_dir = Path(out_cfg.get("quick_tables_dir", "tables/quick_fixture")) if prefix else Path(out_cfg.get("tables_dir", "tables"))
    report_dir = Path(out_cfg.get("quick_reports_dir", "reports/quick_fixture")) if prefix else Path(out_cfg.get("reports_dir", "reports"))
    fig_dir = Path(out_cfg.get("quick_figures_dir", "figures/quick_fixture/discovery")) if prefix else Path(out_cfg.get("figures_dir", "figures/discovery"))
    ensure_dir(table_dir)
    ensure_dir(report_dir)
    ensure_dir(fig_dir)
    return table_dir, report_dir, fig_dir


def seed_list(model_cfg: dict, quick_fixture: bool = False) -> list[int]:
    seeds = [int(s) for s in model_cfg.get("seeds", [7, 17, 23, 42, 99])]
    return seeds[:2] if quick_fixture else seeds


def load_teacher(model_cfg: dict) -> ad.AnnData:
    return ad.read_h5ad(model_cfg["teacher_path"])


def load_teacher_model(adata: ad.AnnData, model_cfg: dict, seed: int = 7) -> SwarmLineageDynamics | None:
    model_path = Path(model_cfg.get("model_dir", "results/swarmlineage/models")) / f"dynamics_seed{seed}.pt"
    if not model_path.exists():
        return None
    x = feature_matrix(adata, np.arange(min(adata.n_obs, 8)), model_cfg)
    latent_dim = adata.obsm[model_cfg.get("latent_key", "X_pca")].shape[1]
    model = SwarmLineageDynamics(x.shape[1], latent_dim, int(model_cfg.get("hidden_dim", 96)))
    try:
        state = torch.load(model_path, map_location="cpu", weights_only=True)
    except TypeError:
        state = torch.load(model_path, map_location="cpu")
    model.load_state_dict(state["teacher_state_dict"])
    model.eval()
    return model


def _fate_entropy(fate: np.ndarray) -> np.ndarray:
    if fate.size == 0:
        return np.zeros(fate.shape[0], dtype=float)
    p = np.clip(fate, 1e-12, 1.0)
    return -np.sum(p * np.log(p), axis=1) / max(np.log(p.shape[1]), 1.0)


def cell_feature_frame(adata: ad.AnnData, model_cfg: dict, seed: int = 7) -> pd.DataFrame:
    idx = np.arange(adata.n_obs)
    z = np.asarray(adata.obsm[model_cfg.get("latent_key", "X_pca")], dtype=float)
    obs = adata.obs.copy()
    entropy = pd.to_numeric(obs.get("ot_transition_entropy", 0.5), errors="coerce").fillna(0.5).to_numpy(dtype=float)
    growth = pd.to_numeric(obs.get("ot_growth", 1.0), errors="coerce").fillna(1.0).to_numpy(dtype=float)
    cycle = pd.to_numeric(obs.get("cell_cycle_score", 0.0), errors="coerce").fillna(0.0).to_numpy(dtype=float)
    cci = pd.to_numeric(obs.get("cci_signal", 0.0), errors="coerce").fillna(0.0).to_numpy(dtype=float)
    fate_cols = [c for c in obs.columns if c.startswith("fate_prob_")]
    fate = obs[fate_cols].to_numpy(dtype=float) if fate_cols else np.zeros((adata.n_obs, 1))
    fate_max = fate.max(axis=1) if fate.size else np.zeros(adata.n_obs)
    fate_entropy = _fate_entropy(fate)
    density = local_density(z)
    velocity = np.asarray(adata.obsm.get("X_ot_velocity", np.zeros_like(z)), dtype=float)
    displacement = np.linalg.norm(velocity, axis=1)
    time_key = model_cfg.get("time_key", "time_numeric")
    frame = pd.DataFrame(
        {
            "cell_id": obs.index.astype(str),
            "seed": int(seed),
            "time_numeric": pd.to_numeric(obs.get(time_key), errors="coerce").to_numpy(dtype=float),
            "lineage": obs.get(model_cfg.get("cell_type_key", "lineage"), "unknown").astype(str).to_numpy(),
            "ot_transition_entropy": entropy,
            "local_density": density,
            "fate_probability_max": fate_max,
            "fate_entropy": fate_entropy,
            "cell_cycle_score": cycle,
            "cci_signal": cci,
            "ot_growth": growth,
            "ot_displacement": displacement,
        }
    )
    model = load_teacher_model(adata, model_cfg, seed=seed)
    if model is not None:
        x = torch.as_tensor(feature_matrix(adata, idx, model_cfg), dtype=torch.float32)
        flags = DynamicsFlags(use_teacher=True, use_birth_death=True, use_diffusion=True, use_cci=True, use_swarm=True, use_memory=True)
        with torch.no_grad():
            frame["learned_sigma"] = model.sigma(x, flags).cpu().numpy()
            frame["birth_hazard"] = model.birth_hazard(x, flags).cpu().numpy()
            frame["death_hazard"] = model.death_hazard(x, flags).cpu().numpy()
            frame["net_growth_hazard"] = frame["birth_hazard"] - frame["death_hazard"]
            frame["learned_displacement"] = np.linalg.norm(model.vector_field(x, flags=flags).cpu().numpy(), axis=1)
    else:
        frame["learned_sigma"] = 0.015 + 0.12 * frame["ot_transition_entropy"]
        frame["birth_hazard"] = np.maximum(np.log(np.maximum(frame["ot_growth"], 1e-3)), 0.0)
        frame["death_hazard"] = np.maximum(-np.log(np.maximum(frame["ot_growth"], 1e-3)), 0.0)
        frame["net_growth_hazard"] = frame["birth_hazard"] - frame["death_hazard"]
        frame["learned_displacement"] = frame["ot_displacement"]
    return frame


def seedwise_feature_frame(adata: ad.AnnData, model_cfg: dict, seeds: list[int]) -> pd.DataFrame:
    return pd.concat([cell_feature_frame(adata, model_cfg, seed=s) for s in seeds], ignore_index=True)


def linear_effects(frame: pd.DataFrame, response: str, predictors: list[str]) -> pd.DataFrame:
    rows = []
    data = frame[[response] + predictors].replace([np.inf, -np.inf], np.nan).dropna()
    if data.shape[0] < len(predictors) + 5:
        return pd.DataFrame(columns=["response", "predictor", "coef", "abs_coef", "r2", "n"])
    y = data[response].to_numpy(dtype=float)
    x = data[predictors].to_numpy(dtype=float)
    x = (x - x.mean(axis=0)) / np.maximum(x.std(axis=0), 1e-8)
    y_center = y - y.mean()
    beta, *_ = np.linalg.lstsq(np.c_[np.ones(x.shape[0]), x], y_center, rcond=None)
    pred = np.c_[np.ones(x.shape[0]), x] @ beta + y.mean()
    r2 = 1.0 - float(np.sum((y - pred) ** 2) / max(np.sum((y - y.mean()) ** 2), 1e-12))
    for p, c in zip(predictors, beta[1:]):
        rows.append({"response": response, "predictor": p, "coef": float(c), "abs_coef": float(abs(c)), "r2": r2, "n": int(data.shape[0])})
    return pd.DataFrame(rows)


def standardized_coef(frame: pd.DataFrame, response: str, predictor: str, covariates: list[str] | None = None) -> float:
    covariates = covariates or []
    stats = linear_effects(frame, response, [predictor] + covariates)
    hit = stats[stats["predictor"] == predictor]
    return float(hit["coef"].iloc[0]) if not hit.empty else float("nan")


def bootstrap_ci(values: list[float] | np.ndarray, seed: int = 7, repeats: int = 500) -> tuple[float, float, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan"), float("nan"), float("nan")
    if arr.size == 1:
        v = float(arr[0])
        return v, v, v
    rng = np.random.default_rng(seed)
    boots = np.array([rng.choice(arr, size=arr.size, replace=True).mean() for _ in range(repeats)])
    return float(arr.mean()), float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))


def permutation_p_value(observed: float, null_values: list[float] | np.ndarray) -> float:
    null = np.asarray(null_values, dtype=float)
    null = null[np.isfinite(null)]
    if not np.isfinite(observed) or null.size == 0:
        return 1.0
    return float((1 + np.sum(np.abs(null) >= abs(observed))) / (null.size + 1))


def bh_q_values(p_values: list[float]) -> list[float]:
    p = np.asarray([1.0 if not np.isfinite(v) else float(v) for v in p_values], dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        idx = order[i]
        val = min(prev, p[idx] * n / (i + 1))
        ranked[idx] = val
        prev = val
    return ranked.tolist()


def seed_stability(effects: list[float] | np.ndarray, min_seed_count: int) -> tuple[bool, float]:
    arr = np.asarray(effects, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size < min_seed_count or arr.size == 0:
        return False, 0.0
    signs = np.sign(arr)
    majority = max(float((signs >= 0).mean()), float((signs <= 0).mean()))
    return bool(majority >= 0.8 and abs(arr.mean()) > 0), majority


def law_tier(
    effect_size: float,
    permutation_q: float,
    negative_control_pass: bool,
    seed_stability_pass: bool,
    rollout_based: bool,
    directly_supervised_or_encoded: bool,
    discovery_cfg: dict,
) -> str:
    cfg = discovery_cfg.get("emergent_law", {})
    acceptable_effect = float(cfg.get("min_effect_size_for_acceptable", 0.01))
    strong_effect = float(cfg.get("min_effect_size_for_strong", 0.03))
    acceptable_q = float(cfg.get("max_permutation_q_for_acceptable", 0.10))
    strong_q = float(cfg.get("max_permutation_q_for_strong", 0.05))
    effect = abs(float(effect_size)) if np.isfinite(effect_size) else 0.0
    if (
        effect >= strong_effect
        and permutation_q <= strong_q
        and negative_control_pass
        and seed_stability_pass
        and rollout_based
        and not directly_supervised_or_encoded
    ):
        return "strong"
    if effect >= acceptable_effect and seed_stability_pass and (permutation_q <= acceptable_q or negative_control_pass):
        return "acceptable"
    if effect >= acceptable_effect * 0.5:
        return "weak"
    return "fail"


def interpretation_level(tier: str, rollout_based: bool, directly_encoded: bool, negative_control_pass: bool) -> str:
    if tier == "fail":
        return "unsupported"
    if directly_encoded:
        return "encoded_control_law_recovery"
    if rollout_based and negative_control_pass and tier in {"acceptable", "strong"}:
        return "rollout_supported_mechanistic_probe" if tier == "acceptable" else "retained_computational_hypothesis"
    if tier == "weak":
        return "exploratory_sensitivity"
    return "demonstration_only"


def gate_record(
    law: str,
    tier: str,
    effect_size: float,
    ci: tuple[float, float, float],
    permutation_p: float,
    permutation_q: float,
    negative_control_pass: bool,
    seed_stability_pass: bool,
    rollout_based: bool,
    directly_supervised_or_encoded: bool,
    table: str,
    report: str,
    status: str = "executed",
) -> dict:
    level = interpretation_level(tier, rollout_based, directly_supervised_or_encoded, negative_control_pass)
    if level not in INTERPRETATION_LEVELS:
        level = "unsupported"
    return {
        "law": law,
        "tier": tier,
        "gate_pass": tier in {"acceptable", "strong"},
        "strong_gate": tier == "strong",
        "effect_size": float(effect_size) if np.isfinite(effect_size) else np.nan,
        "effect_ci_low": float(ci[1]) if np.isfinite(ci[1]) else np.nan,
        "effect_ci_high": float(ci[2]) if np.isfinite(ci[2]) else np.nan,
        "permutation_p": float(permutation_p) if np.isfinite(permutation_p) else 1.0,
        "permutation_q": float(permutation_q) if np.isfinite(permutation_q) else 1.0,
        "negative_control_pass": bool(negative_control_pass),
        "seed_stability_pass": bool(seed_stability_pass),
        "rollout_based": bool(rollout_based),
        "directly_supervised_or_encoded": bool(directly_supervised_or_encoded),
        "interpretation_level": level,
        "table": table,
        "report": report,
        "status": status,
    }


def write_report(path: Path, title: str, body: list[str]) -> None:
    from src.utils.config import write_text

    write_text(path, "\n".join([f"# {title}", "", *body, ""]))
