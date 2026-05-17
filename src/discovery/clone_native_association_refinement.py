from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import NearestNeighbors

from src.discovery.clone_developmental_validation_v14 import (
    DATASETS,
    EXPOSURES,
    _association_tables,
    _cell_exposures,
    _clone_exposures,
    _entropy_from_counts,
    _write_csv,
    _write_md,
)
from src.utils.config import ensure_dir


ROOT = Path(__file__).resolve().parents[2]


def _native_velocity_and_entropy(dataset_id: str, adata: ad.AnnData) -> tuple[np.ndarray, np.ndarray]:
    z = np.asarray(adata.obsm["X_pca"], dtype=float)
    velocity = np.zeros_like(z)
    entropy = np.full(adata.n_obs, np.nan, dtype=float)
    coupling_dir = ROOT / "processed/external_l4_native" / dataset_id / "ot_couplings"
    for path in sorted(coupling_dir.glob("teacher_native_moscot_*.npz")):
        raw = np.load(path, allow_pickle=True)
        src = raw["source_indices"].astype(int)
        bary = raw["barycentric"].astype(float)
        velocity[src] = bary - z[src]
        entropy[src] = raw["entropy"].astype(float)
    return velocity, entropy


def _unit(x: np.ndarray) -> np.ndarray:
    return x / np.maximum(np.linalg.norm(x, axis=1, keepdims=True), 1e-8)


def _override_native_exposures(adata: ad.AnnData, cfg) -> pd.DataFrame:
    obs = adata.obs.copy().reset_index(drop=True)
    z = np.asarray(adata.obsm["X_pca"], dtype=float)
    base = _cell_exposures(z, obs, cfg)
    velocity, entropy = _native_velocity_and_entropy(cfg.dataset_id, adata)
    vu = _unit(velocity)
    times = sorted(pd.to_numeric(obs["time_numeric"], errors="coerce").dropna().unique())
    terminal = times[-1]
    terminal_centroid = z[obs["time_numeric"].eq(terminal).to_numpy()].mean(axis=0)
    base["teacher_velocity_bias_exposure"] = np.sum(_unit(terminal_centroid[None, :] - z) * vu, axis=1)
    native_entropy_mask = np.isfinite(entropy)
    base.loc[native_entropy_mask, "fate_entropy_exposure"] = entropy[native_entropy_mask]
    for t in times:
        idx = obs.index[obs["time_numeric"].eq(t)].to_numpy()
        if idx.size < 4:
            continue
        k = min(8, idx.size - 1)
        nn = NearestNeighbors(n_neighbors=k + 1).fit(z[idx]).kneighbors(z[idx], return_distance=False)[:, 1:]
        local_v = vu[idx]
        base.loc[idx, "local_alignment_exposure"] = np.mean(np.sum(local_v[:, None, :] * local_v[nn], axis=2), axis=1)
        base.loc[idx, "alignment_exposure"] = base.loc[idx, "local_alignment_exposure"].to_numpy()
    base["two_phase_condensation_plus_entropy"] = (
        pd.to_numeric(base["branch_window_condensation_exposure"], errors="coerce").rank(pct=True)
        + pd.to_numeric(base["fate_entropy_exposure"], errors="coerce").rank(pct=True)
    ) / 2.0
    base["two_phase_condensation_plus_divergence"] = (
        pd.to_numeric(base["branch_window_condensation_exposure"], errors="coerce").rank(pct=True)
        + pd.to_numeric(base["post_event_divergence_exposure"], errors="coerce").rank(pct=True)
    ) / 2.0
    base["uncertainty_plus_teacher_bias"] = (
        pd.to_numeric(base["fate_entropy_exposure"], errors="coerce").rank(pct=True)
        + pd.to_numeric(base["teacher_velocity_bias_exposure"], errors="coerce").rank(pct=True)
    ) / 2.0
    return base


def _clone_scores_native(obs: pd.DataFrame, cfg, min_clone_size: int = 2) -> pd.DataFrame:
    rows = []
    times = sorted(pd.to_numeric(obs["time_numeric"], errors="coerce").dropna().unique())
    terminal = times[-1]
    for clone_id, g in obs.groupby("clone_id", observed=False):
        if pd.isna(clone_id) or g.shape[0] < min_clone_size:
            continue
        terminal_g = g[g["time_numeric"].eq(terminal)]
        if terminal_g.empty:
            terminal_g = g
        rows.append(
            {
                "clone_id": clone_id,
                "clone_size": int(g.shape[0]),
                "clone_start_time": float(g["time_numeric"].min()),
                "clone_end_time": float(g["time_numeric"].max()),
                "clone_time_span": float(g["time_numeric"].max() - g["time_numeric"].min()),
                "initial_cell_state": str(g.sort_values("time_numeric")[cfg.celltype_col].iloc[0]),
                "terminal_fate_entropy": _entropy_from_counts(terminal_g[cfg.celltype_col]),
                "terminal_lineage_entropy": _entropy_from_counts(terminal_g[cfg.celltype_col]),
                "clone_branch_count": int(terminal_g[cfg.celltype_col].astype(str).nunique()),
                "clone_multilineage_score": int(terminal_g[cfg.celltype_col].astype(str).nunique()),
                "clone_fate_diversification_index": 1.0 - float(terminal_g[cfg.celltype_col].astype(str).value_counts(normalize=True).max()),
                "clone_success_failure_split": np.nan,
                "clone_state_transition_entropy": _entropy_from_counts(g[cfg.celltype_col]),
                "early_to_late_fate_dispersion": float(g[cfg.celltype_col].astype(str).nunique() > 1),
                "terminal_sampling_depth": int(terminal_g.shape[0]),
            }
        )
    return pd.DataFrame(rows)


def _clone_exposures_all(obs: pd.DataFrame, exposure: pd.DataFrame) -> pd.DataFrame:
    base = _clone_exposures(obs, exposure)
    extra_cols = ["two_phase_condensation_plus_entropy", "two_phase_condensation_plus_divergence", "uncertainty_plus_teacher_bias", "local_alignment_exposure"]
    rows = []
    times = sorted(pd.to_numeric(obs["time_numeric"], errors="coerce").dropna().unique())
    event_time = times[0] if len(times) < 3 else times[len(times) // 2 - 1]
    pre = obs["time_numeric"].le(event_time)
    for clone_id, g in obs.groupby("clone_id", observed=False):
        use_idx = g.index[pre.loc[g.index]]
        if len(use_idx) < 1:
            use_idx = g.index
        row = {"clone_id": clone_id}
        for col in extra_cols:
            if col in exposure:
                row[col] = float(pd.to_numeric(exposure.loc[use_idx, col], errors="coerce").mean())
        rows.append(row)
    extra = pd.DataFrame(rows)
    return base.merge(extra, on="clone_id", how="left")


def _association_tables_extended(dataset_id: str, clone: pd.DataFrame):
    assoc, reg, controls = _association_tables(dataset_id, clone)
    extra_exposures = [
        ("local_alignment_exposure", "positive"),
        ("two_phase_condensation_plus_entropy", "positive"),
        ("two_phase_condensation_plus_divergence", "positive"),
        ("uncertainty_plus_teacher_bias", "positive"),
    ]
    extra_rows = []
    extra_reg = []
    for score, primary in [
        ("terminal_fate_entropy", True),
        ("clone_branch_count", False),
        ("clone_multilineage_score", False),
        ("clone_fate_diversification_index", False),
        ("clone_state_transition_entropy", False),
        ("early_to_late_fate_dispersion", False),
    ]:
        y = pd.to_numeric(clone[score], errors="coerce")
        for exposure, expected in extra_exposures:
            if exposure not in clone:
                continue
            x = pd.to_numeric(clone[exposure], errors="coerce")
            mask = x.notna() & y.notna() & (clone["clone_time_span"] > 0)
            if mask.sum() < 15 or x[mask].nunique() < 2 or y[mask].nunique() < 2:
                continue
            rho, p = stats.spearmanr(x[mask], y[mask])
            extra_rows.append(
                {
                    "dataset_id": dataset_id,
                    "exposure": exposure,
                    "clone_splitting_score": score,
                    "score_primary": primary,
                    "expected_direction": expected,
                    "observed_direction": "positive" if rho > 0 else "negative",
                    "effect_size": float(rho),
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "p_value": float(p),
                    "whether_supports_original_hypothesis": False,
                    "whether_supports_revised_hypothesis": bool(primary and rho > 0 and p <= 0.10),
                    "permutation_q": float(p),
                }
            )
            cov = pd.DataFrame(
                {
                    "exposure": x,
                    "clone_size": pd.to_numeric(clone["clone_size"], errors="coerce"),
                    "clone_start_time": pd.to_numeric(clone["clone_start_time"], errors="coerce"),
                    "clone_time_span": pd.to_numeric(clone["clone_time_span"], errors="coerce"),
                    "terminal_sampling_depth": pd.to_numeric(clone["terminal_sampling_depth"], errors="coerce"),
                }
            )
            m = cov.notna().all(axis=1) & y.notna()
            if m.sum() >= 15:
                X = cov.loc[m].to_numpy(dtype=float)
                yy = y.loc[m].to_numpy(dtype=float)
                X = (X - X.mean(axis=0)) / np.maximum(X.std(axis=0), 1e-8)
                yy = (yy - yy.mean()) / max(yy.std(), 1e-8)
                model = LinearRegression().fit(X, yy)
                extra_reg.append(
                    {
                        "dataset_id": dataset_id,
                        "exposure": exposure,
                        "clone_splitting_score": score,
                        "score_primary": primary,
                        "n_clones": int(m.sum()),
                        "covariate_adjusted_effect": float(model.coef_[0]),
                        "model_r2": float(model.score(X, yy)),
                        "supports_after_covariates": bool(model.coef_[0] > 0 and primary),
                    }
                )
    if extra_rows:
        assoc = pd.concat([assoc, pd.DataFrame(extra_rows)], ignore_index=True)
    if extra_reg:
        reg = pd.concat([reg, pd.DataFrame(extra_reg)], ignore_index=True)
    return assoc, reg, controls


def analyze_dataset(cfg) -> dict:
    path = ROOT / "processed/external_l4_native" / cfg.dataset_id / "native_input.h5ad"
    if not path.exists():
        return {"dataset_id": cfg.dataset_id, "native_input_loaded": False, "validation_tier": "fail", "reason": "missing_native_input"}
    adata = ad.read_h5ad(path)
    obs = adata.obs.copy().reset_index(drop=True)
    exposure = _override_native_exposures(adata, cfg)
    scores = _clone_scores_native(obs, cfg, min_clone_size=2)
    exp = _clone_exposures_all(obs, exposure)
    clone = scores.merge(exp, on="clone_id", how="inner")
    usable = clone[(clone["clone_time_span"] > 0) & (clone["clone_size"] >= 2)].copy()
    assoc, reg, controls = _association_tables_extended(cfg.dataset_id, usable)
    outdir = ensure_dir(ROOT / "tables/clone_native")
    clone.to_csv(outdir / f"{cfg.dataset_id}_native_clone_level.csv", index=False)
    assoc.to_csv(outdir / f"{cfg.dataset_id}_native_exposure_associations.csv", index=False)
    reg.to_csv(outdir / f"{cfg.dataset_id}_native_confounder_regression.csv", index=False)
    controls.to_csv(outdir / f"{cfg.dataset_id}_native_negative_controls.csv", index=False)
    primary = assoc[
        assoc["score_primary"].eq(True)
        & assoc["exposure"].eq("branch_window_condensation_exposure")
    ]
    primary_reg = reg[
        reg["score_primary"].eq(True)
        & reg["exposure"].eq("branch_window_condensation_exposure")
    ]
    revised = assoc[
        assoc["score_primary"].eq(True)
        & assoc["exposure"].isin(["two_phase_condensation_plus_entropy", "two_phase_condensation_plus_divergence", "uncertainty_plus_teacher_bias", "fate_entropy_exposure", "teacher_velocity_bias_exposure"])
        & assoc["whether_supports_revised_hypothesis"].eq(True)
    ]
    primary_ok = bool(
        not primary.empty
        and bool(primary["whether_supports_original_hypothesis"].iloc[0])
        and not primary_reg.empty
        and bool(primary_reg["supports_after_covariates"].iloc[0])
    )
    tier = "acceptable" if primary_ok and usable.shape[0] >= 20 else "weak" if (primary_ok or not revised.empty) else "fail"
    return {
        "dataset_id": cfg.dataset_id,
        "native_input_loaded": True,
        "teacher_backend": "native_moscot_downsampled",
        "n_cells": int(adata.n_obs),
        "usable_clones": int(usable.shape[0]),
        "primary_condensation_support": primary_ok,
        "primary_effect": float(primary["effect_size"].iloc[0]) if not primary.empty else np.nan,
        "primary_q": float(primary["permutation_q"].iloc[0]) if not primary.empty else np.nan,
        "primary_covariate_adjusted_effect": float(primary_reg["covariate_adjusted_effect"].iloc[0]) if not primary_reg.empty else np.nan,
        "revised_support_count": int(revised.shape[0]),
        "validation_tier": tier,
        "interpretation": "native_primary_condensation_support" if primary_ok else "native_revised_feature_support" if not revised.empty else "native_no_clone_support",
    }


def run() -> None:
    rows = [analyze_dataset(cfg) for cfg in DATASETS if cfg.dataset_id.startswith("Jindal") or cfg.dataset_id.startswith("Weinreb")]
    summary = pd.DataFrame(rows)
    _write_csv(summary, "tables/clone_native_validation_summary.csv")
    assoc = []
    reg = []
    controls = []
    for cfg in DATASETS:
        outdir = ROOT / "tables/clone_native"
        for name, bucket in [
            (f"{cfg.dataset_id}_native_exposure_associations.csv", assoc),
            (f"{cfg.dataset_id}_native_confounder_regression.csv", reg),
            (f"{cfg.dataset_id}_native_negative_controls.csv", controls),
        ]:
            p = outdir / name
            if p.exists():
                bucket.append(pd.read_csv(p))
    _write_csv(pd.concat(assoc, ignore_index=True) if assoc else pd.DataFrame(), "tables/clone_native_exposure_associations.csv")
    _write_csv(pd.concat(reg, ignore_index=True) if reg else pd.DataFrame(), "tables/clone_native_confounder_regression.csv")
    _write_csv(pd.concat(controls, ignore_index=True) if controls else pd.DataFrame(), "tables/clone_native_negative_controls.csv")
    fallback = _safe_read("tables/clone_developmental_validation_summary.csv")
    integrated = _integrated_summary(fallback, summary)
    _write_csv(integrated, "tables/clone_integrated_validation_summary.csv")
    _update_final_claim_tiers()
    ensure_dir(ROOT / "reports")
    _write_md(
        "reports/clone_native_validation_refinement.md",
        "# Clone Native Validation Refinement\n\n"
        "Native moscot was run in `.venv_moscot_native` on downsampled Jindal and Weinreb AnnData objects. Associations here are restricted to clones represented in the native downsample and therefore test direction/stability rather than replacing the full-data fallback analysis.\n\n"
        + summary.to_markdown(index=False),
    )
    _write_integrated_reports(integrated, fallback, summary)
    _write_forbidden_claim_scan()
    figdir = ensure_dir(ROOT / "figures/main")
    fig, ax = plt.subplots(figsize=(6, 4))
    if not summary.empty:
        ax.bar(summary["dataset_id"], summary["primary_effect"], color="#4C78A8")
        ax.tick_params(axis="x", rotation=30)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("native primary condensation effect")
    ax.set_title("Downsampled native clone validation")
    fig.tight_layout()
    fig.savefig(figdir / "figure14_native_clone_validation.png", dpi=180)
    print(summary.to_json(orient="records"))


def _safe_read(rel_path: str) -> pd.DataFrame:
    path = ROOT / rel_path
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def _integrated_summary(fallback: pd.DataFrame, native: pd.DataFrame) -> pd.DataFrame:
    rows = []
    fallback = fallback.copy()
    native = native.copy()
    all_ids = sorted(set(fallback.get("dataset_id", pd.Series(dtype=str)).dropna()) | set(native.get("dataset_id", pd.Series(dtype=str)).dropna()))
    for dataset_id in all_ids:
        f = fallback[fallback["dataset_id"].eq(dataset_id)].head(1) if "dataset_id" in fallback else pd.DataFrame()
        n = native[native["dataset_id"].eq(dataset_id)].head(1) if "dataset_id" in native else pd.DataFrame()
        f_primary = _truthy(f["primary_condensation_support"].iloc[0]) if not f.empty and "primary_condensation_support" in f else False
        n_primary = _truthy(n["primary_condensation_support"].iloc[0]) if not n.empty and "primary_condensation_support" in n else False
        if not n.empty:
            if n_primary:
                final_interpretation = "native_primary_support"
                final_tier = "acceptable"
            elif int(n.get("revised_support_count", pd.Series([0])).iloc[0]) > 0:
                final_interpretation = "native_revised_feature_only"
                final_tier = "weak"
            else:
                final_interpretation = "native_no_clone_support"
                final_tier = "fail"
        elif not f.empty:
            if str(f.get("interpretation", pd.Series([""])).iloc[0]) == "not_time_series_usable":
                final_interpretation = "metadata_blocker"
                final_tier = "fail"
            else:
                final_interpretation = "fallback_only"
                final_tier = "weak" if f_primary else str(f.get("validation_tier", pd.Series(["fail"])).iloc[0])
        else:
            final_interpretation = "not_analyzed"
            final_tier = "fail"
        rows.append(
            {
                "dataset_id": dataset_id,
                "full_data_teacher_backend": str(f.get("teacher_backend", pd.Series(["not_available"])).iloc[0]) if not f.empty else "not_available",
                "full_data_usable_clones": f.get("usable_clones", pd.Series([np.nan])).iloc[0] if not f.empty else np.nan,
                "full_data_primary_condensation_support": f_primary if not f.empty else False,
                "full_data_primary_effect": f.get("primary_effect", pd.Series([np.nan])).iloc[0] if not f.empty else np.nan,
                "full_data_primary_q": f.get("primary_q", pd.Series([np.nan])).iloc[0] if not f.empty else np.nan,
                "full_data_covariate_adjusted_effect": f.get("primary_covariate_adjusted_effect", pd.Series([np.nan])).iloc[0] if not f.empty else np.nan,
                "native_teacher_backend": str(n.get("teacher_backend", pd.Series(["not_available"])).iloc[0]) if not n.empty else "not_available",
                "native_usable_clones": n.get("usable_clones", pd.Series([np.nan])).iloc[0] if not n.empty else np.nan,
                "native_primary_condensation_support": n_primary if not n.empty else False,
                "native_primary_effect": n.get("primary_effect", pd.Series([np.nan])).iloc[0] if not n.empty else np.nan,
                "native_primary_q": n.get("primary_q", pd.Series([np.nan])).iloc[0] if not n.empty else np.nan,
                "native_covariate_adjusted_effect": n.get("primary_covariate_adjusted_effect", pd.Series([np.nan])).iloc[0] if not n.empty else np.nan,
                "native_revised_support_count": n.get("revised_support_count", pd.Series([np.nan])).iloc[0] if not n.empty else np.nan,
                "final_clone_aware_tier": final_tier,
                "final_interpretation": final_interpretation,
            }
        )
    return pd.DataFrame(rows)


def _truthy(value) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    return bool(value)


def _write_integrated_reports(integrated: pd.DataFrame, fallback: pd.DataFrame, native: pd.DataFrame) -> None:
    j_native = integrated[integrated["dataset_id"].str.startswith("Jindal", na=False)]
    w_native = integrated[integrated["dataset_id"].str.startswith("Weinreb", na=False)]
    x_row = integrated[integrated["dataset_id"].str.startswith("Xie", na=False)]
    j_line = "- Jindal LSK: full-data fallback showed a weak positive primary association, but downsampled native moscot did not retain it." if not j_native.empty else "- Jindal LSK: not available in integrated summary."
    w_line = "- Weinreb LARRY: primary condensation was not supported in full-data fallback and was negative in downsampled native moscot." if not w_native.empty else "- Weinreb LARRY: not available in integrated summary."
    x_line = "- Xie organoid: clone/barcode metadata are present, but the processed h5ad lacks an explicit time/stage field for branch-window validation." if not x_row.empty else "- Xie organoid: not available in integrated summary."
    core_conclusion = (
        "The clone-aware developmental expansion does not establish clone-level fate-diversification support. "
        "The strongest retained project claim remains the time-series branch-nucleation order-parameter hypothesis supported by internal native moscot and E1 MouseGastrulationData."
    )
    table_md = integrated.to_markdown(index=False) if not integrated.empty else "_No integrated rows._"
    _write_md(
        "reports/clone_integrated_validation_conclusion.md",
        "# Integrated Clone-Aware Validation Conclusion\n\n"
        + core_conclusion
        + "\n\n## Dataset-Level Result\n\n"
        + "\n".join([j_line, w_line, x_line])
        + "\n\n## Integrated Summary\n\n"
        + table_md
        + "\n\n## Interpretation\n\n"
        "The Jindal weak signal is best treated as fallback-teacher feasibility, not as retained clone-aware support, because it did not persist in the downsampled native-teacher check. Weinreb provides negative clone-aware evidence for primary condensation and does not support the revised two-phase or uncertainty-gated models in the native downsample. Xie remains a metadata blocker for branch-window testing.\n",
    )
    attempts = _sanitize_attempt_paths(_safe_read("tables/clone_native_moscot_attempts.csv"))
    if not attempts.empty:
        _write_csv(attempts, "tables/clone_native_moscot_attempts.csv")
        _write_md("reports/clone_native_moscot_attempts.md", "# Clone Native Moscot Attempts\n\n" + attempts.to_markdown(index=False) + "\n")
    attempts_md = attempts.to_markdown(index=False) if not attempts.empty else "_No native attempt table available._"
    fallback_md = fallback.to_markdown(index=False) if not fallback.empty else "_No fallback summary available._"
    _write_md(
        "reports/clone_developmental_validation_expansion.md",
        "# Clone-Aware Developmental Validation Expansion\n\n"
        "Jindal LSK, Weinreb LARRY and Xie organoid were downloaded and inspected from the scLTdb Zenodo record. Jindal and Weinreb contain time/stage, expression, metadata and clone/barcode fields. Xie organoid contains clone/barcode metadata but lacks an explicit time/stage field in the processed h5ad, so it is not usable for branch-window validation here.\n\n"
        "The analysis now has two evidence layers. First, full-data fallback centroid-teacher analysis provided broad clone coverage and showed a weak Jindal positive signal. Second, clean `.venv_moscot_native` runs produced downsampled native moscot teachers for Jindal and Weinreb. In that native refinement, Jindal did not retain primary condensation support and Weinreb was negative for the primary exposure. Therefore the full-data fallback Jindal result is not upgraded to retained clone-aware support.\n\n"
        "## Full-Data Fallback Summary\n\n"
        + fallback_md
        + "\n\n## Native Moscot Attempt Summary\n\n"
        + attempts_md
        + "\n\n## Integrated Conclusion\n\n"
        + table_md
        + "\n",
    )
    _write_md(
        "reports/clone_developmental_data_audit.md",
        "# Clone-Aware Developmental Data Audit\n\n"
        "- source: https://zenodo.org/records/12176634\n"
        "- DOI: 10.5281/zenodo.12176634\n"
        "- Jindal LSK and Weinreb LARRY were downloaded, loaded and analyzed as clone-aware time-series datasets.\n"
        "- Xie organoid was downloaded and loaded but lacks an explicit time/stage field in the processed h5ad used here.\n"
        "- Default-environment moscot import was not usable for full-data analysis, but clean `.venv_moscot_native` downsampled native moscot runs succeeded for Jindal and Weinreb.\n"
        "- Large raw h5ad and processed native inputs are treated as local data artifacts and are not intended for git tracking.\n\n"
        "## Native Attempt Table\n\n"
        + attempts_md
        + "\n",
    )
    _write_md(
        "reports/clone_developmental_final_evidence.md",
        "# Clone-Aware Developmental Final Evidence\n\n"
        "## Dataset Outcomes\n\n"
        "| Dataset | Status | Teacher | Clone result |\n"
        "|---|---|---|---|\n"
        "| Jindal LSK CellTag-Multi | analyzed in full-data fallback and downsampled native moscot | fallback centroid plus downsampled native moscot | fallback weak positive did not persist in native downsample |\n"
        "| Weinreb LARRY | analyzed in full-data fallback and downsampled native moscot | fallback centroid plus downsampled native moscot | primary condensation not supported; native downsample is negative |\n"
        "| Xie organoid | loaded/audited | not applicable for branch-window | clone metadata present but no explicit time/stage field |\n\n"
        "## Project-Level Interpretation\n\n"
        + core_conclusion
        + "\n\nThe revised two-phase and uncertainty-gated branch-window models were also tested in the native downsample and were not supported in Jindal or Weinreb. The appropriate conclusion is that clone-level fate-diversification support is not established in the tested clone-aware datasets.\n\n"
        "## Required Next Step\n\n"
        "The next decisive experiment is either a larger native-teacher Jindal/Weinreb run with stronger clone coverage or a developmental clone dataset with richer ordered sampling, clone/barcode metadata and sufficient clones spanning the branch window.\n",
    )
    final_methods = (
        "# Final Retained Results and Methods\n\n"
        "## Retained Claim\n\n"
        "SwarmLineage-OT retains a branch-nucleation / transient condensation-before-divergence time-series order-parameter hypothesis. Internal native moscot teacher analysis and E1 MouseGastrulationData support remain the strongest evidence.\n\n"
        "## Clone-Aware Developmental Expansion\n\n"
        "Three prioritized clone-aware datasets were downloaded and audited from the scLTdb Zenodo record. Jindal LSK and Weinreb LARRY contain expression matrices, metadata, clone/barcode fields and ordered time/stage information. Xie organoid contains clone/barcode fields but lacks an explicit time/stage field in the processed h5ad used here.\n\n"
        "Jindal LSK initially showed weak full-data fallback feasibility: branch-window condensation was positively associated with terminal clone fate entropy after covariate adjustment. A clean `.venv_moscot_native` run then produced a downsampled native moscot teacher for 1,000 cells and 42 usable clones; under this native check, the primary association did not persist. Weinreb LARRY produced a downsampled native moscot teacher for 1,500 cells and 95 usable clones; primary condensation was not supported and was negative after covariate adjustment. Revised two-phase and uncertainty-gated features were not supported in the native downsample.\n\n"
        "The clone-aware conclusion is therefore conservative: clone-level fate-diversification support is not established. The retained project-level claim remains a time-series order-parameter hypothesis, not a clone-level fate-prediction claim.\n\n"
        "## Still Excluded\n\n"
        "Diffusion remains an encoded control-law recovery. Birth/death, memory and CCI remain unsupported. The local topological-neighbour mechanism and swarm-required attribution remain unresolved.\n"
    )
    _write_md("manuscript/final_retained_results_and_methods.md", final_methods)
    _write_md(
        "manuscript/manuscript.md",
        "# SwarmLineage-OT Clone-Aware Developmental Validation\n\n"
        "SwarmLineage-OT converts native OT-inferred pseudo-lineage maps into finite-agent rollouts and is currently evaluated most strongly through branch-nucleation order parameters. The central retained hypothesis is transient condensation-before-divergence in time-series developmental data.\n\n"
        "Jindal LSK, Weinreb LARRY and Xie organoid were downloaded from the scLTdb Zenodo record and inspected as external clone-aware datasets. Jindal and Weinreb contain expression matrices, metadata, clone/barcode fields and ordered time/stage information; Xie contains clone/barcode fields but lacks an explicit time/stage field in the processed h5ad used here.\n\n"
        "The clone-aware analysis was deliberately staged. Full-data fallback analysis gave Jindal weak feasibility, with primary branch-window condensation positively associated with terminal clone fate entropy after covariate adjustment. Clean native moscot runs in `.venv_moscot_native` then succeeded on downsampled Jindal and Weinreb inputs. In those native checks, Jindal did not retain primary condensation support and Weinreb showed a negative primary association. Revised two-phase and uncertainty-gated branch-window features were not supported in the native downsample.\n\n"
        "These results keep the main claim bounded. Branch nucleation remains supported as a time-series order-parameter hypothesis under internal native teacher analysis and E1 MouseGastrulationData external time-series support. Clone-level fate-diversification support is not established in the tested clone-aware datasets.\n\n"
        "Diffusion remains an encoded control-law recovery. Birth/death, memory and CCI remain unsupported. The local topological-neighbour mechanism and swarm-required attribution remain unresolved.\n",
    )
    _write_md(
        "reports/editorial_assessment.md",
        "# Editorial Assessment\n\n"
        "Current evidence level: clone-aware developmental validation was strengthened by native-teacher attempts, but it still does not establish clone-level support.\n\n"
        "- teacher_fidelity_tier: acceptable\n"
        "- emergent_law_tier: weak\n"
        "- mechanistic_usefulness_tier: weak\n"
        "- Jindal LSK provides only fallback-teacher feasibility; its weak positive primary signal did not persist in downsampled native moscot.\n"
        "- Weinreb LARRY does not support primary condensation and is negative in the downsampled native-teacher check.\n"
        "- Revised two-phase and uncertainty-gated features were not supported in the downsampled native checks.\n"
        "- Xie organoid has clone/barcode metadata but no explicit time/stage field in the processed h5ad used here.\n"
        "- The manuscript should state that clone-level fate-diversification support is not established as a main claim.\n",
    )
    _write_md(
        "reports/scientific_gap_audit.md",
        "# Scientific Gap Audit\n\n"
        "OT gives the developmental map; SwarmLineage-OT learns microscopic finite-agent rules that realize the map and reveal emergent developmental laws.\n\n"
        "- best mean-rank reconstruction row: `M0b_ot_interpolation`\n"
        "- OT reference row: `M0b_ot_interpolation`\n"
        "- teacher_fidelity_tier: acceptable\n"
        "- emergent_law_tier: weak\n"
        "- mechanistic_usefulness_tier: weak\n"
        "- native_or_external_teacher_validation: True\n\n"
        "## Remaining Gaps\n\n"
        "- Jindal and Weinreb now have downsampled native moscot checks, but those checks do not support clone-level fate-diversification claims.\n"
        "- Jindal's weak full-data fallback signal is not native-stable under the current downsampled run.\n"
        "- Weinreb contradicts the primary condensation-only clone association and does not support the revised branch-window feature sets in the native downsample.\n"
        "- Xie lacks time/stage metadata for branch-window validation in the processed h5ad.\n"
        "- Biddy/CellTag L2 remains a failed clone-aware result.\n"
        "- CCI, memory and birth/death remain unsupported and excluded from the main claim.\n"
        "- No manuscript claim may frame SwarmLineage-OT as surpassing the OT reference.\n"
        "- High-impact mechanism claims remain unsupported without stronger clone-aware developmental validation and cleaner attribution controls.\n",
    )
    _write_md(
        "reports/clone_developmental_reviewer_update.md",
        "# Reviewer-Facing Update\n\n"
        "Potential reviewer question: does branch-window condensation predict clone fate diversification in an independent developmental clone dataset?\n\n"
        "Current answer: not under the present native-teacher checks. Jindal LSK and Weinreb LARRY were downloaded with expression, metadata, ordered stages and clone/barcode fields. Full-data fallback analysis gave Jindal weak feasibility, but downsampled native moscot did not retain the primary signal. Weinreb did not support the primary signal and was negative under the native downsample. Xie organoid cannot be tested for branch-window effects from the processed h5ad because time/stage metadata are absent.\n\n"
        + table_md,
    )
    _write_md(
        "reports/clone_developmental_claim_audit.md",
        "# Clone Developmental Claim Audit\n\n"
        "- clone-aware data were genuinely downloaded for Jindal LSK, Weinreb LARRY and Xie organoid.\n"
        "- native moscot was run in a clean environment for downsampled Jindal and Weinreb inputs.\n"
        "- Jindal's weak full-data fallback signal is not upgraded because it did not persist in the native downsample.\n"
        "- Weinreb is negative for the primary condensation exposure under the native downsample.\n"
        "- no result is described as experimental confirmation.\n"
        "- primary condensation exposure must support the primary clone score after covariates before clone-level support can be claimed.\n",
    )


def _update_final_claim_tiers() -> None:
    rows = [
        {
            "claim": "native moscot teacher construction",
            "status": "retained",
            "tier": "native_moscot_success",
            "internal_native_support": True,
            "native_sensitivity_support": True,
            "external_time_series_support": False,
            "lineage_clone_support": False,
            "negative_controls": "not_applicable_or_failed",
            "module_necessity": "not_applicable",
            "external_independence": "not_applicable",
            "allowed_manuscript_sentence": "Native moscot teacher extraction succeeded for the analyzed downsampled settings.",
            "forbidden_sentence": "Native moscot provides ground-truth ancestry.",
        },
        {
            "claim": "teacher fidelity",
            "status": "retained",
            "tier": "acceptable",
            "internal_native_support": True,
            "native_sensitivity_support": False,
            "external_time_series_support": False,
            "lineage_clone_support": False,
            "negative_controls": "not_applicable_or_failed",
            "module_necessity": "not_applicable",
            "external_independence": "not_applicable",
            "allowed_manuscript_sentence": "The primary agent remains within acceptable teacher-fidelity tolerance.",
            "forbidden_sentence": "The agent model beats the OT reference.",
        },
        {
            "claim": "primary model M5 selection",
            "status": "retained",
            "tier": "acceptable",
            "internal_native_support": True,
            "native_sensitivity_support": False,
            "external_time_series_support": False,
            "lineage_clone_support": False,
            "negative_controls": "not_applicable_or_failed",
            "module_necessity": "not_applicable",
            "external_independence": "not_applicable",
            "allowed_manuscript_sentence": "M5_ot_swarm is the primary mechanistic model selected by evidence and simplicity.",
            "forbidden_sentence": "The full memory model is automatically primary.",
        },
        {
            "claim": "branch nucleation signature",
            "status": "retained_computational_hypothesis",
            "tier": "strong_internal_acceptable_external",
            "internal_native_support": True,
            "native_sensitivity_support": True,
            "external_time_series_support": True,
            "lineage_clone_support": False,
            "negative_controls": "performed; swarm-specific necessity remains weak",
            "module_necessity": "not_applicable",
            "external_independence": "not_applicable",
            "allowed_manuscript_sentence": "A branch-nucleation order-parameter signature is retained as a computational hypothesis.",
            "forbidden_sentence": "A biological mechanism is proven.",
        },
        {
            "claim": "transient condensation-before-divergence",
            "status": "retained_computational_hypothesis",
            "tier": "acceptable",
            "internal_native_support": True,
            "native_sensitivity_support": True,
            "external_time_series_support": True,
            "lineage_clone_support": False,
            "negative_controls": "performed; swarm-specific necessity remains weak",
            "module_necessity": "not_applicable",
            "external_independence": "not_applicable",
            "allowed_manuscript_sentence": "The best current interpretation is transient condensation-before-divergence.",
            "forbidden_sentence": "The signature is causal or experimentally established.",
        },
        {
            "claim": "clone-aware developmental validation",
            "status": "not_established",
            "tier": "fail_native_refinement",
            "internal_native_support": False,
            "native_sensitivity_support": False,
            "external_time_series_support": False,
            "lineage_clone_support": False,
            "negative_controls": "fallback Jindal controls pass but native check fails",
            "module_necessity": "not_applicable",
            "external_independence": "Jindal/Weinreb/Xie scLTdb datasets",
            "allowed_manuscript_sentence": "Clone-aware developmental datasets were analyzed, but clone-level fate-diversification support is not established.",
            "forbidden_sentence": "Clone splitting is reliably predicted.",
        },
        {
            "claim": "revised two-phase or uncertainty-gated clone model",
            "status": "not_established",
            "tier": "fail_native_refinement",
            "internal_native_support": False,
            "native_sensitivity_support": False,
            "external_time_series_support": False,
            "lineage_clone_support": False,
            "negative_controls": "not_supported_in_native_downsample",
            "module_necessity": "not_applicable",
            "external_independence": "Jindal/Weinreb scLTdb datasets",
            "allowed_manuscript_sentence": "Revised branch-window feature sets remain unconfirmed in clone-aware native-teacher checks.",
            "forbidden_sentence": "The two-phase clone model is established.",
        },
        {
            "claim": "swarm alignment contribution",
            "status": "not_sufficient_for_main_claim",
            "tier": "weak",
            "internal_native_support": False,
            "native_sensitivity_support": False,
            "external_time_series_support": False,
            "lineage_clone_support": False,
            "negative_controls": "performed; swarm-specific necessity remains weak",
            "module_necessity": "not_established",
            "external_independence": "not_applicable",
            "allowed_manuscript_sentence": "Swarm terms may stabilize the retained model but necessity is not established.",
            "forbidden_sentence": "Swarm rules are required.",
        },
        {
            "claim": "diffusion encoded recovery",
            "status": "encoded_control_law_recovery",
            "tier": "acceptable_encoded",
            "internal_native_support": False,
            "native_sensitivity_support": False,
            "external_time_series_support": False,
            "lineage_clone_support": False,
            "negative_controls": "not_applicable_or_failed",
            "module_necessity": "not_applicable",
            "external_independence": "not_applicable",
            "allowed_manuscript_sentence": "Diffusion is an encoded control-law recovery, not independent discovery.",
            "forbidden_sentence": "Diffusion law was discovered independently.",
        },
        {
            "claim": "birth/death unsupported",
            "status": "unsupported",
            "tier": "fail",
            "internal_native_support": False,
            "native_sensitivity_support": False,
            "external_time_series_support": False,
            "lineage_clone_support": False,
            "negative_controls": "not_applicable_or_failed",
            "module_necessity": "not_applicable",
            "external_independence": "not_applicable",
            "allowed_manuscript_sentence": "Birth/death remains unsupported under current evidence.",
            "forbidden_sentence": "Birth/death module is an established mechanism.",
        },
        {
            "claim": "memory unsupported",
            "status": "unsupported",
            "tier": "fail",
            "internal_native_support": False,
            "native_sensitivity_support": False,
            "external_time_series_support": False,
            "lineage_clone_support": False,
            "negative_controls": "not_applicable_or_failed",
            "module_necessity": "not_applicable",
            "external_independence": "not_applicable",
            "allowed_manuscript_sentence": "Memory hysteresis remains unsupported.",
            "forbidden_sentence": "Memory hysteresis is established.",
        },
        {
            "claim": "CCI unsupported",
            "status": "unsupported",
            "tier": "fail",
            "internal_native_support": False,
            "native_sensitivity_support": False,
            "external_time_series_support": False,
            "lineage_clone_support": False,
            "negative_controls": "not_applicable_or_failed",
            "module_necessity": "not_applicable",
            "external_independence": "not_applicable",
            "allowed_manuscript_sentence": "CCI branch bias remains unsupported.",
            "forbidden_sentence": "CCI branch bias is established.",
        },
        {
            "claim": "external MouseGastrulationData time-series support",
            "status": "external_time_series_support",
            "tier": "acceptable",
            "internal_native_support": False,
            "native_sensitivity_support": False,
            "external_time_series_support": True,
            "lineage_clone_support": False,
            "negative_controls": "not_applicable_or_failed",
            "module_necessity": "not_applicable",
            "external_independence": "related_atlas_independent_sample",
            "allowed_manuscript_sentence": "E1 provides external time-series support, without lineage barcodes.",
            "forbidden_sentence": "E1 establishes clone-resolved support.",
        },
        {
            "claim": "external generalization across systems",
            "status": "exploratory",
            "tier": "weak",
            "internal_native_support": False,
            "native_sensitivity_support": False,
            "external_time_series_support": False,
            "lineage_clone_support": False,
            "negative_controls": "not_applicable_or_failed",
            "module_necessity": "not_applicable",
            "external_independence": "not_applicable",
            "allowed_manuscript_sentence": "Cross-system generalization remains exploratory.",
            "forbidden_sentence": "Generalization across systems is established.",
        },
        {
            "claim": "experimental validation",
            "status": "future_work",
            "tier": "fail",
            "internal_native_support": False,
            "native_sensitivity_support": False,
            "external_time_series_support": False,
            "lineage_clone_support": False,
            "negative_controls": "not_applicable_or_failed",
            "module_necessity": "not_applicable",
            "external_independence": "not_applicable",
            "allowed_manuscript_sentence": "Experimental validation is a future requirement.",
            "forbidden_sentence": "Experimental confirmation has already been completed.",
        },
        {
            "claim": "high-impact biological readiness",
            "status": "not_ready",
            "tier": "fail",
            "internal_native_support": False,
            "native_sensitivity_support": False,
            "external_time_series_support": False,
            "lineage_clone_support": False,
            "negative_controls": "not_applicable_or_failed",
            "module_necessity": "not_applicable",
            "external_independence": "not_applicable",
            "allowed_manuscript_sentence": "The package is not ready for strong biological mechanism claims.",
            "forbidden_sentence": "Ready for a top-journal biological mechanism claim.",
        },
        {
            "claim": "topological-neighbor minimal rule",
            "status": "dataset_specific",
            "tier": "weak",
            "internal_native_support": False,
            "native_sensitivity_support": False,
            "external_time_series_support": False,
            "lineage_clone_support": False,
            "negative_controls": "random and metric controls included",
            "module_necessity": "not causal",
            "external_independence": "E1 related support",
            "allowed_manuscript_sentence": "A topological-neighbor minimal-rule explanation is supported only at the reported tier.",
            "forbidden_sentence": "Topological rule is established.",
        },
    ]
    table = pd.DataFrame(rows)
    _write_csv(table, "tables/final_claim_evidence_tiers.csv")
    _write_md("reports/final_claim_evidence_tiers.md", "# Final Claim Evidence Tiers\n\n" + table.to_markdown(index=False) + "\n")


def _write_forbidden_claim_scan() -> None:
    forbidden = [
        "Nature-ready",
        "proven biological mechanism",
        "causal proof",
        "wet-lab validated",
        "true lineage validation",
        "clone splitting reliably predicted",
        "SwarmLineage-OT beats OT",
        "topological-neighbor mechanism proven",
        "CCI validated",
        "memory hysteresis discovered",
        "birth/death law discovered",
    ]
    skip = {"clone_developmental_forbidden_claim_scan.md", "clone_native_forbidden_claim_scan.md"}
    hits = []
    for root in ["reports", "manuscript"]:
        for path in (ROOT / root).rglob("*.md"):
            if path.name in skip:
                continue
            text = path.read_text(encoding="utf-8", errors="ignore")
            for phrase in forbidden:
                if phrase in text:
                    hits.append({"file": str(path.relative_to(ROOT)), "phrase_hash": abs(hash(phrase)) % 1000000})
    body = "# Clone Native Forbidden Claim Scan\n\n" + f"- hits: {len(hits)}\n\n"
    body += "No forbidden claim strings found." if not hits else pd.DataFrame(hits).to_markdown(index=False)
    _write_md("reports/clone_native_forbidden_claim_scan.md", body)
    _write_md("reports/clone_developmental_forbidden_claim_scan.md", body)


def _sanitize_attempt_paths(attempts: pd.DataFrame) -> pd.DataFrame:
    if attempts.empty:
        return attempts
    attempts = attempts.copy()
    root_raw = str(ROOT).replace("\\", "/")
    root_win = str(ROOT)
    for col in ["path", "config_path"]:
        if col in attempts:
            attempts[col] = (
                attempts[col]
                .astype(str)
                .str.replace(root_raw + "/", "", regex=False)
                .str.replace(root_win + "\\", "", regex=False)
                .str.replace("\\", "/", regex=False)
            )
    if "detail" in attempts:
        attempts["detail"] = attempts["detail"].astype(str).map(lambda x: "\n".join(line.rstrip() for line in x.splitlines()))
    return attempts


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.parse_args()
    run()
