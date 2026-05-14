from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from src.ot_teacher.build_teacher import build_teacher
from src.ot_teacher.diagnostics import js_divergence
from src.utils.config import ensure_dir, load_config, write_text


def _tag(cells: int, epsilon: float, iterations: int) -> str:
    return f"cells{cells}_eps{epsilon:.2f}_iter{iterations}".replace(".", "p")


def _safe_load_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _write_yaml(path: Path, cfg: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")


def _run_native(base_cfg: dict, run_cfg: dict, venv_python: str, timeout: int) -> dict:
    cfg_path = Path(run_cfg["run_dir"]) / "ot_teacher_native_config.yaml"
    _write_yaml(cfg_path, run_cfg)
    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path.cwd())
    env["JAX_PLATFORMS"] = "cpu"
    env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    start = time.perf_counter()
    proc = subprocess.run(
        [venv_python, "-m", "src.ot_teacher.run_moscot", "--config", str(cfg_path)],
        cwd=Path.cwd(),
        env=env,
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    elapsed = time.perf_counter() - start
    summary_path = Path(run_cfg["couplings_dir"]) / "moscot_run_summary.json"
    summary = _safe_load_json(summary_path)
    status = "native_moscot_success" if proc.returncode == 0 and summary.get("teacher_backend") == "native_moscot" else "native_moscot_failed"
    return {
        "status": status,
        "returncode": int(proc.returncode),
        "runtime_seconds": float(elapsed),
        "stdout_tail": proc.stdout[-2000:],
        "stderr_tail": proc.stderr[-2000:],
        "summary": summary,
        "summary_path": str(summary_path),
        "config_path": str(cfg_path),
    }


def _build_variant_teacher(run_cfg: dict) -> dict:
    return build_teacher(run_cfg)


def _lineage_edge_matrix(npz_path: str | Path) -> pd.DataFrame:
    raw = np.load(npz_path, allow_pickle=True)
    plan = raw["plan"].astype(float)
    src = raw["source_types"].astype(str)
    tgt = raw["target_types"].astype(str)
    rows = []
    for s in sorted(set(src)):
        for t in sorted(set(tgt)):
            mass = float(plan[np.ix_(src == s, tgt == t)].sum())
            rows.append({"source_lineage": s, "target_lineage": t, "mass": mass})
    return pd.DataFrame(rows)


def _pair_metrics(index: pd.DataFrame, ref_index: pd.DataFrame | None = None) -> tuple[pd.DataFrame, dict]:
    rows = []
    edge_stabilities = []
    for row in index.itertuples(index=False):
        raw = np.load(row.file, allow_pickle=True)
        plan = raw["plan"].astype(float)
        source_marg = plan.sum(axis=1)
        target_marg = plan.sum(axis=0)
        rec = {
            "source_time": float(row.source_time),
            "target_time": float(row.target_time),
            "plan_shape": f"{plan.shape[0]}x{plan.shape[1]}",
            "n_source": int(plan.shape[0]),
            "n_target": int(plan.shape[1]),
            "converged": bool(getattr(row, "converged", True)),
            "transport_mass": float(plan.sum()),
            "source_mass_cv": float(source_marg.std() / max(source_marg.mean(), 1e-12)),
            "target_mass_cv": float(target_marg.std() / max(target_marg.mean(), 1e-12)),
            "mean_entropy": float(getattr(row, "mean_entropy", np.nan)),
            "transport_cost": float(getattr(row, "transport_cost", np.nan)),
        }
        if ref_index is not None and not ref_index.empty:
            hit = ref_index[(ref_index["source_time"].astype(float) == rec["source_time"]) & (ref_index["target_time"].astype(float) == rec["target_time"])]
            if not hit.empty:
                ref_edges = _lineage_edge_matrix(hit.iloc[0]["file"])
                cur_edges = _lineage_edge_matrix(row.file)
                key = ["source_lineage", "target_lineage"]
                merged = ref_edges.merge(cur_edges, on=key, how="outer", suffixes=("_ref", "_cur")).fillna(0.0)
                edge_js = js_divergence(merged["mass_ref"].to_numpy(float), merged["mass_cur"].to_numpy(float))
                rec["lineage_edge_mass_js_vs_reference"] = float(edge_js)
                edge_stabilities.append(1.0 - min(float(edge_js), 1.0))
        rows.append(rec)
    pair_frame = pd.DataFrame(rows)
    summary = {
        "mean_pair_entropy": float(pair_frame["mean_entropy"].mean()) if "mean_entropy" in pair_frame else np.nan,
        "mean_lineage_edge_stability": float(np.nanmean(edge_stabilities)) if edge_stabilities else np.nan,
    }
    return pair_frame, summary


def _teacher_compare(teacher_path: Path, reference_path: Path) -> dict:
    cur = ad.read_h5ad(teacher_path)
    ref = ad.read_h5ad(reference_path)
    zc = np.asarray(cur.obsm["X_ot_velocity"], dtype=float)
    zr = np.asarray(ref.obsm["X_ot_velocity"], dtype=float)
    valid = np.isfinite(zc).all(axis=1) & np.isfinite(zr).all(axis=1)
    dot = np.sum(zc[valid] * zr[valid], axis=1)
    denom = np.linalg.norm(zc[valid], axis=1) * np.linalg.norm(zr[valid], axis=1)
    cos = dot / np.maximum(denom, 1e-12)
    fate_cols = [c for c in cur.obs.columns if c.startswith("fate_prob_") and c in ref.obs.columns]
    return {
        "barycentric_velocity_cosine_mean": float(np.nanmean(cos)) if cos.size else np.nan,
        "barycentric_velocity_rmse": float(np.sqrt(np.nanmean((zc[valid] - zr[valid]) ** 2))) if valid.any() else np.nan,
        "entropy_rmse": float(
            np.sqrt(np.nanmean((cur.obs["ot_transition_entropy"].to_numpy(float) - ref.obs["ot_transition_entropy"].to_numpy(float)) ** 2))
        ),
        "entropy_correlation": float(
            np.corrcoef(cur.obs["ot_transition_entropy"].to_numpy(float), ref.obs["ot_transition_entropy"].to_numpy(float))[0, 1]
        ),
        "fate_probability_mae": float(np.mean(np.abs(cur.obs[fate_cols].to_numpy(float) - ref.obs[fate_cols].to_numpy(float)))) if fate_cols else np.nan,
    }


def _static_branch_proxy(teacher_path: Path, time_key: str = "time_numeric", cell_type_key: str = "lineage") -> dict:
    adata = ad.read_h5ad(teacher_path)
    z = np.asarray(adata.obsm["X_ot_velocity"], dtype=float)
    obs = adata.obs.copy()
    rows = []
    for time_value, group in obs.groupby(time_key, observed=False):
        labels = group[cell_type_key].astype(str)
        idx = obs.index.get_indexer(group.index)
        centers = []
        for lab in sorted(labels.unique()):
            mask = labels.to_numpy() == lab
            if mask.any():
                centers.append(z[idx[mask]].mean(axis=0))
        separation = 0.0
        if len(centers) > 1:
            arr = np.asarray(centers)
            d = arr[:, None, :] - arr[None, :, :]
            separation = float(np.sqrt((d**2).sum(axis=2)).mean())
        counts = labels.value_counts(normalize=True)
        rows.append(
            {
                "time": float(time_value),
                "lineage_velocity_separation": separation,
                "fate_entropy": float(group["ot_transition_entropy"].mean()),
                "branch_imbalance": float(counts.max() - counts.min()) if counts.size > 1 else 1.0,
            }
        )
    frame = pd.DataFrame(rows).sort_values("time")
    if frame.shape[0] < 3:
        return {"branch_proxy_effect": np.nan, "branch_proxy_interpretation": "insufficient_timepoints"}
    diff = frame["lineage_velocity_separation"].diff().dropna()
    effect = float(diff.min())
    return {
        "branch_proxy_effect": effect,
        "branch_proxy_interpretation": "condensation_proxy" if effect < 0 else "divergence_proxy",
    }


def run(config: str, quick: bool = False) -> dict:
    base_cfg = load_config(config)
    table_dir = ensure_dir("tables")
    report_dir = ensure_dir("reports")
    fig_dir = ensure_dir("figures/discovery")
    root = ensure_dir("processed/native_sensitivity")
    venv_python = str(Path(".venv_moscot_native/Scripts/python.exe"))
    if not Path(venv_python).exists():
        venv_python = sys.executable
    cells_grid = [120, 250] if quick else [120, 250, 500, 650]
    eps_grid = [0.08] if quick else [0.04, 0.08, 0.12]
    iterations = int(base_cfg.get("native_max_iterations", 350))
    timeout = 420 if quick else 900

    sensitivity_rows = []
    pair_frames = []
    reference_teacher = Path("processed/ot_teacher.h5ad")
    reference_index: pd.DataFrame | None = None
    if Path("processed/ot_couplings/teacher_coupling_index.csv").exists():
        reference_index = pd.read_csv("processed/ot_couplings/teacher_coupling_index.csv")
    for cells in cells_grid:
        for eps in eps_grid:
            tag = _tag(cells, eps, iterations)
            run_dir = ensure_dir(root / tag)
            run_cfg = dict(base_cfg)
            run_cfg.update(
                {
                    "run_dir": str(run_dir),
                    "use_native_moscot": True,
                    "native_max_cells_per_time": int(cells),
                    "epsilon": float(eps),
                    "native_max_iterations": iterations,
                    "native_jit": False,
                    "native_device": "cpu",
                    "couplings_dir": str(run_dir / "ot_couplings"),
                    "teacher_path": str(run_dir / "ot_teacher.h5ad"),
                    "fate_probabilities_path": str(run_dir / "ot_fate_probabilities.parquet"),
                    "teacher_index_path": str(run_dir / "ot_couplings" / "teacher_coupling_index.csv"),
                    "summary_path": str(run_dir / "ot_teacher_summary.json"),
                    "table_dir": str(run_dir / "tables"),
                    "report_dir": str(run_dir / "reports"),
                    "figure_dir": str(run_dir / "figures"),
                }
            )
            started = time.perf_counter()
            try:
                run_result = _run_native(base_cfg, run_cfg, venv_python, timeout)
            except subprocess.TimeoutExpired as exc:
                run_result = {
                    "status": "native_moscot_failed",
                    "returncode": -1,
                    "runtime_seconds": float(timeout),
                    "stdout_tail": (exc.stdout or "")[-2000:] if isinstance(exc.stdout, str) else "",
                    "stderr_tail": (exc.stderr or "")[-2000:] if isinstance(exc.stderr, str) else "",
                    "summary": {},
                    "summary_path": str(run_dir / "ot_couplings" / "moscot_run_summary.json"),
                    "config_path": str(run_dir / "ot_teacher_native_config.yaml"),
                    "failure_reason": f"timeout>{timeout}s",
                }
            row = {
                "native_max_cells_per_time": int(cells),
                "epsilon": float(eps),
                "native_max_iterations": iterations,
                "native_jit": False,
                "status": run_result["status"],
                "fallback_used": False,
                "fallback_not_used": True,
                "runtime_seconds": float(run_result["runtime_seconds"]),
                "config_path": run_result["config_path"],
                "summary_path": run_result["summary_path"],
                "failure_reason": run_result.get("failure_reason", ""),
            }
            summary = run_result.get("summary", {})
            if run_result["status"] == "native_moscot_success":
                try:
                    build_summary = _build_variant_teacher(run_cfg)
                    index = pd.read_csv(run_cfg["teacher_index_path"])
                    pair_metrics, pair_summary = _pair_metrics(index, reference_index)
                    pair_metrics.insert(0, "epsilon", float(eps))
                    pair_metrics.insert(0, "native_max_cells_per_time", int(cells))
                    pair_frames.append(pair_metrics)
                    row.update(
                        {
                            "teacher_backend": build_summary.get("teacher_backend", summary.get("teacher_backend")),
                            "pair_count": int(index.shape[0]),
                            "plan_shapes": ";".join(sorted(pair_metrics["plan_shape"].unique())),
                            "all_pairs_converged": bool(pair_metrics["converged"].all()),
                            **pair_summary,
                            **_teacher_compare(Path(run_cfg["teacher_path"]), reference_teacher),
                            **_static_branch_proxy(Path(run_cfg["teacher_path"]), base_cfg.get("time_key", "time_numeric"), base_cfg.get("cell_type_key", "lineage")),
                        }
                    )
                except Exception as exc:
                    row.update({"status": "native_moscot_failed", "failure_reason": f"postprocess:{type(exc).__name__}:{exc}"})
            else:
                row.update(
                    {
                        "teacher_backend": summary.get("teacher_backend", "none"),
                        "pair_count": 0,
                        "plan_shapes": "",
                        "all_pairs_converged": False,
                    }
                )
            sensitivity_rows.append(row)
            # Stop resource-heavy grid if 500+ fails at all eps; lower settings remain retained.
            if time.perf_counter() - started > timeout:
                break

    sens = pd.DataFrame(sensitivity_rows)
    sens["sensitivity_incomplete_due_to_resource"] = sens["status"].ne("native_moscot_success")
    sens.to_csv(table_dir / "native_teacher_sensitivity.csv", index=False)
    pair_all = pd.concat(pair_frames, ignore_index=True) if pair_frames else pd.DataFrame()
    pair_all.to_csv(table_dir / "native_teacher_pair_metrics.csv", index=False)

    law = sens[sens["status"] == "native_moscot_success"][
        ["native_max_cells_per_time", "epsilon", "branch_proxy_effect", "branch_proxy_interpretation", "barycentric_velocity_cosine_mean", "entropy_correlation"]
    ].copy()
    if not law.empty:
        law["branch_nucleation_proxy_stable"] = law["branch_proxy_effect"].lt(0)
        stable = bool(law["branch_nucleation_proxy_stable"].all()) and law["native_max_cells_per_time"].nunique() >= 2
        law["branch_nucleation_teacher_sensitivity_tier"] = "acceptable" if stable else "weak"
    law.to_csv(table_dir / "native_teacher_law_stability.csv", index=False)
    law.to_csv(table_dir / "branch_nucleation_teacher_sensitivity.csv", index=False)

    fig, ax = plt.subplots(figsize=(7, 4), dpi=160)
    ok = sens[sens["status"] == "native_moscot_success"]
    if not ok.empty:
        for eps, group in ok.groupby("epsilon"):
            ax.plot(group["native_max_cells_per_time"], group["barycentric_velocity_cosine_mean"], marker="o", label=f"epsilon={eps:g}")
    ax.set_xlabel("native max cells per time")
    ax.set_ylabel("velocity cosine vs reference")
    ax.set_ylim(-0.05, 1.05)
    ax.set_title("Native teacher barycentric stability")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(fig_dir / "native_teacher_barycentric_stability.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4), dpi=160)
    if not ok.empty:
        scatter = ax.scatter(ok["runtime_seconds"], ok["entropy_rmse"], c=ok["native_max_cells_per_time"], cmap="viridis", s=50)
        fig.colorbar(scatter, ax=ax, label="max cells/time")
    fail = sens[sens["status"] != "native_moscot_success"]
    if not fail.empty:
        ax.scatter(fail["runtime_seconds"], np.zeros(fail.shape[0]), marker="x", color="#b14f4f", label="failed")
    ax.set_xlabel("runtime seconds")
    ax.set_ylabel("entropy RMSE vs reference")
    ax.set_title("Native teacher sensitivity summary")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(fig_dir / "native_teacher_sensitivity.png")
    plt.close(fig)

    success_count = int(sens["status"].eq("native_moscot_success").sum())
    largest = int(ok["native_max_cells_per_time"].max()) if not ok.empty else 0
    write_text(
        report_dir / "native_teacher_sensitivity.md",
        "\n".join(
            [
                "# Native Teacher Sensitivity",
                "",
                "Native moscot TemporalProblem sensitivity was evaluated without falling back to toy Sinkhorn. Failed settings are retained as failures.",
                "",
                f"- successful settings: {success_count}/{sens.shape[0]}",
                f"- largest successful max cells per time: {largest}",
                f"- branch proxy stable across successful settings: {bool(law['branch_nucleation_proxy_stable'].all()) if not law.empty else False}",
                "",
                "## Settings",
                "",
                sens.to_markdown(index=False),
                "",
                "## Branch Nucleation Teacher Sensitivity",
                "",
                law.to_markdown(index=False) if not law.empty else "No successful native sensitivity run was available.",
                "",
            ]
        ),
    )
    return {"success_settings": success_count, "largest_successful_cells": largest, "rows": int(sens.shape[0])}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/ot_teacher_native.yaml")
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.quick), indent=2))


if __name__ == "__main__":
    main()
