"""Run the DevVCell reproducibility pipeline.

The script provides one entry point for regenerating the project outputs used
by the Chinese manuscript. It intentionally exposes a quick mode for validating
an existing run without repeating heavy model training.
"""

from __future__ import annotations

import argparse
import json
import platform
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


CRITICAL_OUTPUTS = [
    "results/tables/stage_vulnerability.csv",
    "results/tables/perturbation_priority.csv",
    "results/tables/ot_grn_tf_system_developmental_impact.csv",
    "data/processed/cell_level_subset_v1.h5ad",
    "results/cell_level_v1/tables/cell_level_transition_metrics.csv",
    "results/cell_level_v1/tables/transition_statistical_summary.csv",
    "results/cell_level_v1/tables/transition_paired_bootstrap_differences.csv",
    "results/cell_level_v1/figures/transition_bootstrap_ci.png",
    "results/cell_level_v1/models/transition_context_residual_mlp.pt",
    "results/cell_level_v1/tables/cell_level_tf_grn_stimulus_summary.csv",
    "data/external/scperturb/scperturb_benchmark_manifest.json",
    "results/external_perturbation_v1/tables/external_perturbation_metrics.csv",
    "results/external_perturbation_v1/tables/external_perturbation_condition_metrics.csv",
    "results/external_perturbation_v1/figures/external_perturbation_benchmark_mse.png",
    "results/external_perturbation_v1/external_perturbation_summary.json",
    "config/external_perturbation_benchmark.json",
    "results/ablation_v1/tables/transition_ablation_summary.csv",
    "results/ablation_v1/tables/stimulus_ablation_summary.csv",
    "results/nature_virtual_cell_evidence/evidence_manifest.json",
    "results/nature_virtual_cell_evidence/paper_methods_evidence_pipeline.md",
    "results/nature_virtual_cell_evidence/tables/input_schema_validation.csv",
    "results/nature_virtual_cell_evidence/tables/claim_gate_matrix.csv",
    "results/nature_virtual_cell_evidence/tables/fate_recovery_virtual_screen.csv",
    "results/nature_virtual_cell_evidence/tables/external_perturbation_bias_control_summary.csv",
    "results/nature_virtual_cell_evidence/figures/figure_competence_window_response.png",
    "manuscript/figures/devvcell_framework_cn_ai.png",
    "manuscript/main_cn.tex",
    "manuscript/main_cn.pdf",
    "docs/COMPLETION_AUDIT_CN.md",
    "docs/EXTERNAL_PERTURBATION_BENCHMARK_CN.md",
    "docs/PEER_METHOD_COMPARISON_CN.md",
    "docs/REPRODUCIBILITY_CN.md",
    "docs/NATURE_DELIVERY_ROADMAP_CN.md",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=["quick", "main", "full"],
        default="quick",
        help=(
            "quick validates existing outputs and recompiles the paper; main reruns "
            "the primary analyses; full also reruns all ablations."
        ),
    )
    parser.add_argument("--cell-config", default="config/cell_level_baseline.json")
    parser.add_argument("--ablation-config", default="config/ablation_suite.json")
    parser.add_argument("--external-config", default="config/external_perturbation_benchmark.json")
    parser.add_argument("--force-subset", action="store_true", help="Re-export the cell-level subset.")
    parser.add_argument("--skip-ablation", action="store_true", help="Skip ablation aggregation/runs.")
    parser.add_argument("--skip-external-benchmark", action="store_true", help="Skip external scPerturb ingest and benchmark.")
    parser.add_argument("--no-compile-paper", action="store_true", help="Do not run xelatex.")
    return parser.parse_args()


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT_ROOT)).replace("\\", "/")
    except ValueError:
        return str(path)


def run_step(name: str, command: list[str], cwd: Path = PROJECT_ROOT) -> dict[str, object]:
    print(f"[DevVCell] {name}: {' '.join(command)}")
    started = time.perf_counter()
    completed = subprocess.run(command, cwd=cwd, check=True)
    elapsed = time.perf_counter() - started
    return {
        "name": name,
        "command": command,
        "cwd": rel(cwd),
        "returncode": int(completed.returncode),
        "elapsed_seconds": round(elapsed, 3),
    }


def output_inventory(paths: list[str]) -> list[dict[str, object]]:
    rows = []
    for item in paths:
        path = PROJECT_ROOT / item
        rows.append(
            {
                "path": item,
                "exists": path.exists(),
                "size_bytes": int(path.stat().st_size) if path.exists() and path.is_file() else None,
                "modified_utc": (
                    datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat()
                    if path.exists()
                    else None
                ),
            }
        )
    return rows


def compile_paper(steps: list[dict[str, object]]) -> None:
    xelatex = shutil.which("xelatex")
    if xelatex is None:
        raise RuntimeError("xelatex is not available; install TeX Live or run with --no-compile-paper.")
    manuscript_dir = PROJECT_ROOT / "manuscript"
    command = [xelatex, "-interaction=nonstopmode", "-halt-on-error", "main_cn.tex"]
    steps.append(run_step("compile paper pass 1", command, cwd=manuscript_dir))
    steps.append(run_step("compile paper pass 2", command, cwd=manuscript_dir))


def main() -> None:
    args = parse_args()
    steps: list[dict[str, object]] = []

    if args.mode in {"main", "full"}:
        steps.append(run_step("prototype metrics", [sys.executable, "scripts/devvcell_lite.py"]))
        steps.append(run_step("OT+GRN evidence", [sys.executable, "scripts/ot_grn_developmental_impact.py"]))

        subset_path = PROJECT_ROOT / "data" / "processed" / "cell_level_subset_v1.h5ad"
        if args.force_subset or not subset_path.exists():
            command = [sys.executable, "scripts/export_cell_level_subset.py", "--config", args.cell_config]
            if args.force_subset:
                command.append("--force")
            steps.append(run_step("export cell-level subset", command))

        steps.append(
            run_step(
                "train cell-level transition baseline",
                [sys.executable, "scripts/train_cell_transition_baseline.py", "--config", args.cell_config],
            )
        )
        steps.append(
            run_step(
                "run TF/GRN stimulus head",
                [sys.executable, "scripts/run_stimulus_response_head.py", "--config", args.cell_config],
            )
        )

    if not args.skip_external_benchmark:
        steps.append(
            run_step(
                "external scPerturb benchmark ingest",
                [
                    sys.executable,
                    "scripts/prepare_external_perturbation_benchmark.py",
                    "--record-id",
                    "10044268",
                    "--file-name",
                    "DatlingerBock2021.h5ad",
                    "--output-dir",
                    "data/external/scperturb",
                ],
            )
        )
        steps.append(
            run_step(
                "external scPerturb guide-transfer benchmark",
                [
                    sys.executable,
                    "scripts/run_external_perturbation_benchmark.py",
                    "--config",
                    args.external_config,
                ],
            )
        )

    steps.append(
        run_step(
            "transition statistics",
            [
                sys.executable,
                "scripts/summarize_transition_statistics.py",
                "--metrics",
                "results/cell_level_v1/tables/cell_level_transition_metrics.csv",
                "--output-dir",
                "results/cell_level_v1",
                "--reference-model",
                "context_residual_mlp",
            ],
        )
    )

    if not args.skip_ablation:
        ablation_command = [sys.executable, "scripts/run_ablation_suite.py", "--config", args.ablation_config]
        if args.mode != "full":
            ablation_command.extend(["--skip-transition", "--skip-stimulus"])
        steps.append(run_step("ablation suite", ablation_command))

    steps.append(
        run_step(
            "Formal Nature virtual-cell evidence pipeline",
            [sys.executable, "scripts/run_nature_virtual_cell_evidence_pipeline.py", "--config", "config/nature_virtual_cell_evidence.json"],
        )
    )

    if not args.no_compile_paper:
        compile_paper(steps)

    inventory = output_inventory(CRITICAL_OUTPUTS)
    missing = [row["path"] for row in inventory if not row["exists"]]
    manifest = {
        "analysis": "devvcell_reproducibility_pipeline",
        "mode": args.mode,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "steps": steps,
        "critical_outputs": inventory,
        "missing_critical_outputs": missing,
    }
    manifest_path = PROJECT_ROOT / "results" / "reproducibility_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)

    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    if missing:
        raise RuntimeError(f"Missing critical outputs: {missing}")


if __name__ == "__main__":
    main()
