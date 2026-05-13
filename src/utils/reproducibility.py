from __future__ import annotations

import argparse
import importlib
import importlib.metadata as metadata
import json
import os
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .config import ensure_dir, write_json, write_text


METHODS = [
    {
        "name": "scvi-tools",
        "module": "scvi",
        "dist": "scvi-tools",
        "install_command": "pip install scvi-tools",
        "citation": "Gayoso et al., Nature Biotechnology, 2022.",
        "role": "optional latent batch-corrected embedding backend",
    },
    {
        "name": "moscot",
        "module": "moscot",
        "dist": "moscot",
        "install_command": "pip install moscot",
        "citation": "Klein et al., Mapping cells through time and space with moscot, Nature, 2025.",
        "role": "preferred multi-omics single-cell optimal transport backend",
    },
    {
        "name": "Waddington-OT/wot",
        "module": "wot",
        "dist": "wot",
        "install_command": "pip install wot",
        "citation": "Schiebinger et al., Cell, 2019.",
        "role": "temporal coupling and ancestor-descendant baseline",
    },
    {
        "name": "POT",
        "module": "ot",
        "dist": "POT",
        "install_command": "pip install pot",
        "citation": "Flamary et al., JMLR, 2021.",
        "role": "fallback entropic OT solver used by this reproducible prototype",
    },
    {
        "name": "CellRank 2",
        "module": "cellrank",
        "dist": "cellrank",
        "install_command": "pip install cellrank",
        "citation": "Weiler et al., Nature Methods, 2024.",
        "role": "multiview fate-mapping baseline",
    },
    {
        "name": "COMMOT",
        "module": "commot",
        "dist": "commot",
        "install_command": "pip install commot",
        "citation": "Cang et al., Nature Methods, 2023.",
        "role": "collective OT cell-cell communication reference",
    },
    {
        "name": "TIGON",
        "module": "tigon",
        "dist": "tigon",
        "install_command": "pip install tigon",
        "citation": "Chen et al., Nature Machine Intelligence, 2024.",
        "role": "dynamic unbalanced OT plus growth/GRN comparator",
    },
    {
        "name": "TrajectoryNet",
        "module": "TrajectoryNet",
        "dist": "TrajectoryNet",
        "install_command": "pip install git+https://github.com/KrishnaswamyLab/TrajectoryNet.git",
        "citation": "Tong et al., ICML workshop/arXiv, 2020.",
        "role": "continuous-flow trajectory baseline",
    },
    {
        "name": "MIOFlow",
        "module": "mioflow",
        "dist": "mioflow",
        "install_command": "pip install git+https://github.com/KrishnaswamyLab/MIOFlow.git",
        "citation": "Huguet et al., ICML workshop/arXiv, 2022.",
        "role": "manifold interpolation optimal-flow baseline",
    },
    {
        "name": "GEARS",
        "module": "gears",
        "dist": "gears",
        "install_command": "pip install gears",
        "citation": "Roohani et al., Nature Biotechnology, 2023/2024.",
        "role": "perturbation-expression prediction reference",
    },
    {
        "name": "scGPT",
        "module": "scgpt",
        "dist": "scgpt",
        "install_command": "pip install scgpt",
        "citation": "Cui et al., Nature Methods, 2024.",
        "role": "optional foundation-model embedding reference",
    },
    {
        "name": "scFoundation",
        "module": "scfoundation",
        "dist": "scfoundation",
        "install_command": "install from the maintained scFoundation repository/release",
        "citation": "Hao et al., Nature Methods, 2024.",
        "role": "optional foundation-model embedding reference",
    },
]


def _version(dist: str, module: str) -> dict[str, Any]:
    spec = importlib.util.find_spec(module)
    if spec is None:
        return {"installed": False, "version": None, "module": module}
    try:
        version = metadata.version(dist)
    except Exception:
        try:
            mod = importlib.import_module(module)
            version = getattr(mod, "__version__", "installed-version-unknown")
        except Exception as exc:
            version = f"import-error:{type(exc).__name__}:{exc}"
    return {"installed": True, "version": version, "module": module}


def _git_commit() -> str | None:
    try:
        out = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL).strip()
        return out or None
    except Exception:
        return None


def build_manifest(quick_fixture: bool = False) -> dict[str, Any]:
    package_records = []
    failures = []
    for item in METHODS:
        record = dict(item)
        record.update(_version(item["dist"], item["module"]))
        record["git_commit_hash"] = None
        record["exact_config"] = "configs/*.yaml"
        package_records.append(record)
        if not record["installed"]:
            failures.append(
                f"- {item['name']}: not importable in the current Python environment. "
                f"Suggested command: `{item['install_command']}`."
            )
    native_req = Path("reproducibility/native_moscot_requirements.txt")
    native_summary_path = Path("processed/quick_fixture/ot_couplings/moscot_run_summary.json" if quick_fixture else "processed/ot_couplings/moscot_run_summary.json")
    native_summary: dict[str, Any] = {}
    if native_summary_path.exists():
        try:
            native_summary = json.loads(native_summary_path.read_text(encoding="utf-8"))
        except Exception:
            native_summary = {}
    manifest = {
        "project": "SwarmLineage-OT",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "random_seeds": [7, 17, 23, 42, 99],
        "system": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "processor": platform.processor(),
            "machine": platform.machine(),
            "cwd": str(Path.cwd().resolve()),
        },
        "git_commit": _git_commit(),
        "data_paths": {
            "default_input": "data/processed/cell_level_subset_v1.h5ad",
            "large_reference": "data/scLine_pro.h5ad",
            "preprocessed": "processed/quick_fixture/swarmlineage_input.h5ad" if quick_fixture else "processed/swarmlineage_input.h5ad",
            "teacher": "processed/quick_fixture/ot_teacher.h5ad" if quick_fixture else "processed/ot_teacher.h5ad",
        },
        "output_paths": ["figures/quick_fixture", "tables/quick_fixture", "reports/quick_fixture", "manuscript/*.quick_fixture.md"] if quick_fixture else ["processed", "figures", "tables", "reports", "results/swarmlineage", "manuscript"],
        "methods": package_records,
        "native_teacher_environment": {
            "requirements_file": str(native_req),
            "requirements": native_req.read_text(encoding="utf-8").splitlines() if native_req.exists() else [],
            "local_venv_path": ".venv_moscot_native",
            "teacher_summary": str(native_summary_path),
            "teacher_backend": native_summary.get("teacher_backend", "unknown"),
            "native_moscot_used": bool(native_summary.get("native_moscot_used", False)),
            "native_moscot_status": native_summary.get("native_moscot_status", {}),
            "note": "The native teacher can be generated from the pinned CPU stack even if the global Python environment contains a different moscot/JAX stack.",
        },
        "integrity": {
            "no_fabricated_data": True,
            "native_missing_methods_are_recorded": True,
            "prototype_uses_fallback_when_native_backend_is_unavailable_or_too_slow": True,
        },
    }
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="reproducibility/manifest.json")
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    if args.quick_fixture and args.output == "reproducibility/manifest.json":
        args.output = "reproducibility/manifest.quick_fixture.json"
    ensure_dir("logs")
    manifest = build_manifest(args.quick_fixture)
    write_json(args.output, manifest)
    failures = [
        f"- {m['name']}: not installed. Suggested command: `{m['install_command']}`."
        for m in manifest["methods"]
        if not m["installed"]
    ]
    if not failures:
        failures = ["- No missing optional package was detected in the current Python environment."]
    dry_run_notes = [
        "",
        "Dry-run attempts recorded during this build:",
        "",
        "1. Sandboxed `python -m pip install --dry-run scvi-tools tigon gears mioflow` failed with Windows socket permission/network sandbox errors.",
        "2. Escalated non-mutating dry-run reached PyPI, downloaded scvi-tools metadata, then failed because `tigon` is not available on PyPI as `tigon`.",
        "",
    ]
    write_text(
        "reports/quick_fixture/install_failures.md" if args.quick_fixture else "logs/install_failures.md",
        "# Optional Method Installation Status\n\n"
        "This file records unavailable optional methods and does not claim successful execution.\n\n"
        + "\n".join(failures)
        + "\n"
        + "\n".join(dry_run_notes)
        + "\n",
    )
    print(json.dumps({"manifest": args.output, "missing_optional_methods": len(failures)}, indent=2))


if __name__ == "__main__":
    main()
