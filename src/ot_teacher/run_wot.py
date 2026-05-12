from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path

from src.ot_teacher.run_moscot import run_ot
from src.utils.config import load_config, write_json


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/ot_teacher.yaml")
    args = parser.parse_args()
    cfg = load_config(args.config)
    wot_available = importlib.util.find_spec("wot") is not None
    summary = run_ot(cfg, label="wot")
    summary["native_wot_available"] = bool(wot_available)
    summary["native_wot_used"] = False
    summary["note"] = "WOT-compatible temporal entropic couplings were computed with the same auditable fallback solver to support sensitivity comparison."
    write_json(Path(cfg.get("couplings_dir", "processed/ot_couplings")) / "wot_run_summary.json", summary)
    print(json.dumps({"wot_fallback_pairs": len(summary["pairs"]), "native_wot_available": wot_available}, indent=2))


if __name__ == "__main__":
    main()

