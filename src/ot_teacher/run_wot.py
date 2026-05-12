from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path

from src.utils.config import load_config, write_json


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/ot_teacher.yaml")
    parser.add_argument("--quick-fixture", action="store_true")
    args = parser.parse_args()
    cfg = load_config(args.config)
    if args.quick_fixture:
        cfg = dict(cfg)
        cfg["couplings_dir"] = "processed/quick_fixture/ot_couplings"
    wot_available = importlib.util.find_spec("wot") is not None
    summary = {
        "native_wot_available": bool(wot_available),
        "native_wot_used": False,
        "teacher_backend": "native_wot" if False else "status_only",
        "pairs": [],
        "note": "WOT is not executed by reusing the toy Sinkhorn teacher. Native WOT execution must be implemented separately before WOT can be reported as a compared method.",
    }
    write_json(Path(cfg.get("couplings_dir", "processed/ot_couplings")) / "wot_run_summary.json", summary)
    print(json.dumps({"wot_status_only": True, "native_wot_available": wot_available}, indent=2))


if __name__ == "__main__":
    main()
