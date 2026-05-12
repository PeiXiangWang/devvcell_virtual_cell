from __future__ import annotations

import argparse
import importlib
from pathlib import Path

from src.utils.config import ensure_dir, load_config, write_text


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/baselines.yaml")
    args = parser.parse_args()
    cfg = load_config(args.config)
    ensure_dir(cfg.get("output_dir", "results/baselines"))
    available = importlib.util.find_spec("tigon") is not None
    text = "# TIGON Baseline\n\n"
    text += f"- package importable: {available}\n"
    text += "- status: not executed natively in the quick run; TIGON is treated as a required external comparator before any high-impact claim.\n"
    write_text(Path(cfg.get("output_dir", "results/baselines")) / "tigon_status.md", text)
    print({"tigon_available": available})


if __name__ == "__main__":
    main()

