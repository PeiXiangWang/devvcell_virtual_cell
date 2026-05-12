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
    available = importlib.util.find_spec("mioflow") is not None
    text = "# MIOFlow Baseline\n\n"
    text += f"- package importable: {available}\n"
    text += "- status: not executed natively in this quick run; it remains a required comparator for a publishable benchmark.\n"
    write_text(Path(cfg.get("output_dir", "results/baselines")) / "mioflow_status.md", text)
    print({"mioflow_available": available})


if __name__ == "__main__":
    main()

