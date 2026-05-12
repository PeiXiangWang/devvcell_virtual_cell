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
    available = importlib.util.find_spec("cellrank") is not None
    text = "# CellRank 2 Baseline\n\n"
    text += f"- package importable: {available}\n"
    text += "- status: availability recorded; full CellRank2 estimator execution is deferred because the quick prototype uses OT-teacher couplings and latent held-out metrics.\n"
    write_text(Path(cfg.get("output_dir", "results/baselines")) / "cellrank_status.md", text)
    print({"cellrank_available": available})


if __name__ == "__main__":
    main()

