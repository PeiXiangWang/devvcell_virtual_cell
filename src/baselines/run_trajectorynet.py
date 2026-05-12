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
    available = importlib.util.find_spec("TrajectoryNet") is not None
    text = "# TrajectoryNet Baseline\n\n"
    text += f"- package importable: {available}\n"
    text += "- status: not executed natively in this quick run; use the recorded environment and config to add it as a strong continuous-flow baseline.\n"
    write_text(Path(cfg.get("output_dir", "results/baselines")) / "trajectorynet_status.md", text)
    print({"trajectorynet_available": available})


if __name__ == "__main__":
    main()

