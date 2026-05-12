from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def load_config(path: str | Path | None) -> dict[str, Any]:
    if path is None:
        return {}
    config_path = Path(path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file does not exist: {config_path}")
    text = config_path.read_text(encoding="utf-8")
    if config_path.suffix.lower() in {".json"}:
        return json.loads(text)
    try:
        import yaml

        data = yaml.safe_load(text)
        return data or {}
    except Exception as exc:
        raise RuntimeError(f"Failed to parse config {config_path}: {exc}") from exc


def ensure_parent(path: str | Path) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    return out


def ensure_dir(path: str | Path) -> Path:
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def write_json(path: str | Path, payload: Any) -> None:
    out = ensure_parent(path)
    def default(obj: Any) -> Any:
        try:
            import numpy as np

            if isinstance(obj, np.integer):
                return int(obj)
            if isinstance(obj, np.floating):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
        except Exception:
            pass
        return str(obj)

    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=default), encoding="utf-8")


def write_text(path: str | Path, text: str) -> None:
    out = ensure_parent(path)
    out.write_text(text, encoding="utf-8")


def as_path(path: str | Path) -> Path:
    return Path(path).expanduser().resolve()
