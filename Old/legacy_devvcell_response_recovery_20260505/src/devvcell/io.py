"""I/O helpers shared by DevVCell scripts."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def resolve_project_path(path_like: str | Path) -> Path:
    """Resolve a path relative to the project root unless it is absolute."""

    path = Path(path_like)
    if path.is_absolute():
        return path
    return PROJECT_ROOT / path


def load_json(path_like: str | Path) -> dict[str, Any]:
    path = resolve_project_path(path_like)
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path_like: str | Path, payload: dict[str, Any]) -> None:
    path = resolve_project_path(path_like)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
