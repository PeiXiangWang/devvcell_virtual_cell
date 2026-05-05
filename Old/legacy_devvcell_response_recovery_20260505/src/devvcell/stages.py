"""Stage parsing helpers for Theiler-stage developmental analyses."""

from __future__ import annotations

import re


_STAGE_RE = re.compile(r"(\d+)")


def stage_number(label: object) -> int:
    """Return the numeric Theiler stage from a label.

    Examples
    --------
    >>> stage_number("Theiler stage 15")
    15
    >>> stage_number("TS15")
    15
    """

    match = _STAGE_RE.findall(str(label))
    if not match:
        raise ValueError(f"Cannot parse Theiler stage number from {label!r}")
    return int(match[-1])


def canonical_stage(label: object) -> str:
    """Return a normalized stage label."""

    return f"Theiler stage {stage_number(label)}"
