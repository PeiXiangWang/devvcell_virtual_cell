"""Markdown reporting helpers."""

from __future__ import annotations

import pandas as pd


def results_report(
    *,
    timeseries: pd.DataFrame,
    spectral_summary: pd.DataFrame,
    reconstruction_summary: pd.DataFrame,
    fingerprint: pd.DataFrame,
) -> str:
    lines = [
        "# DevSpectrum Results Report",
        "",
        "## Inputs",
        "",
        f"- time-series rows: {timeseries.shape[0]}",
        f"- lineages: {timeseries['lineage'].nunique() if 'lineage' in timeseries else 0}",
        f"- modules/features: {timeseries['module_name'].nunique() if 'module_name' in timeseries else 0}",
        "",
        "## Spectral fit summary",
        "",
        spectral_summary.head(20).to_markdown(index=False) if not spectral_summary.empty else "No spectral summary available.",
        "",
        "## Missing-stage reconstruction",
        "",
        reconstruction_summary.to_markdown(index=False) if not reconstruction_summary.empty else "No reconstruction summary available.",
        "",
        "## Perturbation spectral fingerprint",
        "",
        fingerprint.to_markdown(index=False) if not fingerprint.empty else "No perturbation fingerprint available.",
        "",
        "## Limitations",
        "",
        "- GSE212050 provides a short five-point control time-course, so DCT is intentionally low-order.",
        "- Tal1/T chimera data are treated as E8.5 endpoint projections, not full perturbation time spectra.",
        "- Rescue candidates are in silico spectral hypotheses.",
        "",
    ]
    return "\n".join(lines)
