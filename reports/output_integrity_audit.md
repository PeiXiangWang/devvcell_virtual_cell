# Output Integrity Audit

- L2 raw and processed h5ad files are under ignored data/processed paths.
- Committed artifacts are code, reports and CSV summaries.
- E1 and internal teacher backends remain labelled separately from L2 results.
- L2 reports the actual validation tier and does not promote failed/weak results.
- Outcome-preserving Jindal/Weinreb native inputs and moscot couplings are written under ignored `processed/external_l6_outcome_preserving/` paths and are not intended for git tracking.
- Developmental atlas native inputs and moscot couplings are written under ignored `processed/developmental_atlas/` paths and are not intended for git tracking.
- The downloaded CellRank/Farrell zebrafish AnnData is written under ignored `data/external_developmental/` and is not intended for git tracking; only the source URL, config, tables and reports are committed.
- E5 native epsilon sensitivity outputs are written under ignored `processed/developmental_atlas_sensitivity/`; the committed table records backend, pair count and plan shapes without committing transport matrices.
- Committed outcome-preserving artifacts are limited to code, small configs, CSV summaries, reports and one small figure.
- `python -m src.train.evaluate --config configs/train.yaml` now reads the latest clone audit table and preserves the outcome-preserving clone boundary in final retained documents.
- `python -m src.train.evaluate --config configs/train.yaml` also reads the developmental branch-window atlas summary and preserves its weak-tier boundary in final retained documents.
- Bash/WSL quick-fixture validation was attempted but unavailable in this Windows environment; the PowerShell quick-fixture path passed.
