# Output Integrity Audit

- Main internal outputs remain under top-level reports/tables.
- External E1 outputs are under `processed/external`, `tables/external*`, `reports/external*`, and `figures/external`.
- External generated h5ad/couplings are not committed because processed data are gitignored.
- Registry-only fallback candidates are not described as validation results.
- Native and fallback teacher backends are reported explicitly.
