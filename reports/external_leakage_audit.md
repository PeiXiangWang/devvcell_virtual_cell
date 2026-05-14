# External Leakage Audit

- External preprocessing was fit only on the selected external MouseGastrulationData component.
- Internal SwarmLineage-OT data, internal moscot teacher and internal learned embeddings were not used to fit the external PCA.
- `split_role` is set to train for all external cells because E1 tests branch-nucleation order parameters, not held-out prediction.
- No lineage barcodes were present in this external component; lineage validation is not claimed.
