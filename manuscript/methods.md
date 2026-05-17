# Methods

The final analysis uses a pre-registered branch-window detector based on lineage-separation contraction, post-event divergence, local velocity alignment, fate entropy, branch imbalance and local density. Native moscot is attempted first for usable time-series datasets; fallback or unusable status is recorded explicitly.

For the final sprint, GSE154572 was downloaded from GEO and converted to AnnData using WT embryoid-body cells, four ordered time points and unsupervised cluster labels as a proxy taxonomy. Because curated lineage labels are absent, this analysis is capped at weak support. STDS0000074/GSE123187 was verified through STOMICS and a public h5ad file was inspected; it did not provide a cell-level, multi-stage, cell-type annotated matrix suitable for the branch-window detector in the inspected form.

## GRN / Regulon Audit

GRN analysis was run on the existing internal, E1, E3, E5, GSE154572 and E2 AnnData objects. The pipeline first checks whether pySCENIC/GRNBoost2-style and CellOracle-style packages are importable; they were not claimed as successfully run unless importable. The executable fallback uses a curated developmental TF list, matches TF symbols to each dataset, infers positive TF-target co-expression modules, and computes per-cell fallback regulon activities. Branch-window regulatory order parameters include regulon condensation, regulon alignment, regulon divergence, TF switching, branch-specific regulon entropy, susceptibility and OT/time-lag consistency proxies. Negative controls include time shuffling, lineage-label shuffling, regulon-activity shuffling and random TF programs. Perturbation rows are in silico activity-rescaling sensitivity probes only, not experimental validation.
