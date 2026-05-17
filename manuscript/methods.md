# Methods

The final analysis uses a pre-registered branch-window detector based on lineage-separation contraction, post-event divergence, local velocity alignment, fate entropy, branch imbalance and local density. Native moscot is attempted first for usable time-series datasets; fallback or unusable status is recorded explicitly.

For the final sprint, GSE154572 was downloaded from GEO and converted to AnnData using WT embryoid-body cells, four ordered time points and unsupervised cluster labels as a proxy taxonomy. Because curated lineage labels are absent, this analysis is capped at weak support. STDS0000074/GSE123187 was verified through STOMICS and a public h5ad file was inspected; it did not provide a cell-level, multi-stage, cell-type annotated matrix suitable for the branch-window detector in the inspected form.
