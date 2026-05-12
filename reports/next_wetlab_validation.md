# Next Wet-Lab Validation

Minimal validation should test model-predicted shifts in fate, growth and diffusion rather than only expression.

| priority | perturbation | expected readout | assay | decision criterion |
|---:|---|---|---|---|
| 1 | FGF/FGFR modulation in gastruloid or embryo-derived culture | mesoderm/neural fate ratio, MKI67/TOP2A, branch entropy | perturb-and-profile scRNA-seq time course | direction matches SwarmLineage-OT and negative LR control fails |
| 2 | CXCL12/CXCR4 niche-axis perturbation | migration/dispersion and lineage bias | spatial transcriptomics or barcoded live endpoint | spatial/latent dispersion changes in predicted direction |
| 3 | WNT/FZD perturbation | primitive streak/mesoderm branch probability | short time-course scRNA-seq | branch probability and proliferation marker shift |
| 4 | density titration | density-dependent birth/death hazard | controlled aggregate size series | hazard changes with density after cell-cycle adjustment |
| 5 | lineage barcoding follow-up | ancestor-descendant calibration | CellTag/CRISPR barcode with snapshots | OT teacher and simulator improve over interpolation |
