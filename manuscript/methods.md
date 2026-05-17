# Methods

The final analysis uses a pre-registered branch-window detector based on lineage-separation contraction, post-event divergence, local velocity alignment, fate entropy, branch imbalance and local density. Native moscot is attempted first for usable time-series datasets; fallback or unusable status is recorded explicitly.

For the final sprint, GSE154572 was downloaded from GEO and converted to AnnData using WT embryoid-body cells, four ordered time points and unsupervised cluster labels as a proxy taxonomy. Because curated lineage labels are absent, this analysis is capped at weak support. STDS0000074/GSE123187 was verified through STOMICS and a public h5ad file was inspected; it did not provide a cell-level, multi-stage, cell-type annotated matrix suitable for the branch-window detector in the inspected form.

## GRN / Regulon Audit

GRN analysis was run on the existing internal, E1, E3, E5, GSE154572 and E2 AnnData objects. The pipeline first checks whether pySCENIC/GRNBoost2-style and CellOracle-style packages are importable; they were not claimed as successfully run unless importable. The executable analysis reuses PathwayFinder's local OmniPath CollecTRI and DoRothEA gene-symbol TF-target priors, which had been downloaded and audited in the PathwayFinder knowledge-graph workflow. For each matched TF, prior targets are used to compute signed regulon activity; TF-target co-expression is used only when no local prior targets overlap the dataset. Branch-window regulatory order parameters include regulon condensation, regulon alignment, regulon divergence, TF switching, branch-specific regulon entropy, susceptibility and OT/time-lag consistency proxies. Negative controls include time shuffling, lineage-label shuffling, regulon-activity shuffling and random TF programs. Perturbation rows are in silico activity-rescaling sensitivity probes only, not experimental validation.

## Breakthrough Sprint Synthesis

The breakthrough synthesis pre-specified ten candidate upgrades and evaluated each against internal native, E1, developmental-atlas, clone-stress, topological-neighbor and GRN audit outputs. Each direction was assigned hypothesis, success criteria, failure criteria, tier and manuscript role. The final selected story was chosen by evidence tier rather than by architectural ambition: branch-window taxonomy, agent rollout interpretability and failure-boundary auditing were retained; GRN mechanism, clone prediction, topological specificity and swarm-required causality were not retained.
