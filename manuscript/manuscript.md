# SwarmLineage-OT: lineage-supervised swarm virtual cells from optimal transport pseudo-lineages

## Abstract

Optimal transport can infer pseudo-lineage maps from destructive single-cell snapshots, but it does not by itself produce an executable virtual cell population or expose microscopic developmental control laws. We introduce SwarmLineage-OT, a prototype lineage-supervised swarm virtual-cell model in which finite cellular agents use local rules, density-dependent birth-death, adaptive diffusion, memory fields and cell-cell communication to realize OT-inferred couplings and probe emergent developmental laws.

## Results

We audited local single-cell resources, built a stage-based OT teacher, trained a finite-agent simulator, and evaluated the model with tiered evidence gates: teacher fidelity, emergent-law robustness and mechanistic usefulness. OT interpolation is treated as an oracle-like teacher/reference interpolation, not as a competitor that the agent model must outperform.

In the current discovery-hardened run, teacher fidelity is acceptable, emergent-law evidence is weak overall, and mechanistic usefulness is weak. Diffusion is retained only as an encoded control-law recovery; branch nucleation and phase regimes remain exploratory rollout-based probes; birth/death, memory hysteresis and CCI branch bias are unsupported in the current evidence table.

## Discussion

The central contribution is the conversion of OT pseudo-lineage into executable finite-agent supervision and a falsifiable mechanism-discovery audit. The current teacher backend is native moscot on a downsampled main AnnData run, which removes the toy-fallback blocker for teacher construction. The main limitation is that validation remains computational, and strong biological claims still require external lineage, spatial, perturbation or wet-lab evidence.
