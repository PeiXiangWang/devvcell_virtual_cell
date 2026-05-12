# SwarmLineage-OT: lineage-supervised swarm virtual cells from optimal transport pseudo-lineages

## Abstract

Optimal transport can infer pseudo-lineage maps from destructive single-cell snapshots, but it does not by itself produce an executable virtual cell population or expose microscopic developmental control laws. We introduce SwarmLineage-OT, a prototype lineage-supervised swarm virtual-cell model in which finite cellular agents use local rules, density-dependent birth-death, adaptive diffusion, memory fields and cell-cell communication to realize OT-inferred couplings and probe emergent developmental laws.

## Results

We audited local single-cell resources, built a stage-based OT teacher, trained a finite-agent simulator, and evaluated the model with three gates: teacher fidelity, emergent-law robustness and mechanistic usefulness. OT interpolation is treated as an oracle-like teacher/reference interpolation, not as a competitor that the agent model must outperform.

Discovery analyses estimate fate-uncertainty-driven diffusion, density-dependent birth/death, branch nucleation order parameters, CCI-mediated branch bias, memory-dependent hysteresis and finite-agent phase regimes.

## Discussion

The central contribution is the conversion of OT pseudo-lineage into executable finite-agent supervision and mechanistic law discovery. The main limitation is that validation remains computational and strong biological claims require native moscot/WOT or external teacher validation plus external lineage, spatial or perturbation evidence.
