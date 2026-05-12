# SwarmLineage-OT: lineage-supervised swarm virtual cells from optimal transport pseudo-lineages

## Abstract

Optimal transport can infer pseudo-lineage maps from destructive single-cell snapshots, but it does not by itself produce an executable virtual cell population. We introduce SwarmLineage-OT, a prototype lineage-supervised swarm virtual-cell model in which finite cellular agents use local rules, density-dependent birth-death, adaptive diffusion and cell-cell communication to generate developmental trajectories constrained by OT-inferred couplings.

## Results

We audited local single-cell resources, built a real-data stage-based OT teacher, trained a minimal finite-agent simulator, and evaluated held-out stage reconstruction across eleven ablations and negative controls. The current quick-run result is a research prototype, not a final Nature-level result.

## Discussion

The central contribution is the conversion of OT pseudo-lineage into executable finite-agent supervision. The main limitation is that validation remains computational and relies on fallback OT couplings in the quick run.
