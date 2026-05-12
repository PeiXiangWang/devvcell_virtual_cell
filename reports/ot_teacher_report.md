# OT Teacher Report

The teacher was built from adjacent developmental-stage entropic OT couplings. If native moscot is unavailable, couplings are explicitly labelled toy_sinkhorn_fallback and cannot support high-level moscot claims.

- cells: 8000
- latent dimensions: 30
- terminal pseudo-fates: neural, erythroid, mesoderm_muscle
- mean transition entropy: 0.6218
- mean growth proxy: 1.0000
- teacher backend: toy_sinkhorn_fallback
- teacher reliability score, heuristic not a publication claim: 0.6889

## Caveats

- These couplings are OT-inferred pseudo-lineage, not true lineage tracing.
- `stage_num`/Theiler stage is treated as real ordered developmental stage for this dataset; it is still coarser than dense experimental time.
- Native moscot/WOT baselines should be rerun without quick fallback before any high-impact claim.
