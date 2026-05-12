# OT Teacher Report

The teacher was built from adjacent developmental-stage entropic OT couplings. Native moscot availability is recorded, but quick execution uses a POT/SciPy fallback plan when native moscot is too slow or unavailable.

- cells: 8000
- latent dimensions: 30
- terminal pseudo-fates: neural, erythroid, mesoderm_muscle
- mean transition entropy: 0.6197
- mean growth proxy: 1.0000
- teacher reliability score, heuristic not a publication claim: 0.6900

## Caveats

- These couplings are OT-inferred pseudo-lineage, not true lineage tracing.
- `stage_num`/Theiler stage is treated as real ordered developmental stage for this dataset; it is still coarser than dense experimental time.
- Native moscot/WOT baselines should be rerun without quick fallback before any high-impact claim.

## WOT-Style Sensitivity

|   source_time |   target_time |   mass_js |
|--------------:|--------------:|----------:|
|            12 |            13 |         0 |
|            13 |            14 |         0 |
|            14 |            15 |         0 |
|            15 |            16 |         0 |
|            16 |            17 |         0 |
|            17 |            18 |         0 |
|            18 |            19 |         0 |
