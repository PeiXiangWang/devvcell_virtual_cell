# Discovery Diffusion Law

The analysis regresses learned diffusion scale against OT transition entropy, local density, fate commitment, cell-cycle score and CCI signal.

## Top Effects

| response      | predictor             |        coef |   abs_coef |       r2 |    n |
|:--------------|:----------------------|------------:|-----------:|---------:|-----:|
| learned_sigma | ot_transition_entropy |  0.0286115  | 0.0286115  | 0.983058 | 8000 |
| learned_sigma | fate_probability_max  | -0.0034916  | 0.0034916  | 0.983058 | 8000 |
| learned_sigma | local_density         |  0.00101922 | 0.00101922 | 0.983058 | 8000 |

- diffusion_law_gate: True
- Interpretation is mechanistic only if the signal is stable across seeds and not driven by teacher fallback artifacts.
