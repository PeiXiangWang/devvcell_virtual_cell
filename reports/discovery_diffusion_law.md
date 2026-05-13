# Discovery Diffusion Law

Diffusion is hardened with standard regression, entropy-shuffle negative control, no-entropy proxy analysis and seed-wise coefficient stability.

The current training objective directly includes a sigma-to-entropy calibration term. Therefore entropy-associated learned sigma is classified as encoded control-law recovery unless independently supported by rollout displacement.

## Tier

- tier: acceptable
- effect_size: 0.0265607
- 95% CI: [0.0262808, 0.0268474]
- permutation_p: 0.00990099
- permutation_q: 0.00990099
- seed_stability_pass: True (sign consistency=1.000)
- negative_control_pass: True
- directly_supervised_or_encoded: True

## Top Regression Terms

| response             | predictor             |       coef |   abs_coef |       r2 |     n | analysis_mode                |
|:---------------------|:----------------------|-----------:|-----------:|---------:|------:|:-----------------------------|
| learned_displacement | fate_entropy          |  1.99224   |  1.99224   | 0.381384 | 40000 | empirical_displacement_proxy |
| learned_displacement | fate_probability_max  |  1.707     |  1.707     | 0.381384 | 40000 | empirical_displacement_proxy |
| learned_displacement | ot_transition_entropy |  0.94625   |  0.94625   | 0.381384 | 40000 | empirical_displacement_proxy |
| learned_displacement | local_density         | -0.441477  |  0.441477  | 0.381384 | 40000 | empirical_displacement_proxy |
| learned_displacement | cci_signal            |  0.25106   |  0.25106   | 0.381384 | 40000 | empirical_displacement_proxy |
| learned_sigma        | fate_entropy          |  0.082005  |  0.082005  | 0.332279 | 40000 | no_entropy_input_proxy       |
| learned_sigma        | fate_probability_max  |  0.0706911 |  0.0706911 | 0.332279 | 40000 | no_entropy_input_proxy       |
| learned_sigma        | ot_transition_entropy |  0.0265607 |  0.0265607 | 0.977212 | 40000 | standard                     |
