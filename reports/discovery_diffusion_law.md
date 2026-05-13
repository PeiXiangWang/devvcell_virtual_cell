# Discovery Diffusion Law

Diffusion is hardened with standard regression, entropy-shuffle negative control, no-entropy proxy analysis and seed-wise coefficient stability.

The current training objective directly includes a sigma-to-entropy calibration term. Therefore entropy-associated learned sigma is classified as encoded control-law recovery unless independently supported by rollout displacement.

## Tier

- tier: acceptable
- effect_size: 0.0288975
- 95% CI: [0.0286935, 0.0290633]
- permutation_p: 0.00990099
- permutation_q: 0.00990099
- seed_stability_pass: True (sign consistency=1.000)
- negative_control_pass: True
- directly_supervised_or_encoded: True

## Top Regression Terms

| response             | predictor             |       coef |   abs_coef |       r2 |     n | analysis_mode                |
|:---------------------|:----------------------|-----------:|-----------:|---------:|------:|:-----------------------------|
| learned_displacement | ot_transition_entropy |  1.3229    |  1.3229    | 0.521076 | 40000 | empirical_displacement_proxy |
| learned_displacement | local_density         | -0.866772  |  0.866772  | 0.521076 | 40000 | empirical_displacement_proxy |
| learned_displacement | fate_entropy          |  0.527551  |  0.527551  | 0.521076 | 40000 | empirical_displacement_proxy |
| learned_displacement | fate_probability_max  |  0.284236  |  0.284236  | 0.521076 | 40000 | empirical_displacement_proxy |
| learned_displacement | cci_signal            |  0.197473  |  0.197473  | 0.521076 | 40000 | empirical_displacement_proxy |
| learned_sigma        | fate_entropy          |  0.0728326 |  0.0728326 | 0.285082 | 40000 | no_entropy_input_proxy       |
| learned_sigma        | fate_probability_max  |  0.0612388 |  0.0612388 | 0.285082 | 40000 | no_entropy_input_proxy       |
| learned_sigma        | ot_transition_entropy |  0.0288975 |  0.0288975 | 0.980809 | 40000 | standard                     |
