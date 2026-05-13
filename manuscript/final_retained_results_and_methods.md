# Final Retained Results and Methods

## Central Claim

OT gives the developmental map; SwarmLineage-OT learns microscopic finite-agent rules that realize the map and reveal emergent developmental laws.

`M0b_ot_interpolation` is an oracle-like OT teacher/reference interpolation. The finite-agent model is evaluated by teacher fidelity, emergent-law robustness and mechanistic usefulness, not by beating the OT reference.

## Tier Summary

| teacher_fidelity_tier   | emergent_law_tier   | mechanistic_usefulness_tier   | mechanistic_gate_pass   | strong_gate   |   laws_at_least_acceptable |   laws_strong | native_or_external_teacher_validation   |
|:------------------------|:--------------------|:------------------------------|:------------------------|:--------------|---------------------------:|--------------:|:----------------------------------------|
| acceptable              | weak                | weak                          | False                   | False         |                          2 |             1 | True                                    |

## Retained Computational Hypotheses

| law               | tier       | interpretation_level              | rollout_based   | directly_supervised_or_encoded   |
|:------------------|:-----------|:----------------------------------|:----------------|:---------------------------------|
| diffusion         | acceptable | encoded_control_law_recovery      | False           | True                             |
| branch_nucleation | strong     | retained_computational_hypothesis | True            | False                            |

## Exploratory / Demonstration Only

| law           | tier   | interpretation_level    | rollout_based   | directly_supervised_or_encoded   |
|:--------------|:-------|:------------------------|:----------------|:---------------------------------|
| phase_diagram | weak   | exploratory_sensitivity | True            | False                            |

## Unsupported Claims

| law               | tier   | status   |
|:------------------|:-------|:---------|
| birth_death       | fail   | executed |
| memory_hysteresis | fail   | executed |
| cci_branch_bias   | fail   | executed |

## Core Metrics

| model                                 |   sinkhorn |    mmd_rbf |   energy |   celltype_composition_rmse |
|:--------------------------------------|-----------:|-----------:|---------:|----------------------------:|
| M0_linear_label_interpolation         |   0.281661 | 0.00822039 | 0.292374 |                  0.00396321 |
| M0b_ot_interpolation                  |   0.255076 | 0.00878499 | 0.335883 |                  0.00348103 |
| M10_shuffled_time_ot                  |   0.304366 | 0.00526067 | 0.213518 |                  0.014373   |
| M11_random_lr_labels                  |   0.300051 | 0.0135742  | 0.450108 |                  0.0202327  |
| M1_intrinsic_neural                   |   0.309283 | 0.0277509  | 0.917    |                  0.00420263 |
| M2_ot_teacher_force                   |   0.29088  | 0.0120656  | 0.38817  |                  0.0039636  |
| M3_ot_birth_death                     |   0.307424 | 0.0121     | 0.40647  |                  0.0210994  |
| M4_ot_adaptive_diffusion              |   0.291023 | 0.0120433  | 0.387111 |                  0.00384322 |
| M5_ot_swarm                           |   0.288679 | 0.013716   | 0.439041 |                  0.00402532 |
| M6_ot_swarm_birth_death               |   0.306787 | 0.0132087  | 0.451684 |                  0.0224666  |
| M7_ot_swarm_birth_death_diffusion     |   0.307358 | 0.014038   | 0.46045  |                  0.0219325  |
| M8_ot_swarm_birth_death_diffusion_cci |   0.307362 | 0.0140396  | 0.460488 |                  0.0219325  |
| M9_full_memory                        |   0.307017 | 0.0142543  | 0.468093 |                  0.0217414  |

## Limitations

- Current results are computational hypotheses.
- Some laws are encoded control-law recoveries and must not be written as independent biological discoveries.
- Native moscot teacher extraction removes the toy-fallback blocker for teacher construction, but not the need for external validation.
- No wet-lab validation or causal mechanism is claimed.
