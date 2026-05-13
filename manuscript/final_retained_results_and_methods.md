# Final Retained Results and Methods

## Central Claim

OT gives the developmental map; SwarmLineage-OT learns microscopic finite-agent rules that realize the map and reveal emergent developmental laws.

`M0b_ot_interpolation` is an oracle-like OT teacher/reference interpolation. The finite-agent model is evaluated by teacher fidelity, emergent-law robustness and mechanistic usefulness, not by beating the OT reference.

## Tier Summary

| teacher_fidelity_tier   | emergent_law_tier   | mechanistic_usefulness_tier   | mechanistic_gate_pass   | strong_gate   |   laws_at_least_acceptable |   laws_strong | native_or_external_teacher_validation   |
|:------------------------|:--------------------|:------------------------------|:------------------------|:--------------|---------------------------:|--------------:|:----------------------------------------|
| acceptable              | weak                | weak                          | False                   | False         |                          1 |             0 | False                                   |

## Retained Computational Hypotheses

| law       | tier       | interpretation_level         | rollout_based   | directly_supervised_or_encoded   |
|:----------|:-----------|:-----------------------------|:----------------|:---------------------------------|
| diffusion | acceptable | encoded_control_law_recovery | False           | True                             |

## Exploratory / Demonstration Only

| law               | tier   | interpretation_level    | rollout_based   | directly_supervised_or_encoded   |
|:------------------|:-------|:------------------------|:----------------|:---------------------------------|
| branch_nucleation | weak   | exploratory_sensitivity | True            | False                            |
| phase_diagram     | weak   | exploratory_sensitivity | True            | False                            |

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
| M10_shuffled_time_ot                  |   0.31024  | 0.00516795 | 0.204379 |                  0.018496   |
| M11_random_lr_labels                  |   0.289617 | 0.0172297  | 0.628864 |                  0.0233767  |
| M1_intrinsic_neural                   |   0.309506 | 0.0276649  | 0.915083 |                  0.00387904 |
| M2_ot_teacher_force                   |   0.274306 | 0.0167648  | 0.568815 |                  0.00386852 |
| M3_ot_birth_death                     |   0.288529 | 0.0164293  | 0.590722 |                  0.0168219  |
| M4_ot_adaptive_diffusion              |   0.274514 | 0.0167285  | 0.566489 |                  0.00420092 |
| M5_ot_swarm                           |   0.273951 | 0.0178781  | 0.613508 |                  0.00448128 |
| M6_ot_swarm_birth_death               |   0.292737 | 0.017475   | 0.630969 |                  0.0209927  |
| M7_ot_swarm_birth_death_diffusion     |   0.284228 | 0.0183164  | 0.632902 |                  0.0200099  |
| M8_ot_swarm_birth_death_diffusion_cci |   0.284235 | 0.0183179  | 0.632942 |                  0.0200099  |
| M9_full_memory                        |   0.28393  | 0.0183946  | 0.637903 |                  0.0200099  |

## Limitations

- Current results are computational hypotheses.
- Some laws are encoded control-law recoveries and must not be written as independent biological discoveries.
- toy_sinkhorn_fallback is not native moscot/WOT.
- No wet-lab validation or causal mechanism is claimed.
