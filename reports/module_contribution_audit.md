# Module Contribution Audit

OT interpolation is the teacher/reference map, not a competitor to beat.

- teacher_fidelity_tier: acceptable
- emergent_law_tier: weak
- mechanistic_usefulness_tier: weak

## Reconstruction Context

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

## Discovery Tiers

| law               | tier       | gate_pass   | strong_gate   |   effect_size |   effect_ci_low |   effect_ci_high |   permutation_p |   permutation_q | negative_control_pass   | seed_stability_pass   | rollout_based   | directly_supervised_or_encoded   | interpretation_level              | table                                   | report                                 | status   |
|:------------------|:-----------|:------------|:--------------|--------------:|----------------:|-----------------:|----------------:|----------------:|:------------------------|:----------------------|:----------------|:---------------------------------|:----------------------------------|:----------------------------------------|:---------------------------------------|:---------|
| diffusion         | acceptable | True        | False         |   0.0265607   |     0.0262808   |      0.0268474   |      0.00990099 |        0.029703 | True                    | True                  | False           | True                             | encoded_control_law_recovery      | tables\diffusion_law_regression.csv     | reports\discovery_diffusion_law.md     | executed |
| birth_death       | fail       | False       | False         |   7.89754e-05 |    -0.000477091 |      0.000584699 |      0.019802   |        0.039604 | True                    | False                 | True            | False                            | unsupported                       | tables\birth_death_law.csv              | reports\discovery_birth_death_law.md   | executed |
| branch_nucleation | strong     | True        | True          |  -0.226249    |    -0.259899    |     -0.193552    |      0.00990099 |        0.029703 | True                    | True                  | True            | False                            | retained_computational_hypothesis | tables\swarm_order_parameters.csv       | reports\discovery_branch_nucleation.md | executed |
| memory_hysteresis | fail       | False       | False         |   0           |     0           |      0           |      1          |        1        | False                   | False                 | True            | False                            | unsupported                       | tables\memory_hysteresis_experiment.csv | reports\discovery_memory_hysteresis.md | executed |
| cci_branch_bias   | fail       | False       | False         |   0           |     0           |      0           |      1          |        1        | False                   | False                 | True            | False                            | unsupported                       | tables\cci_branch_bias.csv              | reports\discovery_cci_branch_bias.md   | executed |
| phase_diagram     | weak       | False       | False         |   1           |     1           |      1           |      1          |        1        | False                   | True                  | True            | False                            | exploratory_sensitivity           | tables\phase_diagram.csv                | reports\discovery_phase_diagram.md     | executed |
