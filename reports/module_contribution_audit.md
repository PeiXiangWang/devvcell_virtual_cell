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

## Discovery Tiers

| law               | tier       | gate_pass   | strong_gate   |   effect_size |   effect_ci_low |   effect_ci_high |   permutation_p |   permutation_q | negative_control_pass   | seed_stability_pass   | rollout_based   | directly_supervised_or_encoded   | interpretation_level         | table                                   | report                                 | status   |
|:------------------|:-----------|:------------|:--------------|--------------:|----------------:|-----------------:|----------------:|----------------:|:------------------------|:----------------------|:----------------|:---------------------------------|:-----------------------------|:----------------------------------------|:---------------------------------------|:---------|
| diffusion         | acceptable | True        | False         |   0.0288975   |     0.0286935   |      0.0290633   |      0.00990099 |       0.0594059 | True                    | True                  | False           | True                             | encoded_control_law_recovery | tables\diffusion_law_regression.csv     | reports\discovery_diffusion_law.md     | executed |
| birth_death       | fail       | False       | False         |   9.33802e-05 |    -0.000395111 |      0.000551136 |      0.019802   |       0.0594059 | True                    | False                 | True            | False                            | unsupported                  | tables\birth_death_law.csv              | reports\discovery_birth_death_law.md   | executed |
| branch_nucleation | weak       | False       | False         |  -0.0572253   |    -0.0830656   |     -0.0253689   |      0.148515   |       0.29703   | False                   | False                 | True            | False                            | exploratory_sensitivity      | tables\swarm_order_parameters.csv       | reports\discovery_branch_nucleation.md | executed |
| memory_hysteresis | fail       | False       | False         |   0           |     0           |      0           |      1          |       1         | False                   | False                 | True            | False                            | unsupported                  | tables\memory_hysteresis_experiment.csv | reports\discovery_memory_hysteresis.md | executed |
| cci_branch_bias   | fail       | False       | False         |   0           |     0           |      0           |      1          |       1         | False                   | False                 | True            | False                            | unsupported                  | tables\cci_branch_bias.csv              | reports\discovery_cci_branch_bias.md   | executed |
| phase_diagram     | weak       | False       | False         |   1           |     1           |      1           |      1          |       1         | False                   | True                  | True            | False                            | exploratory_sensitivity      | tables\phase_diagram.csv                | reports\discovery_phase_diagram.md     | executed |
