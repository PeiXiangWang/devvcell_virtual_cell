# Module Contribution Audit

This audit no longer treats OT interpolation as a method that the agent model must beat. OT interpolation is the teacher/reference map.

- primary agent: `M9_full_memory`
- `M9_full_memory` fidelity: relative_sinkhorn=1.113, relative_mmd=2.094, branch_fate_error=0.0200.
- mechanistic_usefulness_gate: True

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

## Emergent-Law Contributions

| law               | gate   | table                                   | status   |
|:------------------|:-------|:----------------------------------------|:---------|
| diffusion         | True   | tables\diffusion_law_regression.csv     | executed |
| birth_death       | True   | tables\birth_death_law.csv              | executed |
| branch_nucleation | True   | tables\swarm_order_parameters.csv       | executed |
| memory_hysteresis | True   | tables\memory_hysteresis_experiment.csv | executed |
| cci_branch_bias   | True   | tables\cci_branch_bias.csv              | executed |
| phase_diagram     | True   | tables\phase_diagram.csv                | executed |

## Interpretation

- `M0b_ot_interpolation` can remain the lowest-error reconstruction because it is the reference interpolation derived from the teacher.
- Swarm, birth/death, diffusion, CCI and memory are evaluated by fidelity plus whether they expose robust order parameters, hazards, phase regimes or perturbation sensitivities.
- Modules that degrade fidelity without yielding stable laws should remain exploratory or be disabled in retained biological claims.
