# Ablation Interpretation

This report is generated directly from `tables/ablation_metrics.csv`; it is not a claim of publishable superiority.

## Mean Metrics

| ('model', '')                         |   ('sinkhorn', 'mean') |   ('sinkhorn', 'std') |   ('mmd_rbf', 'mean') |   ('mmd_rbf', 'std') |   ('energy', 'mean') |   ('energy', 'std') |   ('celltype_composition_rmse', 'mean') |   ('celltype_composition_rmse', 'std') |
|:--------------------------------------|-----------------------:|----------------------:|----------------------:|---------------------:|---------------------:|--------------------:|----------------------------------------:|---------------------------------------:|
| M0_linear_label_interpolation         |               0.281661 |            0.00953524 |            0.00822039 |          0.00061027  |             0.292374 |           0.0112253 |                              0.00396321 |                             0.00253493 |
| M0b_ot_interpolation                  |               0.255076 |            0.0111819  |            0.00878499 |          0.00047848  |             0.335883 |           0.0142804 |                              0.00348103 |                             0.0019157  |
| M10_shuffled_time_ot                  |               0.31024  |            0.0146304  |            0.00516795 |          0.000847604 |             0.204379 |           0.0306884 |                              0.018496   |                             0.00861311 |
| M11_random_lr_labels                  |               0.289617 |            0.0125825  |            0.0172297  |          0.00226102  |             0.628864 |           0.0428832 |                              0.0233767  |                             0.00445816 |
| M1_intrinsic_neural                   |               0.309506 |            0.00982942 |            0.0276649  |          0.00152121  |             0.915083 |           0.013263  |                              0.00387904 |                             0.00193728 |
| M2_ot_teacher_force                   |               0.274306 |            0.00931908 |            0.0167648  |          0.000799382 |             0.568815 |           0.0139597 |                              0.00386852 |                             0.0021581  |
| M3_ot_birth_death                     |               0.288529 |            0.0110321  |            0.0164293  |          0.00251527  |             0.590722 |           0.0443668 |                              0.0168219  |                             0.00997471 |
| M4_ot_adaptive_diffusion              |               0.274514 |            0.00918233 |            0.0167285  |          0.000816995 |             0.566489 |           0.0140748 |                              0.00420092 |                             0.00233245 |
| M5_ot_swarm                           |               0.273951 |            0.00930118 |            0.0178781  |          0.000872036 |             0.613508 |           0.0144048 |                              0.00448128 |                             0.00237237 |
| M6_ot_swarm_birth_death               |               0.292737 |            0.00923603 |            0.017475   |          0.00240024  |             0.630969 |           0.0384701 |                              0.0209927  |                             0.00931139 |
| M7_ot_swarm_birth_death_diffusion     |               0.284228 |            0.0151806  |            0.0183164  |          0.00343481  |             0.632902 |           0.0727929 |                              0.0200099  |                             0.00872083 |
| M8_ot_swarm_birth_death_diffusion_cci |               0.284235 |            0.0151803  |            0.0183179  |          0.00343563  |             0.632942 |           0.0728207 |                              0.0200099  |                             0.00872083 |
| M9_full_memory                        |               0.28393  |            0.0152408  |            0.0183946  |          0.00344891  |             0.637903 |           0.0729498 |                              0.0200099  |                             0.00872083 |

## Primary Paired Tests

| metric                    | baseline            | challenger     |   n |   baseline_mean |   challenger_mean |   effect_baseline_minus_challenger |   effect_ci_low |   effect_ci_high |   p_value |   q_value |
|:--------------------------|:--------------------|:---------------|----:|----------------:|------------------:|-----------------------------------:|----------------:|-----------------:|----------:|----------:|
| sinkhorn                  | M2_ot_teacher_force | M9_full_memory |   5 |      0.274306   |         0.28393   |                        -0.00962397 |     -0.0262105  |      0.00642452  |   0.90625 |         1 |
| mmd_rbf                   | M2_ot_teacher_force | M9_full_memory |   5 |      0.0167648  |         0.0183946 |                        -0.00162977 |     -0.00467704 |      0.000949307 |   0.84375 |         1 |
| energy                    | M2_ot_teacher_force | M9_full_memory |   5 |      0.568815   |         0.637903  |                        -0.0690878  |     -0.131981   |     -0.0212908   |   1       |         1 |
| celltype_composition_rmse | M2_ot_teacher_force | M9_full_memory |   5 |      0.00386852 |         0.0200099 |                        -0.0161414  |     -0.0236141  |     -0.00847339  |   1       |         1 |

## Current Interpretation

Strongest non-reference baseline for paired diagnostic: `M2_ot_teacher_force`. `M0b_ot_interpolation` is retained as the OT teacher/reference interpolation, not as a competitor that the finite-agent model must beat. Modules are evaluated by teacher fidelity plus whether they provide stable mechanistic diagnostics and emergent-law signals.
