# Ablation Interpretation

This report is generated directly from `tables/ablation_metrics.csv`; it is not a claim of publishable superiority.

## Mean Metrics

| ('model', '')                         |   ('sinkhorn', 'mean') |   ('sinkhorn', 'std') |   ('mmd_rbf', 'mean') |   ('mmd_rbf', 'std') |   ('energy', 'mean') |   ('energy', 'std') |   ('celltype_composition_rmse', 'mean') |   ('celltype_composition_rmse', 'std') |
|:--------------------------------------|-----------------------:|----------------------:|----------------------:|---------------------:|---------------------:|--------------------:|----------------------------------------:|---------------------------------------:|
| M0_ot_interpolation                   |               0.28892  |            0.00936219 |            0.0097388  |          0.00075215  |             0.345831 |           0.0115977 |                              0.00381493 |                             0.00262295 |
| M10_shuffled_time_ot                  |               0.316176 |            0.00747013 |            0.00513459 |          0.000904152 |             0.218758 |           0.0329254 |                              0.0108743  |                             0.00222133 |
| M11_random_lr_labels                  |               0.301723 |            0.0104156  |            0.0277966  |          0.00378125  |             0.94278  |           0.0765006 |                              0.00983581 |                             0.00551719 |
| M1_intrinsic_neural                   |               0.29825  |            0.00876728 |            0.0243012  |          0.00128358  |             0.809875 |           0.0246339 |                              0.00408543 |                             0.00133439 |
| M2_ot_teacher_force                   |               0.29825  |            0.00876728 |            0.0243012  |          0.00128358  |             0.809875 |           0.0246339 |                              0.00408543 |                             0.00133439 |
| M3_ot_birth_death                     |               0.294263 |            0.00867486 |            0.0257157  |          0.00300539  |             0.849029 |           0.0612029 |                              0.00905709 |                             0.00460614 |
| M4_ot_adaptive_diffusion              |               0.298341 |            0.00868134 |            0.024264   |          0.00127684  |             0.807488 |           0.0248584 |                              0.00408543 |                             0.00133439 |
| M5_ot_swarm                           |               0.301489 |            0.00904874 |            0.0267753  |          0.00143666  |             0.909166 |           0.0269617 |                              0.00390988 |                             0.00210723 |
| M6_ot_swarm_birth_death               |               0.297042 |            0.00876958 |            0.028327   |          0.00321477  |             0.951595 |           0.0637824 |                              0.00915502 |                             0.00496029 |
| M7_ot_swarm_birth_death_diffusion     |               0.297119 |            0.00882763 |            0.0282685  |          0.00321027  |             0.948511 |           0.0641559 |                              0.00973161 |                             0.00388418 |
| M8_ot_swarm_birth_death_diffusion_cci |               0.297228 |            0.00881111 |            0.0283344  |          0.0032256   |             0.950664 |           0.0645807 |                              0.00973161 |                             0.00388418 |
| M9_full_pheromone                     |               0.302977 |            0.0090176  |            0.0309501  |          0.00344266  |             1.05423  |           0.06745   |                              0.0096355  |                             0.00391603 |

## Primary Paired Tests

| metric                    | baseline            | challenger        |   n |   baseline_mean |   challenger_mean |   effect_baseline_minus_challenger |   effect_ci_low |   effect_ci_high |   p_value |   q_value |
|:--------------------------|:--------------------|:------------------|----:|----------------:|------------------:|-----------------------------------:|----------------:|-----------------:|----------:|----------:|
| sinkhorn                  | M0_ot_interpolation | M9_full_pheromone |   5 |      0.28892    |         0.302977  |                        -0.0140572  |     -0.0242044  |      -0.00391005 |   0.96875 |         1 |
| mmd_rbf                   | M0_ot_interpolation | M9_full_pheromone |   5 |      0.0097388  |         0.0309501 |                        -0.0212113  |     -0.0241942  |      -0.0183683  |   1       |         1 |
| energy                    | M0_ot_interpolation | M9_full_pheromone |   5 |      0.345831   |         1.05423   |                        -0.7084     |     -0.773996   |      -0.651966   |   1       |         1 |
| celltype_composition_rmse | M0_ot_interpolation | M9_full_pheromone |   5 |      0.00381493 |         0.0096355 |                        -0.00582057 |     -0.00733134 |      -0.00431828 |   1       |         1 |

## Current Interpretation

Modules are retained only when they improve held-out reconstruction or provide a mechanistic diagnostic. If the full model does not beat the strongest baseline on at least two primary metrics, the manuscript must state that the Nature-level claim is not supported.
