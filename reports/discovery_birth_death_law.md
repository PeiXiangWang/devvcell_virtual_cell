# Discovery Birth Death Law

Birth/death hardening uses empirical stochastic event logs, predicted hazards, density permutation controls and seed-wise event stability.

## Tier

- tier: fail
- density effect: 9.33802e-05 [-0.000395111, 0.000551136]
- permutation_q: 0.019802
- seed_stability_pass: False (sign consistency=0.600)
- empirical_hazard_correlation: nan
- calibration_slope: nan
- calibration_rmse: nan
- carrying_capacity_like_threshold: False
- rollout_based_event_log: True

## Hazard Regression

| response          | predictor        |         coef |    abs_coef |        r2 |     n |
|:------------------|:-----------------|-------------:|------------:|----------:|------:|
| net_growth_hazard | local_density    |  9.33802e-05 | 9.33802e-05 | 0.0325562 | 40000 |
| net_growth_hazard | ot_growth        |  0.00081594  | 0.00081594  | 0.0325562 | 40000 |
| net_growth_hazard | cell_cycle_score | -2.1684e-19  | 2.1684e-19  | 0.0325562 | 40000 |
| net_growth_hazard | fate_entropy     |  0.000795409 | 0.000795409 | 0.0325562 | 40000 |
| net_growth_hazard | cci_signal       |  1.8602e-06  | 1.8602e-06  | 0.0325562 | 40000 |

## Event Counts By Lineage

| lineage         |   birth |   death |
|:----------------|--------:|--------:|
| erythroid       |    4118 |    4506 |
| mesoderm_muscle |    4640 |    5127 |
| neural          |    4509 |    5284 |

## Event Counts By Time

|   time |   birth |   death |
|-------:|--------:|--------:|
|  14.25 |    3206 |    3871 |
|  14.5  |    3621 |    3591 |
|  14.75 |    3012 |    3790 |
|  15    |    3428 |    3665 |
