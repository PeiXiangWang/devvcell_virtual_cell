# Discovery Birth Death Law

Birth/death hardening uses empirical stochastic event logs, predicted hazards, density permutation controls and seed-wise event stability.

## Tier

- tier: fail
- density effect: 7.89754e-05 [-0.000477091, 0.000584699]
- permutation_q: 0.019802
- seed_stability_pass: False (sign consistency=0.600)
- empirical_hazard_correlation: nan
- calibration_slope: nan
- calibration_rmse: nan
- carrying_capacity_like_threshold: False
- rollout_based_event_log: True

## Hazard Regression

| response          | predictor        |        coef |    abs_coef |        r2 |     n |
|:------------------|:-----------------|------------:|------------:|----------:|------:|
| net_growth_hazard | local_density    | 7.89754e-05 | 7.89754e-05 | 0.0259699 | 40000 |
| net_growth_hazard | ot_growth        | 0.000833337 | 0.000833337 | 0.0259699 | 40000 |
| net_growth_hazard | cell_cycle_score | 9.75782e-19 | 9.75782e-19 | 0.0259699 | 40000 |
| net_growth_hazard | fate_entropy     | 0.000577423 | 0.000577423 | 0.0259699 | 40000 |
| net_growth_hazard | cci_signal       | 1.80344e-05 | 1.80344e-05 | 0.0259699 | 40000 |

## Event Counts By Lineage

| lineage         |   birth |   death |
|:----------------|--------:|--------:|
| erythroid       |    4228 |    4485 |
| mesoderm_muscle |    4795 |    5116 |
| neural          |    4282 |    5326 |

## Event Counts By Time

|   time |   birth |   death |
|-------:|--------:|--------:|
|  14.25 |    3192 |    3871 |
|  14.5  |    3616 |    3568 |
|  14.75 |    3056 |    3813 |
|  15    |    3441 |    3675 |
