# Topological Rule Specificity v1.3

- conclusion: topological_and_metric_both_supported_but_not_specific
- Interpretation: topological kNN can reproduce the branch signature, but metric-radius neighbours can also reproduce it and random controls are not fully clean. This blocks a strong topological-specific claim.

| neighbor_rule   |   internal_condensation_rate |   E1_condensation_rate |   L2_condensation_rate |   downsample_condensation_rate |   mean_internal_effect |   mean_E1_effect | negative_control_clean   | conclusion                                             |
|:----------------|-----------------------------:|-----------------------:|-----------------------:|-------------------------------:|-----------------------:|-----------------:|:-------------------------|:-------------------------------------------------------|
| topological     |                          1   |                      1 |                    0   |                       1        |             -0.70972   |        -0.376148 | True                     | topological_and_metric_both_supported_but_not_specific |
| metric          |                          1   |                      1 |                    0   |                       0.966667 |             -0.424855  |        -0.286628 | True                     | topological_and_metric_both_supported_but_not_specific |
| random          |                          0.4 |                      1 |                    0   |                       0.733333 |              0.0105956 |        -0.10529  | False                    | topological_and_metric_both_supported_but_not_specific |
| label           |                          0   |                      1 |                    0   |                     nan        |              0         |        -0.683283 | True                     | topological_and_metric_both_supported_but_not_specific |
| mixed           |                          1   |                      1 |                    0.7 |                     nan        |             -0.750032  |        -0.293695 | True                     | topological_and_metric_both_supported_but_not_specific |
