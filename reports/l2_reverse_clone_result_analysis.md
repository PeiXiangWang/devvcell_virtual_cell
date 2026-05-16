# L2 Reverse Clone Result Analysis

- interpretation: l2_primary_condensation_fails_with_alternative_exposure_confounding
- Biddy/CellTag remains a failed clone-aware test under the current operationalization. Multiple exposure definitions and covariate-adjusted regressions do not justify converting it into support.

## Exposure Sensitivity

| exposure_definition   |   n_clones |   spearman |     p_value | direction_supports_hypothesis   |
|:----------------------|-----------:|-----------:|------------:|:--------------------------------|
| condensation_exposure |        214 |  -0.352496 | 1.1747e-07  | False                           |
| alignment_exposure    |        214 |   0.14476  | 0.0343094   | True                            |
| entropy_exposure      |        214 |   0.651124 | 3.38582e-27 | True                            |
| density_exposure      |        214 |  -0.283499 | 2.55789e-05 | False                           |
| post_divergence       |        214 |   0.187033 | 0.00606411  | True                            |
| absolute_condensation |        214 |   0.352496 | 1.1747e-07  | True                            |

## Covariate-Adjusted Regression

| exposure_definition   |   n_clones |   standardized_exposure_coefficient |       r2 | direction_supports_hypothesis_after_covariates   |
|:----------------------|-----------:|------------------------------------:|---------:|:-------------------------------------------------|
| condensation_exposure |        214 |                          -0.18182   | 0.378602 | False                                            |
| alignment_exposure    |        214 |                           0.0707876 | 0.355704 | True                                             |
| entropy_exposure      |        214 |                           0.557559  | 0.63532  | True                                             |
| density_exposure      |        214 |                          -0.119392  | 0.363662 | False                                            |
| post_divergence       |        214 |                           0.125238  | 0.366547 | True                                             |
| absolute_condensation |        214 |                           0.18182   | 0.378602 | True                                             |
