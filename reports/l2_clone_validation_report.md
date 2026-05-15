# L2 Clone-Aware Validation Report

Clone branch splitting score is terminal lineage entropy scaled by terminal clone size. Condensation exposure is the negative mean pre-event latent distance from the time-specific centroid, so larger values indicate stronger pre-event condensation.

| dataset_id           | l2_validation_tier   | teacher_backend   | clone_metadata_loaded   |   usable_clone_count |   condensation_spearman |   condensation_spearman_p |   condensation_pearson |   condensation_pearson_p |   condensation_permutation_p | negative_controls_pass   | sensitivity_stable   | condensation_predicts_clone_branch_splitting   | interpretation                                                                                      |
|:---------------------|:---------------------|:------------------|:------------------------|---------------------:|------------------------:|--------------------------:|-----------------------:|-------------------------:|-----------------------------:|:-------------------------|:---------------------|:-----------------------------------------------|:----------------------------------------------------------------------------------------------------|
| L2_biddy_2018_nature | fail                 | fallback_sinkhorn | True                    |                  214 |               -0.352496 |                1.1747e-07 |              -0.351001 |              1.33941e-07 |                  0.000999001 | True                     | True                 | False                                          | Clone-aware validation does not support the hypothesis in Biddy_2018 under this operationalization. |

Regression:

| term                                  |   coefficient |   std_error |   t_stat |     p_value |
|:--------------------------------------|--------------:|------------:|---------:|------------:|
| intercept                             |  -1.373       | 0.195616    | -7.01886 | 3.07123e-11 |
| clone_pre_event_condensation_exposure |  -0.0529902   | 0.0173596   | -3.0525  | 0.00256401  |
| clone_size                            |   0.000246553 | 7.77967e-05 |  3.1692  | 0.00175789  |
| clone_time_span                       |   0.0335028   | 0.00479162  |  6.99195 | 3.58844e-11 |
| clone_start_time                      |   0.0525975   | 0.0066043   |  7.96413 | 1.05471e-13 |

Controls:

| control                 |   mean_null_spearman |      ci_low |      ci_high | control_signal_abs_ge_observed   |
|:------------------------|---------------------:|------------:|-------------:|:---------------------------------|
| clone_id_shuffle        |         -0.00767232  | -0.0156725  |  0.000506378 | False                            |
| time_shuffle            |         -0.230554    | -0.237246   | -0.223681    | False                            |
| branch_label_shuffle    |         -0.000731523 | -0.0108095  |  0.00839916  | False                            |
| order_parameter_shuffle |          0.00172947  | -0.00809716 |  0.0105411   | False                            |
| random_teacher_velocity |         -0.00217895  | -0.0120193  |  0.00745346  | False                            |

Sensitivity:

| sensitivity        |   usable_clones |   spearman |     p_value |
|:-------------------|----------------:|-----------:|------------:|
| min_clone_size_5   |             398 |  -0.196441 | 7.9745e-05  |
| min_clone_size_10  |             214 |  -0.352496 | 1.1747e-07  |
| min_clone_size_20  |             122 |  -0.365042 | 3.56334e-05 |
| min_clone_size_50  |              39 |  -0.214348 | 0.190081    |
| downsample_seed_7  |             179 |  -0.334967 | 4.57773e-06 |
| downsample_seed_17 |             176 |  -0.318338 | 1.66465e-05 |
| downsample_seed_23 |             175 |  -0.337601 | 4.90307e-06 |
| downsample_seed_31 |             183 |  -0.317774 | 1.16804e-05 |
| downsample_seed_43 |             182 |  -0.28341  | 0.000105816 |
