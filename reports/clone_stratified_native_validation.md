# Clone-Stratified Native Validation

This analysis tests whether the weak Jindal full-data fallback association was a teacher artifact, downsample artifact or clone-coverage artifact. Jindal and Weinreb were rerun under five native sampling strategies designed to increase clone coverage while retaining native moscot temporal couplings.

## Final Interpretation

- status: weinreb_sampling_specific_condensation_signal
- interpretation: Jindal fallback positivity does not survive clone-stratified native sampling, while Weinreb shows a sampling-specific native condensation signal; this is mixed weak evidence rather than general clone support.

## Strategy Summary

| dataset_id                              | sampling_strategy          |   n_cells |   usable_clones |   clone_coverage_ratio | teacher_backend   | native_moscot_success   | best_model           | best_model_tier   |   best_model_effect | condensation_only_tier   | uncertainty_plus_teacher_bias_tier   |
|:----------------------------------------|:---------------------------|----------:|----------------:|-----------------------:|:------------------|:------------------------|:---------------------|:------------------|--------------------:|:-------------------------|:-------------------------------------|
| Jindal_2023_NatureBiotechnology_LSK_RNA | cell_random_native         |      1000 |              42 |              0.0303468 | native_moscot     | True                    | condensation_only    | fail              |           0.107362  | fail                     | fail                                 |
| Jindal_2023_NatureBiotechnology_LSK_RNA | clone_stratified_native    |      1000 |             500 |              0.361272  | native_moscot     | True                    | condensation_only    | fail              |         nan         | fail                     | fail                                 |
| Jindal_2023_NatureBiotechnology_LSK_RNA | clone_time_balanced_native |      1000 |             250 |              0.180636  | native_moscot     | True                    | condensation_only    | fail              |           0.0211439 | fail                     | fail                                 |
| Jindal_2023_NatureBiotechnology_LSK_RNA | clone_fate_balanced_native |      1000 |              59 |              0.0426301 | native_moscot     | True                    | post_divergence_only | acceptable        |           0.236541  | fail                     | fail                                 |
| Jindal_2023_NatureBiotechnology_LSK_RNA | max_feasible_native        |      1300 |             325 |              0.234827  | native_moscot     | True                    | post_divergence_only | fail              |           0.0237771 | fail                     | fail                                 |
| Weinreb_2020_Science                    | cell_random_native         |      1500 |              95 |              0.0244782 | native_moscot     | True                    | post_divergence_only | acceptable        |           0.379708  | fail                     | fail                                 |
| Weinreb_2020_Science                    | clone_stratified_native    |      1500 |             500 |              0.128833  | native_moscot     | True                    | condensation_only    | fail              |         nan         | fail                     | fail                                 |
| Weinreb_2020_Science                    | clone_time_balanced_native |      1500 |             258 |              0.0664777 | native_moscot     | True                    | condensation_only    | acceptable        |           0.132927  | acceptable               | fail                                 |
| Weinreb_2020_Science                    | clone_fate_balanced_native |      1500 |             181 |              0.0466375 | native_moscot     | True                    | fate_entropy_only    | fail              |           0.108658  | fail                     | fail                                 |
| Weinreb_2020_Science                    | max_feasible_native        |      1950 |             337 |              0.0868333 | native_moscot     | True                    | condensation_only    | acceptable        |           0.172378  | acceptable               | fail                                 |

## Native Runs

| dataset_id                              | sampling_strategy          | native_moscot_success   | teacher_backend   |   runtime_seconds |   n_pairs | plan_shapes           | failure_reason   |
|:----------------------------------------|:---------------------------|:------------------------|:------------------|------------------:|----------:|:----------------------|:-----------------|
| Jindal_2023_NatureBiotechnology_LSK_RNA | cell_random_native         | True                    | native_moscot     |           9.8042  |         1 | (500, 500)            |                  |
| Jindal_2023_NatureBiotechnology_LSK_RNA | clone_stratified_native    | True                    | native_moscot     |          10.1674  |         1 | (500, 500)            |                  |
| Jindal_2023_NatureBiotechnology_LSK_RNA | clone_time_balanced_native | True                    | native_moscot     |          11.262   |         1 | (500, 500)            |                  |
| Jindal_2023_NatureBiotechnology_LSK_RNA | clone_fate_balanced_native | True                    | native_moscot     |          10.0139  |         1 | (500, 500)            |                  |
| Jindal_2023_NatureBiotechnology_LSK_RNA | max_feasible_native        | True                    | native_moscot     |           9.69835 |         1 | (650, 650)            |                  |
| Weinreb_2020_Science                    | cell_random_native         | True                    | native_moscot     |           9.64921 |         2 | (500, 500);(500, 500) |                  |
| Weinreb_2020_Science                    | clone_stratified_native    | True                    | native_moscot     |          12.8431  |         2 | (500, 500);(500, 500) |                  |
| Weinreb_2020_Science                    | clone_time_balanced_native | True                    | native_moscot     |          13.1148  |         2 | (500, 500);(500, 500) |                  |
| Weinreb_2020_Science                    | clone_fate_balanced_native | True                    | native_moscot     |          11.85    |         2 | (500, 500);(500, 500) |                  |
| Weinreb_2020_Science                    | max_feasible_native        | True                    | native_moscot     |          12.6366  |         2 | (650, 650);(650, 650) |                  |

## Pre-Registered Models

Models tested: condensation_only, alignment_only, fate_entropy_only, teacher_bias_only, post_divergence_only, condensation_plus_entropy, condensation_plus_post_divergence, uncertainty_plus_teacher_bias and full_branch_window_model. Outcomes tested: terminal_fate_entropy as the primary outcome plus terminal_lineage_entropy, clone_multilineage_score, clone_branch_count, clone_fate_diversification_index and clone_transition_entropy as sensitivity outcomes.

Clone-stratified strategies increased clone coverage, but one-cell-per-terminal-clone sampling can make terminal fate entropy uninformative; time-balanced and max-feasible strategies partially address that tradeoff. No result is interpreted as experimental, causal or wet-lab validation. Clone-aware support requires native/validated teacher support plus covariate-adjusted and matched positive association across datasets, not a single strategy-specific signal.

## Supported or Weak Primary-Outcome Rows

| dataset_id                              | sampling_strategy          | model                     | outcome               |   n_clones |   raw_effect |   bootstrap_ci_low |   bootstrap_ci_high |   covariate_adjusted_effect |   within_start_time_bin_effect |   clone_size_matched_effect |   time_span_matched_effect |   permutation_q | support_tier   |
|:----------------------------------------|:---------------------------|:--------------------------|:----------------------|-----------:|-------------:|-------------------:|--------------------:|----------------------------:|-------------------------------:|----------------------------:|---------------------------:|----------------:|:---------------|
| Jindal_2023_NatureBiotechnology_LSK_RNA | clone_fate_balanced_native | post_divergence_only      | terminal_fate_entropy |         59 |    0.236541  |        -0.0315075  |            0.490271 |                   0.200793  |                      0.236541  |                   0.207369  |                  0.266975  |      0.049505   | acceptable     |
| Weinreb_2020_Science                    | cell_random_native         | post_divergence_only      | terminal_fate_entropy |         95 |    0.379708  |         0.228851   |            0.530261 |                   0.157744  |                      0.264917  |                   0.432709  |                  0.313484  |      0.00990099 | acceptable     |
| Weinreb_2020_Science                    | clone_time_balanced_native | condensation_only         | terminal_fate_entropy |        258 |    0.132927  |         0.0287587  |            0.272554 |                   0.122107  |                      0.132927  |                   0.0995088 |                  0.123518  |      0.019802   | acceptable     |
| Weinreb_2020_Science                    | clone_time_balanced_native | condensation_plus_entropy | terminal_fate_entropy |        258 |    0.115309  |        -0.00621189 |            0.268728 |                   0.0982432 |                      0.115309  |                   0.0917962 |                  0.103601  |      0.0792079  | weak           |
| Weinreb_2020_Science                    | max_feasible_native        | condensation_only         | terminal_fate_entropy |        337 |    0.172378  |         0.0654015  |            0.259299 |                   0.160313  |                      0.172378  |                   0.162821  |                  0.171507  |      0.00990099 | acceptable     |
| Weinreb_2020_Science                    | max_feasible_native        | fate_entropy_only         | terminal_fate_entropy |        337 |    0.0893027 |        -0.0325401  |            0.170883 |                   0.0523473 |                      0.0893027 |                   0.0839121 |                  0.0874455 |      0.0792079  | weak           |
| Weinreb_2020_Science                    | max_feasible_native        | condensation_plus_entropy | terminal_fate_entropy |        337 |    0.140343  |         0.0309879  |            0.224415 |                   0.109432  |                      0.140343  |                   0.130832  |                  0.134652  |      0.019802   | acceptable     |
