# External Validation Summary

E1 remains the only acceptable external time-series support. L1 clone-aware testing was attempted but did not support the clone-level hypothesis. E2 is weak/proxy feasibility only.

## E1

| dataset_id                               | external_validation_tier   | external_teacher_backend   | branch_event_detected   | condensation_before_divergence_reproduced   | lineage_validated   |   lineage_separation_effect |   alignment_effect |   entropy_effect | negative_control_pass   | seed_stability_pass   | interpretation                           | status   |   blocker |
|:-----------------------------------------|:---------------------------|:---------------------------|:------------------------|:--------------------------------------------|:--------------------|----------------------------:|-------------------:|-----------------:|:------------------------|:----------------------|:-----------------------------------------|:---------|----------:|
| E1_mouse_gastrulation_wt_chimera_sample1 | acceptable                 | native_moscot              | True                    | True                                        | False               |                    -1.12489 |          0.0195907 |                0 | True                    | True                  | transient_condensation_before_divergence | analyzed |       nan |

## L1

| dataset_id              | lineage_validation_tier   | teacher_backend     |   clone_count_tested |   condensation_to_clone_splitting_spearman |   condensation_permutation_p |   density_to_clone_splitting_spearman |   density_spearman_p | negative_controls_pass   | interpretation                                                                                                                | status                  |
|:------------------------|:--------------------------|:--------------------|---------------------:|-------------------------------------------:|-----------------------------:|--------------------------------------:|---------------------:|:-------------------------|:------------------------------------------------------------------------------------------------------------------------------|:------------------------|
| L1_kim_2020_cellreports | fail                      | not_run_clone_proxy |                   21 |                                  -0.104584 |                     0.668663 |                             -0.517459 |            0.0162863 | True                     | Clone metadata were analyzed, but condensation exposure did not predict clone branch splitting under this operationalization. | analyzed_not_supportive |

## E2

| dataset_id              | e2_validation_tier   | teacher_backend        |   cells_analyzed |   time_points | branch_event_detected   | condensation_direction_observed   | lineage_validated   | interpretation                                                                                    |
|:------------------------|:---------------------|:-----------------------|-----------------:|--------------:|:------------------------|:----------------------------------|:--------------------|:--------------------------------------------------------------------------------------------------|
| E2_GSE212050_gastruloid | weak                 | not_run_temporal_proxy |             1800 |             5 | True                    | False                             | False               | Local time-series feasibility support only; no native teacher or clone-level test was run for E2. |
