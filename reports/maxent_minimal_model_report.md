# Maximum-Entropy Minimal Model

- maxent_model_tier: acceptable
- selected_k: 2

The model uses only pairwise topological label similarity and velocity alignment. It is not trained on branch-event labels. Because it is a prototype, a weak or failed result is interpreted as limited support for the minimal-model explanation rather than as a biological negative result.

## Fit

| dataset   |   k |   J_pairwise_label |   h_velocity_alignment |   pair_count |   same_label_rate |   random_same_label_rate |
|:----------|----:|-------------------:|-----------------------:|-------------:|------------------:|-------------------------:|
| internal  |   2 |            5.04614 |              0.0665986 |        10000 |          0.9874   |                0.3352    |
| E1        |   2 |            4.12307 |              0.464853  |         3600 |          0.803056 |                0.0619444 |
| L2        |   2 |            1.2472  |              0.483257  |         6400 |          0.9225   |                0.77375   |
| E2        |   2 |            1.2696  |              0.892242  |         3600 |          0.739444 |                0.443611  |

## Prediction

| dataset   |   k |   observed_lineage_separation_effect |   predicted_lineage_separation_effect |   observed_alignment_effect |   predicted_alignment_effect | condensation_direction_match   | alignment_direction_match   |
|:----------|----:|-------------------------------------:|--------------------------------------:|----------------------------:|-----------------------------:|:-------------------------------|:----------------------------|
| internal  |   2 |                             -0.38623 |                             -0.834605 |                   0.0062109 |                    0.0434396 | True                           | True                        |
| E1        |   2 |                             -1.11797 |                             -0.804805 |                   0.0882369 |                    0.0452177 | True                           | True                        |
| L2        |   2 |                              1.62105 |                             -0.555002 |                  -0.325421  |                    0.0330632 | False                          | False                       |
| E2        |   2 |                              1.50861 |                             -0.559394 |                   0         |                    0.0418286 | False                          | False                       |
