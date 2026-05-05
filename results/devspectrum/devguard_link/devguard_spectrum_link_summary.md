# DevGuard-DevSpectrum Link Summary

## DTI correlations

| cohort       | devguard_metric   | spectral_metric          |   spearman_rho |   p_value |   n_lineages |
|:-------------|:------------------|:-------------------------|---------------:|----------:|-------------:|
| Tal1_chimera | DTI               | global_spectral_distance |     -0.103755  | 0.637557  |           23 |
| Tal1_chimera | DTI               | mean_absolute_residual   |      0.282609  | 0.191353  |           23 |
| Tal1_chimera | DTI               | max_absolute_residual    |     -0.0494071 | 0.822858  |           23 |
| T_chimera    | DTI               | global_spectral_distance |     -0.398128  | 0.0324358 |           29 |
| T_chimera    | DTI               | mean_absolute_residual   |     -0.369303  | 0.0486549 |           29 |
| T_chimera    | DTI               | max_absolute_residual    |     -0.299088  | 0.114999  |           29 |

## Top failure-mode spectral signatures

| cohort       | normality_class     | module_name             |   mean_residual |   mean_abs_residual |   median_abs_residual |   max_abs_residual |   n_records |   n_cells |
|:-------------|:--------------------|:------------------------|----------------:|--------------------:|----------------------:|-------------------:|------------:|----------:|
| Tal1_chimera | within_stage_normal | erythroid               |      -0.798747  |           0.798747  |             0.775184  |           1.77418  |          21 |     16286 |
| Tal1_chimera | abnormal_off_normal | erythroid               |      -0.787538  |           0.787538  |             0.795395  |           1.13726  |          18 |      4553 |
| Tal1_chimera | fate_deviation      | erythroid               |      -0.774738  |           0.774738  |             0.781186  |           1.13308  |          19 |      7385 |
| T_chimera    | abnormal_off_normal | erythroid               |      -0.0843081 |           0.32621   |             0.270778  |           0.885862 |           8 |       168 |
| T_chimera    | within_stage_normal | erythroid               |      -0.226463  |           0.227363  |             0.195166  |           1.18062  |          27 |     14920 |
| T_chimera    | fate_deviation      | erythroid               |      -0.18713   |           0.212061  |             0.180735  |           0.717586 |          18 |       627 |
| T_chimera    | abnormal_off_normal | extraembryonic_mesoderm |       0.0319163 |           0.0989667 |             0.0724152 |           0.313394 |           8 |       168 |
| T_chimera    | abnormal_off_normal | cell_cycle              |      -0.046766  |           0.0961796 |             0.0988274 |           0.185513 |           8 |       168 |
| T_chimera    | fate_deviation      | cell_cycle              |      -0.0736219 |           0.0897555 |             0.073136  |           0.229079 |          18 |       627 |
| T_chimera    | abnormal_off_normal | cardiac                 |       0.0357476 |           0.085696  |             0.0395676 |           0.401161 |           8 |       168 |
