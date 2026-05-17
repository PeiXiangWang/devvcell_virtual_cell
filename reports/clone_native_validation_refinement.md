# Clone Native Validation Refinement

Native moscot was run in `.venv_moscot_native` on downsampled Jindal and Weinreb AnnData objects. Associations here are restricted to clones represented in the native downsample and therefore test direction/stability rather than replacing the full-data fallback analysis.

| dataset_id                              | native_input_loaded   | teacher_backend           |   n_cells |   usable_clones | primary_condensation_support   |   primary_effect |   primary_q |   primary_covariate_adjusted_effect |   revised_support_count | validation_tier   | interpretation          |
|:----------------------------------------|:----------------------|:--------------------------|----------:|----------------:|:-------------------------------|-----------------:|------------:|------------------------------------:|------------------------:|:------------------|:------------------------|
| Jindal_2023_NatureBiotechnology_LSK_RNA | True                  | native_moscot_downsampled |      1000 |              42 | False                          |        0.0524734 |   0.961786  |                            0.181401 |                       0 | fail              | native_no_clone_support |
| Weinreb_2020_Science                    | True                  | native_moscot_downsampled |      1500 |              95 | False                          |       -0.32322   |   0.0126828 |                           -0.331556 |                       0 | fail              | native_no_clone_support |
