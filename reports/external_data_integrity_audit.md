# External Data Integrity Audit

| dataset_id                               | ann_data_exists   | expression_matrix_shape   | time_stage_column                         | lineage_column               | clone_or_lineage_barcode   | invented_clone_labels   | internal_data_leakage_detected   | lineage_validated   | external_time_series_support   |
|:-----------------------------------------|:------------------|:--------------------------|:------------------------------------------|:-----------------------------|:---------------------------|:------------------------|:---------------------------------|:--------------------|:-------------------------------|
| E1_mouse_gastrulation_wt_chimera_sample1 | True              | 1800x2000                 | time_point/time_numeric from stage.mapped | lineage from celltype.mapped | absent                     | False                   | False                            | False               | True                           |

No clone/barcode field was used or invented for E1. E1 remains time-series support only.
