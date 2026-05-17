# Integrated Clone-Aware Validation Conclusion

The clone-aware developmental expansion does not establish clone-level fate-diversification support. The strongest retained project claim remains the time-series branch-nucleation order-parameter hypothesis supported by internal native moscot and E1 MouseGastrulationData.

## Dataset-Level Result

- Jindal LSK: full-data fallback showed a weak positive primary association, but downsampled native moscot did not retain it.
- Weinreb LARRY: primary condensation was not supported in full-data fallback and was negative in downsampled native moscot.
- Xie organoid: clone/barcode metadata are present, but the processed h5ad lacks an explicit time/stage field for branch-window validation.

## Integrated Summary

| dataset_id                              | full_data_teacher_backend   |   full_data_usable_clones | full_data_primary_condensation_support   |   full_data_primary_effect |   full_data_primary_q |   full_data_covariate_adjusted_effect | native_teacher_backend    |   native_usable_clones | native_primary_condensation_support   |   native_primary_effect |   native_primary_q |   native_covariate_adjusted_effect |   native_revised_support_count | final_clone_aware_tier   | final_interpretation    |
|:----------------------------------------|:----------------------------|--------------------------:|:-----------------------------------------|---------------------------:|----------------------:|--------------------------------------:|:--------------------------|-----------------------:|:--------------------------------------|------------------------:|-------------------:|-----------------------------------:|-------------------------------:|:-------------------------|:------------------------|
| Jindal_2023_NatureBiotechnology_LSK_RNA | fallback_centroid_teacher   |                       379 | True                                     |                  0.340983  |           1.00683e-10 |                              0.417369 | native_moscot_downsampled |                     42 | False                                 |               0.0524734 |          0.961786  |                           0.181401 |                              0 | fail                     | native_no_clone_support |
| Weinreb_2020_Science                    | fallback_centroid_teacher   |                      1373 | False                                    |                  0.0265663 |           0.419426    |                              0.144969 | native_moscot_downsampled |                     95 | False                                 |              -0.32322   |          0.0126828 |                          -0.331556 |                              0 | fail                     | native_no_clone_support |
| Xie_2023_NatureMethods_Organoid         | fallback_centroid_teacher   |                       nan | False                                    |                nan         |         nan           |                            nan        | not_available             |                    nan | False                                 |             nan         |        nan         |                         nan        |                            nan | fail                     | metadata_blocker        |

## Interpretation

The Jindal weak signal is best treated as fallback-teacher feasibility, not as retained clone-aware support, because it did not persist in the downsampled native-teacher check. Weinreb provides negative clone-aware evidence for primary condensation and does not support the revised two-phase or uncertainty-gated models in the native downsample. Xie remains a metadata blocker for branch-window testing.
