# Leakage Audit

- split_mode: strict_time_holdout
- native teacher excludes obs rows marked `eval_holdout` in `run_native_moscot_teacher`.
- holdout gap bridge edge is labelled below and must not be described as an ordinary adjacent observed edge.

|   source_time |   target_time | teacher_backend   | edge_type              |
|--------------:|--------------:|:------------------|:-----------------------|
|            13 |            14 | native_moscot     | adjacent_observed_edge |
|            12 |            13 | native_moscot     | adjacent_observed_edge |
|            14 |            16 | native_moscot     | holdout_gap_bridge     |
|            18 |            19 | native_moscot     | adjacent_observed_edge |
|            17 |            18 | native_moscot     | adjacent_observed_edge |
|            16 |            17 | native_moscot     | adjacent_observed_edge |
