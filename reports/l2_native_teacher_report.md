# L2 Native Teacher Report

- teacher_backend: fallback_sinkhorn
- native_moscot_available: False
- native_moscot_detail: import timed out after 90s
- native_moscot_used: False
- fallback_used: True
- native_failure_reason: native import unavailable: import timed out after 90s
- runtime_seconds: 0.10

| method_label   |   source_time |   target_time | file                                                      |   n_source |   n_target | solver         | teacher_backend       |   transport_cost |   mean_entropy |   mean_growth |   runtime_seconds | native_moscot_available   | native_moscot_used   | fallback_used   |   native_max_cells_per_time |   epsilon |
|:---------------|--------------:|--------------:|:----------------------------------------------------------|-----------:|-----------:|:---------------|:----------------------|-----------------:|---------------:|--------------:|------------------:|:--------------------------|:---------------------|:----------------|----------------------------:|----------:|
| teacher        |             6 |             9 | processed\external_l2\ot_couplings\teacher_t6_to_t9.npz   |         87 |        120 | numpy_sinkhorn | toy_sinkhorn_fallback |         0.615251 |       0.607323 |      6.20669  |         0.0993623 | False                     | False                | True            |                         120 |      0.08 |
| teacher        |             9 |            12 | processed\external_l2\ot_couplings\teacher_t9_to_t12.npz  |        120 |        120 | numpy_sinkhorn | toy_sinkhorn_fallback |         0.5874   |       0.596861 |      1.13072  |         0.0993623 | False                     | False                | True            |                         120 |      0.08 |
| teacher        |            12 |            15 | processed\external_l2\ot_couplings\teacher_t12_to_t15.npz |        120 |        120 | numpy_sinkhorn | toy_sinkhorn_fallback |         0.540161 |       0.576608 |      0.996557 |         0.0993623 | False                     | False                | True            |                         120 |      0.08 |
| teacher        |            15 |            21 | processed\external_l2\ot_couplings\teacher_t15_to_t21.npz |        120 |        120 | numpy_sinkhorn | toy_sinkhorn_fallback |         0.86324  |       0.821129 |      1.67503  |         0.0993623 | False                     | False                | True            |                         120 |      0.08 |
| teacher        |            21 |            28 | processed\external_l2\ot_couplings\teacher_t21_to_t28.npz |        120 |        120 | numpy_sinkhorn | toy_sinkhorn_fallback |         0.645112 |       0.666878 |      1.59734  |         0.0993623 | False                     | False                | True            |                         120 |      0.08 |
