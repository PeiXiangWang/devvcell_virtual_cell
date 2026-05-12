# Scientific Gap Audit

- Best mean-rank model: `M0b_ot_interpolation`.
- Strongest baseline: `M0b_ot_interpolation`.
- Full model passes the predefined gate of beating the strongest baseline on at least two core metrics: False.
- High-level claims are allowed only with native_moscot/native_wot or externally validated teacher. The current fallback teacher remains toy_sinkhorn_fallback unless reports say otherwise.

## Baseline Execution Matrix

| baseline                              | available   | executed   | status                                                                                                                                                                                                                             | module        |
|:--------------------------------------|:------------|:-----------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------------|
| M0_linear_label_interpolation         | True        | True       | evaluated_in_this_run                                                                                                                                                                                                              | nan           |
| M0b_ot_interpolation                  | True        | True       | evaluated_in_this_run                                                                                                                                                                                                              | nan           |
| M10_shuffled_time_ot                  | True        | True       | evaluated_in_this_run                                                                                                                                                                                                              | nan           |
| M11_random_lr_labels                  | True        | True       | evaluated_in_this_run                                                                                                                                                                                                              | nan           |
| M1_intrinsic_neural                   | True        | True       | evaluated_in_this_run                                                                                                                                                                                                              | nan           |
| M2_ot_teacher_force                   | True        | True       | evaluated_in_this_run                                                                                                                                                                                                              | nan           |
| M3_ot_birth_death                     | True        | True       | evaluated_in_this_run                                                                                                                                                                                                              | nan           |
| M4_ot_adaptive_diffusion              | True        | True       | evaluated_in_this_run                                                                                                                                                                                                              | nan           |
| M5_ot_swarm                           | True        | True       | evaluated_in_this_run                                                                                                                                                                                                              | nan           |
| M6_ot_swarm_birth_death               | True        | True       | evaluated_in_this_run                                                                                                                                                                                                              | nan           |
| M7_ot_swarm_birth_death_diffusion     | True        | True       | evaluated_in_this_run                                                                                                                                                                                                              | nan           |
| M8_ot_swarm_birth_death_diffusion_cci | True        | True       | evaluated_in_this_run                                                                                                                                                                                                              | nan           |
| M9_full_memory                        | True        | True       | evaluated_in_this_run                                                                                                                                                                                                              | nan           |
| CellRank2                             | False       | False      | import_timeout>45s                                                                                                                                                                                                                 | cellrank      |
| TrajectoryNet                         | False       | False      | CalledProcessError:Command '['C:\\Users\\14915\\AppData\\Local\\Programs\\Python\\Python311\\python.exe', '-c', "import TrajectoryNet; print(getattr(TrajectoryNet, '__version__', 'unknown'))"]' returned non-zero exit status 1. | TrajectoryNet |
| MIOFlow                               | False       | False      | CalledProcessError:Command '['C:\\Users\\14915\\AppData\\Local\\Programs\\Python\\Python311\\python.exe', '-c', "import mioflow; print(getattr(mioflow, '__version__', 'unknown'))"]' returned non-zero exit status 1.             | mioflow       |
| TIGON                                 | False       | False      | CalledProcessError:Command '['C:\\Users\\14915\\AppData\\Local\\Programs\\Python\\Python311\\python.exe', '-c', "import tigon; print(getattr(tigon, '__version__', 'unknown'))"]' returned non-zero exit status 1.                 | tigon         |

## Required Next Steps Before Strong Claims

- Run native moscot TemporalProblem and extract native transport plans, or validate the fallback teacher externally.
- Add real external held-out dataset validation and lineage/perturbation validation.
- Tune trainable swarm, event and memory coefficients only on training times, then re-run the strict holdout gate.
- Keep negative controls in the report; shuffled time/LR controls must degrade performance.
