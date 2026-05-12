# Scientific Gap Audit

Scientific goal after re-scoping: OT gives the developmental map; SwarmLineage-OT learns microscopic rules that realize the map and reveal emergent developmental laws.

- best mean-rank reconstruction row: `M0b_ot_interpolation`
- OT reference row: `M0b_ot_interpolation`
- teacher_fidelity_gate: True
- emergent_law_gate: True
- mechanistic_usefulness_gate: True

## Baseline/Reference Execution Matrix

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

## Remaining Gaps

- Native moscot TemporalProblem or externally validated teacher is still required for strong claims beyond a toy fallback.
- Emergent laws must be checked across seeds, held-out times and at least one external developmental dataset.
- CCI and memory laws are computational hypotheses unless supported by matched spatial, perturbation or lineage-tracing data.
- The manuscript must not state that SwarmLineage-OT outperforms OT; the correct claim is teacher fidelity plus mechanistic discovery.
