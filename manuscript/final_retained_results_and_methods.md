# Final Retained Results and Methods

## Central Claim

OT gives the developmental map; SwarmLineage-OT learns microscopic finite-agent rules that realize the map and reveal emergent developmental laws.

The retained manuscript must not claim that SwarmLineage-OT outperforms OT interpolation. `M0b_ot_interpolation` is an oracle-like teacher/reference interpolation.

## Evaluation Gates

- teacher_fidelity_gate: True
- emergent_law_gate: True
- mechanistic_usefulness_gate: True

## Retained Metrics

- teacher fidelity metrics: `tables/teacher_fidelity_metrics.csv`
- emergent law gates: `tables/emergent_law_gate_summary.csv`

| model                                 |   sinkhorn |    mmd_rbf |   energy |   celltype_composition_rmse |
|:--------------------------------------|-----------:|-----------:|---------:|----------------------------:|
| M0_linear_label_interpolation         |   0.281661 | 0.00822039 | 0.292374 |                  0.00396321 |
| M0b_ot_interpolation                  |   0.255076 | 0.00878499 | 0.335883 |                  0.00348103 |
| M10_shuffled_time_ot                  |   0.31024  | 0.00516795 | 0.204379 |                  0.018496   |
| M11_random_lr_labels                  |   0.289617 | 0.0172297  | 0.628864 |                  0.0233767  |
| M1_intrinsic_neural                   |   0.309506 | 0.0276649  | 0.915083 |                  0.00387904 |
| M2_ot_teacher_force                   |   0.274306 | 0.0167648  | 0.568815 |                  0.00386852 |
| M3_ot_birth_death                     |   0.288529 | 0.0164293  | 0.590722 |                  0.0168219  |
| M4_ot_adaptive_diffusion              |   0.274514 | 0.0167285  | 0.566489 |                  0.00420092 |
| M5_ot_swarm                           |   0.273951 | 0.0178781  | 0.613508 |                  0.00448128 |
| M6_ot_swarm_birth_death               |   0.292737 | 0.017475   | 0.630969 |                  0.0209927  |
| M7_ot_swarm_birth_death_diffusion     |   0.284228 | 0.0183164  | 0.632902 |                  0.0200099  |
| M8_ot_swarm_birth_death_diffusion_cci |   0.284235 | 0.0183179  | 0.632942 |                  0.0200099  |
| M9_full_memory                        |   0.28393  | 0.0183946  | 0.637903 |                  0.0200099  |

## Methods Retained

- strict time holdout and teacher-edge holdout support leakage-resistant evaluation.
- OT teacher construction records backend status; toy fallback is not presented as native moscot/WOT.
- `SwarmLineageDynamics` represents trainable intrinsic, teacher, swarm, birth/death, adaptive diffusion, CCI and memory components.
- Stochastic birth/death uses event simulation and writes event logs.
- Discovery modules estimate diffusion, growth, branch nucleation, memory hysteresis, CCI branch bias and phase-regime laws.

## Interpretation

If `M0b_ot_interpolation` has the lowest reconstruction error, this is expected for a teacher/reference. The agent model is retained when it stays close enough to the teacher and yields stable mechanistic laws.

## Limitations

- Current emergent laws are computational hypotheses, not validated biological mechanisms.
- Strong biological claims require native moscot/WOT or external teacher validation plus external or perturbation validation.
- If teacher fidelity is poor and no emergent-law gates are stable, the current scientific hypothesis should be reported as unsupported.
