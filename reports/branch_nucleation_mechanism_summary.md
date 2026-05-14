# Branch Nucleation Mechanism Summary

- branch_nucleation_tier: strong
- best_interpretation: transient_condensation_before_divergence
- primary_model_hint: M5_ot_swarm
- unsupported modules to exclude from main claim: birth/death, CCI, memory.
- architectural controls can show related condensation signals, so the current evidence supports an order-parameter signature, not proof that swarm or teacher terms are necessary by themselves.
- Negative controls include temporal, velocity, lineage, fate, no-swarm, no-teacher and random-teacher controls.

## Model Comparison

| variant                               | branch_nucleation_tier   |   lineage_separation_effect |   effect_ci_low |   effect_ci_high |   permutation_p |   permutation_q | seed_stability_pass   |   sign_consistency |   n_seed_windows |   local_velocity_alignment_A_effect |   branch_cohesion_C_effect |   lineage_separation_S_effect |   fate_entropy_H_effect |   branch_imbalance_B_effect |   local_density_mean_effect |   n_agents_effect | best_interpretation                      |
|:--------------------------------------|:-------------------------|----------------------------:|----------------:|-----------------:|----------------:|----------------:|:----------------------|-------------------:|-----------------:|------------------------------------:|---------------------------:|------------------------------:|------------------------:|----------------------------:|----------------------------:|------------------:|:-----------------------------------------|
| M5_ot_swarm                           | strong                   |                   -0.265822 |       -0.294104 |        -0.242814 |      0.00990099 |      0.00990099 | True                  |                  1 |                5 |                           0.0369052 |                -0.00247282 |                     -0.265822 |               0         |                   0         |                 -0.00225804 |               0   | transient_condensation_before_divergence |
| M7_ot_swarm_birth_death_diffusion     | strong                   |                   -0.21572  |       -0.28427  |        -0.132548 |      0.00990099 |      0.00990099 | True                  |                  1 |                5 |                           0.0448192 |                -0.00766641 |                     -0.21572  |               0.0094626 |                   0.0252695 |                  0.010676   |             -27.8 | transient_condensation_before_divergence |
| M8_ot_swarm_birth_death_diffusion_cci | strong                   |                   -0.215772 |       -0.284331 |        -0.132587 |      0.00990099 |      0.00990099 | True                  |                  1 |                5 |                           0.0448182 |                -0.00766538 |                     -0.215772 |               0.0094626 |                   0.0252695 |                  0.0106693  |             -27.8 | transient_condensation_before_divergence |
| M9_full_memory                        | strong                   |                   -0.207681 |       -0.276311 |        -0.124319 |      0.00990099 |      0.00990099 | True                  |                  1 |                5 |                           0.0473033 |                -0.00771597 |                     -0.207681 |               0.0094626 |                   0.0252695 |                  0.0104832  |             -27.8 | transient_condensation_before_divergence |

## Negative Controls

| control                     |   effect_size |   effect_ci_low |   effect_ci_high | seed_stability_pass   | gate_tier   | reason                                                                                |   permutation_p |   permutation_q | variant              |
|:----------------------------|--------------:|----------------:|-----------------:|:----------------------|:------------|:--------------------------------------------------------------------------------------|----------------:|----------------:|:---------------------|
| shuffled_temporal_order     |     0.0251302 |      0.0172486  |        0.032537  | False                 | fail        | permutation/shuffle null control, expected not to reproduce retained branch signature |             0.8 |        0.933333 | nan                  |
| shuffled_velocity           |    -0.226249  |     -0.226249   |       -0.226249  | False                 | fail        | permutation/shuffle null control, expected not to reproduce retained branch signature |             0.6 |        0.84     | nan                  |
| shuffled_lineage_labels     |     0.0111182 |      0.00504945 |        0.0179362 | False                 | fail        | permutation/shuffle null control, expected not to reproduce retained branch signature |             1   |        1        | nan                  |
| shuffled_fate_probabilities |    -0.230224  |     -0.230864   |       -0.229587  | False                 | fail        | permutation/shuffle null control, expected not to reproduce retained branch signature |             0.4 |        0.84     | nan                  |
| no_swarm_model              |    -0.275622  |     -0.305429   |       -0.252852  | True                  | weak        | architectural negative/control comparator                                             |             0.2 |        0.7      | M2_ot_teacher_force  |
| no_teacher_model            |    -0.475492  |     -0.500823   |       -0.450543  | True                  | weak        | architectural negative/control comparator                                             |             0.2 |        0.7      | M1_intrinsic_neural  |
| random_teacher_velocity     |     0.065666  |     -0.0152877  |        0.136086  | True                  | weak        | architectural negative/control comparator                                             |             0.6 |        0.84     | M10_shuffled_time_ot |

## External Experiment E1

Branch nucleation, interpreted as transient condensation-before-divergence, is supported internally under native moscot teacher sensitivity and receives external time-series support in MouseGastrulationData WT chimera sample 1. This remains computational evidence, not experimental validation.

- selected external dataset: MouseGastrulationData WT chimera sample 1.
- external teacher backend: native_moscot.
- external validation tier: acceptable.
- no experimental lineage tracing or experimental validation is claimed; causality and high-impact readiness are not established.
- diffusion remains encoded control-law recovery; birth/death, memory and CCI remain unsupported.
