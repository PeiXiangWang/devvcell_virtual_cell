# Discovery CCI Branch Bias

CCI hardening now attempts full population simulator rerollouts using `cci_perturbation` configs: remove LR edge, remove receiver, shuffle receiver, random LR and zero CCI signal.

## Tier

- tier: fail
- effect_size: 0 [0, 0]
- permutation_q: 1
- seed_stability_pass: False (sign consistency=1.000)
- negative_control_pass: False
- rollout_based: True
- status: executed

|   seed | perturbation    | sender          | receiver        |   branch_probability_shift |   fate_entropy_shift |   birth_hazard_shift |   death_hazard_shift |   diffusion_shift |   branch_composition_shift |   event_count_shift |   sender_receiver_specificity | rerollout_type                      |
|-------:|:----------------|:----------------|:----------------|---------------------------:|---------------------:|---------------------:|---------------------:|------------------:|---------------------------:|--------------------:|------------------------------:|:------------------------------------|
|      7 | remove_lr_edge  | neural          | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       0.000379171 |                          0 |                   0 |                      0.645833 | full_population_simulator_rerollout |
|      7 | remove_receiver | neural          | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       0.000379171 |                          0 |                   0 |                      0.645833 | full_population_simulator_rerollout |
|     17 | remove_lr_edge  | neural          | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       0           |                          0 |                   0 |                      0.645833 | full_population_simulator_rerollout |
|     17 | remove_receiver | neural          | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       0           |                          0 |                   0 |                      0.645833 | full_population_simulator_rerollout |
|     23 | remove_lr_edge  | neural          | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       0           |                          0 |                   0 |                      0.645833 | full_population_simulator_rerollout |
|     23 | remove_receiver | neural          | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       0           |                          0 |                   0 |                      0.645833 | full_population_simulator_rerollout |
|     42 | remove_lr_edge  | neural          | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       3.75577e-05 |                          0 |                   0 |                      0.645833 | full_population_simulator_rerollout |
|     42 | remove_receiver | neural          | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       3.75577e-05 |                          0 |                   0 |                      0.645833 | full_population_simulator_rerollout |
|     99 | remove_lr_edge  | neural          | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       0.000626405 |                          0 |                   0 |                      0.645833 | full_population_simulator_rerollout |
|     99 | remove_receiver | neural          | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       0.000555784 |                          0 |                   0 |                      0.645833 | full_population_simulator_rerollout |
|      7 | remove_lr_edge  | mesoderm_muscle | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       0.000379171 |                          0 |                   0 |                      0.310205 | full_population_simulator_rerollout |
|      7 | remove_receiver | mesoderm_muscle | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       0.000379171 |                          0 |                   0 |                      0.310205 | full_population_simulator_rerollout |
|     17 | remove_lr_edge  | mesoderm_muscle | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       0           |                          0 |                   0 |                      0.310205 | full_population_simulator_rerollout |
|     17 | remove_receiver | mesoderm_muscle | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       0           |                          0 |                   0 |                      0.310205 | full_population_simulator_rerollout |
|     23 | remove_lr_edge  | mesoderm_muscle | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       0           |                          0 |                   0 |                      0.310205 | full_population_simulator_rerollout |
|     23 | remove_receiver | mesoderm_muscle | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       0           |                          0 |                   0 |                      0.310205 | full_population_simulator_rerollout |
|     42 | remove_lr_edge  | mesoderm_muscle | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       3.75577e-05 |                          0 |                   0 |                      0.310205 | full_population_simulator_rerollout |
|     42 | remove_receiver | mesoderm_muscle | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       3.75577e-05 |                          0 |                   0 |                      0.310205 | full_population_simulator_rerollout |
|     99 | remove_lr_edge  | mesoderm_muscle | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       0.000555784 |                          0 |                   0 |                      0.310205 | full_population_simulator_rerollout |
|     99 | remove_receiver | mesoderm_muscle | mesoderm_muscle |                          0 |                    0 |                    0 |                    0 |       0.000555784 |                          0 |                   0 |                      0.310205 | full_population_simulator_rerollout |
