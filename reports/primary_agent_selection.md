# Primary Agent Selection

Full model is not automatically the primary model. Unsupported modules are excluded from retained main claims.
The primary model is selected by teacher fidelity plus branch-nucleation evidence, not by architectural completeness.
Architectural controls can retain related condensation signals; therefore primary selection identifies the minimal retained mechanistic model, not a proven causal necessity claim.

- primary_mechanistic_model: M5_ot_swarm
- reason: best fidelity/mechanism score after penalizing unsupported modules; unsupported burden=0.

| model                                 | teacher_fidelity_tier   |   relative_sinkhorn |   relative_mmd |   composition_rmse | branch_nucleation_tier   |   branch_nucleation_effect | branch_nucleation_seed_stability   |   unsupported_module_burden |   complexity_penalty | uses_unsupported_modules   |   selection_score |   mean_sinkhorn | recommendation                     |
|:--------------------------------------|:------------------------|--------------------:|---------------:|-------------------:|:-------------------------|---------------------------:|:-----------------------------------|----------------------------:|---------------------:|:---------------------------|------------------:|----------------:|:-----------------------------------|
| M5_ot_swarm                           | acceptable              |             1.13174 |        1.56129 |         0.00402532 | strong                   |                  -0.265822 | True                               |                           0 |                    0 | False                      |           13.7326 |        0.288679 | primary_mechanistic_model          |
| M9_full_memory                        | acceptable              |             1.20363 |        1.62257 |         0.0217414  | strong                   |                  -0.207681 | True                               |                           1 |                    1 | True                       |           11.8403 |        0.307017 | secondary_exploratory_model        |
| M7_ot_swarm_birth_death_diffusion     | acceptable              |             1.20497 |        1.59795 |         0.0219325  | strong                   |                  -0.21572  | True                               |                           1 |                    1 | True                       |           11.8397 |        0.307358 | models_not_retained_for_main_claim |
| M8_ot_swarm_birth_death_diffusion_cci | acceptable              |             1.20498 |        1.59814 |         0.0219325  | strong                   |                  -0.215772 | True                               |                           2 |                    2 | True                       |           10.3396 |        0.307362 | models_not_retained_for_main_claim |
