# Discovery CCI Branch Bias

CCI hardening adds sender-receiver specificity, shuffle/zero/random LR controls and feature-recomputed perturbation proxies.

Full population graph-perturbation rerollout is not yet implemented; therefore this law is capped at weak/demonstration level.

## Tier

- tier: weak
- effect_size: -0.0290759 [-0.0290759, -0.0290759]
- permutation_q: 0.0625
- seed_stability_pass: True (sign consistency=1.000)
- negative_control_pass: True
- rollout_based: False

|   seed | perturbation         | sender          | receiver        |   branch_probability_shift |   fate_entropy_shift |   birth_hazard_shift |   death_hazard_shift |   diffusion_shift |   branch_composition_shift |   event_count_shift |   sender_receiver_specificity |
|-------:|:---------------------|:----------------|:----------------|---------------------------:|---------------------:|---------------------:|---------------------:|------------------:|---------------------------:|--------------------:|------------------------------:|
|      7 | remove_lr_edge_proxy | neural          | mesoderm_muscle |                 -0.0559294 |           0.00664494 |          -0.00981631 |           0.00654039 |       0.00169357  |                 -0.0559294 |         -0.0163567  |                      0.645833 |
|     17 | remove_lr_edge_proxy | neural          | mesoderm_muscle |                 -0.0559294 |           0.00664494 |          -0.0123717  |           0.00824971 |       0.00169714  |                 -0.0559294 |         -0.0206214  |                      0.645833 |
|     23 | remove_lr_edge_proxy | neural          | mesoderm_muscle |                 -0.0559294 |           0.00664494 |          -0.0108527  |           0.00724016 |       0.00172327  |                 -0.0559294 |         -0.0180929  |                      0.645833 |
|     42 | remove_lr_edge_proxy | neural          | mesoderm_muscle |                 -0.0559294 |           0.00664494 |          -0.00939129 |           0.0062675  |       0.0017208   |                 -0.0559294 |         -0.0156588  |                      0.645833 |
|     99 | remove_lr_edge_proxy | neural          | mesoderm_muscle |                 -0.0559294 |           0.00664494 |          -0.0128503  |           0.00857296 |       0.00171768  |                 -0.0559294 |         -0.0214232  |                      0.645833 |
|      7 | remove_lr_edge_proxy | mesoderm_muscle | mesoderm_muscle |                 -0.0268639 |           0.00319168 |          -0.00471495 |           0.00314147 |       0.00081345  |                 -0.0268639 |         -0.00785641 |                      0.310205 |
|     17 | remove_lr_edge_proxy | mesoderm_muscle | mesoderm_muscle |                 -0.0268639 |           0.00319168 |          -0.00594235 |           0.00396248 |       0.000815168 |                 -0.0268639 |         -0.00990484 |                      0.310205 |
|     23 | remove_lr_edge_proxy | mesoderm_muscle | mesoderm_muscle |                 -0.0268639 |           0.00319168 |          -0.00521276 |           0.00347758 |       0.000827717 |                 -0.0268639 |         -0.00869034 |                      0.310205 |
|     42 | remove_lr_edge_proxy | mesoderm_muscle | mesoderm_muscle |                 -0.0268639 |           0.00319168 |          -0.0045108  |           0.00301039 |       0.000826532 |                 -0.0268639 |         -0.00752119 |                      0.310205 |
|     99 | remove_lr_edge_proxy | mesoderm_muscle | mesoderm_muscle |                 -0.0268639 |           0.00319168 |          -0.00617221 |           0.00411774 |       0.000825032 |                 -0.0268639 |         -0.01029    |                      0.310205 |
|      7 | remove_lr_edge_proxy | mesoderm_muscle | neural          |                 -0.0201858 |           0.00237743 |          -0.00337734 |           0.00224749 |       0.000609651 |                 -0.0201858 |         -0.00562483 |                      0.233333 |
|     17 | remove_lr_edge_proxy | mesoderm_muscle | neural          |                 -0.0201858 |           0.00237743 |          -0.00449643 |           0.00299627 |       0.000610441 |                 -0.0201858 |         -0.0074927  |                      0.233333 |
|     23 | remove_lr_edge_proxy | mesoderm_muscle | neural          |                 -0.0201858 |           0.00237743 |          -0.00368466 |           0.00245666 |       0.000613425 |                 -0.0201858 |         -0.00614132 |                      0.233333 |
|     42 | remove_lr_edge_proxy | mesoderm_muscle | neural          |                 -0.0201858 |           0.00237743 |          -0.00337024 |           0.00224916 |       0.000615696 |                 -0.0201858 |         -0.0056194  |                      0.233333 |
|     99 | remove_lr_edge_proxy | mesoderm_muscle | neural          |                 -0.0201858 |           0.00237743 |          -0.00455725 |           0.00303813 |       0.000618306 |                 -0.0201858 |         -0.00759538 |                      0.233333 |
