# Discovery CCI Branch Bias

LR edge knockout/shuffle is evaluated by branch probability, fate entropy, birth hazard, diffusion and sender-receiver specificity shifts.

| perturbation        | sender          | receiver        |   branch_probability_shift |   fate_entropy_shift |   birth_hazard_shift |   diffusion_shift |   sender_receiver_specificity |
|:--------------------|:----------------|:----------------|---------------------------:|---------------------:|---------------------:|------------------:|------------------------------:|
| lr_edge_knockout    | neural          | mesoderm_muscle |                -0.0559294  |           0.0322917  |          -0.0258333  |        0.019375   |                      0.645833 |
| lr_edge_knockout    | mesoderm_muscle | mesoderm_muscle |                -0.0268639  |           0.0155103  |          -0.0124082  |        0.00930615 |                      0.310205 |
| lr_edge_knockout    | mesoderm_muscle | neural          |                -0.0201858  |           0.0116667  |          -0.00933333 |        0.007      |                      0.233333 |
| lr_edge_knockout    | neural          | neural          |                -0.0133247  |           0.0077012  |          -0.00616096 |        0.00462072 |                      0.154024 |
| lr_receiver_shuffle | all             | shuffled        |                -0.00671698 |           0.00335849 |          -0.00335849 |        0.00167924 |                      0.21627  |

- cci_branch_bias_gate: True
