# Local Perturbation Propagation v1.3

- conclusion: localized_response_without_topological_specificity
- Interpretation: perturbation propagation remains an in silico graph diagnostic. It should not be described as experimental intervention evidence.

| dataset   | graph_rule   | model_proxy     |   perturbation_strength |   affected_fraction |   response_attenuation |   mean_seed_velocity_norm | localized   |
|:----------|:-------------|:----------------|------------------------:|--------------------:|-----------------------:|--------------------------:|:------------|
| internal  | topological  | M5_graph_swarm  |                    0.05 |           0.0938889 |             0.00457531 |                   7.33179 | True        |
| internal  | topological  | M5_graph_swarm  |                    0.1  |           0.0938889 |             0.00915063 |                   7.33179 | True        |
| internal  | topological  | M5_graph_swarm  |                    0.2  |           0.0938889 |             0.0183013  |                   7.33179 | True        |
| internal  | topological  | M2_teacher_only |                    0.05 |           0.0111111 |             5e-06      |                   7.33179 | True        |
| internal  | topological  | M2_teacher_only |                    0.1  |           0.0111111 |             1e-05      |                   7.33179 | True        |
| internal  | topological  | M2_teacher_only |                    0.2  |           0.0111111 |             2e-05      |                   7.33179 | True        |
| internal  | metric       | M5_graph_swarm  |                    0.05 |           0.0938889 |             0.00457531 |                   7.33179 | True        |
| internal  | metric       | M5_graph_swarm  |                    0.1  |           0.0938889 |             0.00915063 |                   7.33179 | True        |
| internal  | metric       | M5_graph_swarm  |                    0.2  |           0.0938889 |             0.0183013  |                   7.33179 | True        |
| internal  | metric       | M2_teacher_only |                    0.05 |           0.0111111 |             5e-06      |                   7.33179 | True        |
| internal  | metric       | M2_teacher_only |                    0.1  |           0.0111111 |             1e-05      |                   7.33179 | True        |
| internal  | metric       | M2_teacher_only |                    0.2  |           0.0111111 |             2e-05      |                   7.33179 | True        |
| internal  | random       | M5_graph_swarm  |                    0.05 |           0.256111  |             0.00457531 |                   7.33179 | True        |
| internal  | random       | M5_graph_swarm  |                    0.1  |           0.256111  |             0.00915063 |                   7.33179 | True        |
| internal  | random       | M5_graph_swarm  |                    0.2  |           0.256111  |             0.0183013  |                   7.33179 | True        |
| internal  | random       | M2_teacher_only |                    0.05 |           0.0111111 |             5e-06      |                   7.33179 | True        |
| internal  | random       | M2_teacher_only |                    0.1  |           0.0111111 |             1e-05      |                   7.33179 | True        |
| internal  | random       | M2_teacher_only |                    0.2  |           0.0111111 |             2e-05      |                   7.33179 | True        |
| E1        | topological  | M5_graph_swarm  |                    0.05 |           0.0760095 |             0.00457531 |                   8.23924 | True        |
| E1        | topological  | M5_graph_swarm  |                    0.1  |           0.0760095 |             0.00915063 |                   8.23924 | True        |
| E1        | topological  | M5_graph_swarm  |                    0.2  |           0.0760095 |             0.0183013  |                   8.23924 | True        |
| E1        | topological  | M2_teacher_only |                    0.05 |           0.0158353 |             5e-06      |                   8.23924 | True        |
| E1        | topological  | M2_teacher_only |                    0.1  |           0.0158353 |             1e-05      |                   8.23924 | True        |
| E1        | topological  | M2_teacher_only |                    0.2  |           0.0158353 |             2e-05      |                   8.23924 | True        |
