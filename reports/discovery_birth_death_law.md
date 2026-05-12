# Discovery Birth Death Law

Net growth hazard is analysed against local density, OT mass expansion, cell-cycle score, fate entropy and CCI signal.

## Regression

| response          | predictor        |         coef |    abs_coef |      r2 |    n |
|:------------------|:-----------------|-------------:|------------:|--------:|-----:|
| net_growth_hazard | local_density    | -0.00083333  | 0.00083333  | 0.05934 | 8000 |
| net_growth_hazard | ot_growth        |  0.000759121 | 0.000759121 | 0.05934 | 8000 |
| net_growth_hazard | cell_cycle_score |  1.6263e-18  | 1.6263e-18  | 0.05934 | 8000 |
| net_growth_hazard | fate_entropy     | -0.000464614 | 0.000464614 | 0.05934 | 8000 |
| net_growth_hazard | cci_signal       |  3.99305e-05 | 3.99305e-05 | 0.05934 | 8000 |

- carrying_capacity_like_threshold: True
- birth_death_law_gate: True
- event counts: births=13267, deaths=14917
