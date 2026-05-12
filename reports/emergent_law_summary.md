# Emergent Law Summary

The core scientific objective is to learn microscopic finite-agent rules that realize the OT-inferred pseudo-lineage and expose laws that OT interpolation does not directly return.

- passed discovery gates: 6/6
- emergent_law_gate: True

## Discovery Gates

| law               | gate   | table                                   | status   |
|:------------------|:-------|:----------------------------------------|:---------|
| diffusion         | True   | tables\diffusion_law_regression.csv     | executed |
| birth_death       | True   | tables\birth_death_law.csv              | executed |
| branch_nucleation | True   | tables\swarm_order_parameters.csv       | executed |
| memory_hysteresis | True   | tables\memory_hysteresis_experiment.csv | executed |
| cci_branch_bias   | True   | tables\cci_branch_bias.csv              | executed |
| phase_diagram     | True   | tables\phase_diagram.csv                | executed |

## Result Interpretation

- These analyses are mechanistic probes of the trained simulator, not wet-lab validation.
- Stable laws can support the research prototype even when the OT reference has the lowest reconstruction error.
- If both teacher fidelity and all emergent-law gates fail, the current scientific hypothesis is not supported.
