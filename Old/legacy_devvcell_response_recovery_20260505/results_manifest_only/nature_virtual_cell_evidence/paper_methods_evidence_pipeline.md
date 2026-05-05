# Paper Methods: DevVCell Formal Evidence Pipeline

The DevVCell evidence pipeline was run with configuration `config/nature_virtual_cell_evidence.json`. The pipeline validates all required input files, checks table schemas, verifies required primary and baseline models, computes paired statistical comparisons, records full provenance, and writes a manifest linking every claim to output tables and figures.

## Input Validation

For each required result table, the pipeline checks file existence, SHA-256 hash, minimum row count, required columns, and numeric validity for required metrics. Required numeric columns must be present, parseable, non-missing, and finite unless explicitly listed under `schema.nullable_numeric_columns`; nullable columns must still contain at least one finite value and no invalid tokens. Missing inputs, schema errors, missing required models, and invalid required metrics are treated as hard failures according to the configured failure policy.

## Developmental Transition Evidence

Heldout cell-level transition metrics are grouped by model and evaluated using paired system-stage strata. The primary model is `context_residual_mlp`. For each baseline, the pipeline computes mean paired MSE difference, bootstrap 95% confidence interval with 10000 resamples, relative improvement, the number of strata in which the primary model is better, and a paired Wilcoxon signed-rank p-value with Benjamini-Hochberg correction.

## Competence-Window Evidence

TF/GRN stimulus-response rows are stratified by the configured competence window, TS14--TS19. For response amplitude, fate displacement, recovery cost, inverse recovery probability, and alignment with normal transition, the pipeline reports group means, medians, bootstrap confidence intervals, window-versus-outside contrasts, Mann-Whitney tests, and false-discovery-rate adjusted q-values. It also tests Pearson and Spearman correlations between stage vulnerability and stimulus-response metrics.

## Fate and Recovery Virtual Screen

The pipeline constructs a TF-system-stage table containing response amplitude, fate displacement, recovery cost, recovery probability, priority score, and rescue candidate annotations. Rescue candidates are carried forward only as computational hypotheses unless supported by external perturbation or wet-lab rescue assays.

## External Perturbation Calibration

External scPerturb guide-transfer metrics are evaluated by condition-paired comparisons. The primary external model is `gene_context_ridge_residual`. The pipeline computes paired improvements against configured baselines, bootstrap confidence intervals, Wilcoxon signed-rank p-values, effect-cosine summaries, and condition-centered metrics that remove condition-level average difficulty before comparing models.

## Claim Gates

Claim gates are recorded in `tables/claim_gate_matrix.csv`. They are designed to make the current evidence boundary explicit. By default, claim-gate failures are reported but do not abort the pipeline unless `fail_on_claim_threshold_failure` is enabled.

## Interpretation Boundary

The current TF/GRN stimulus, fate-displacement, and recovery-cost outputs are computational proxy readouts. They are suitable for formal hypothesis generation and manuscript-facing evidence tracking, but they do not replace direct CellOT/GEARS/scGen/CPA/foundation-model baselines or stage-dependent perturbation experiments.
