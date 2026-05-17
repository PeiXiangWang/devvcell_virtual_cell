# Breakthrough Sprint Evidence Map

Final selected story: **computational branch-window framework with taxonomy and failure-boundary audit**.

No new biological mechanism is established; the defensible paper story is an OT-guided virtual-cell framework that classifies developmental transition windows and reports strict evidence boundaries.

## Direction Tiers

| direction                              | hypothesis                                                                               | success_criteria                                                              | failure_criteria                                         | tier       | interpretation                                                                                  | enters_main_story   |
|:---------------------------------------|:-----------------------------------------------------------------------------------------|:------------------------------------------------------------------------------|:---------------------------------------------------------|:-----------|:------------------------------------------------------------------------------------------------|:--------------------|
| branch_window_taxonomy                 | Branch windows are better treated as a taxonomy rather than a single condensation claim. | internal/E1 typed cleanly and failures become interpretable boundary types    | controls produce same types or taxonomy is arbitrary     | acceptable | Supported as a useful reframing; independent support remains weak.                              | True                |
| grn_regulon_transition_taxonomy        | PathwayFinder-prior regulon features annotate branch-window types.                       | known TF programs add interpretable labels with clean controls                | random TF/regulon controls reproduce the same scores     | fail       | GRN mechanism fails; regulon annotation remains candidate-only.                                 | False               |
| known_developmental_biology_recovery   | The framework can recover established developmental TF programs.                         | multiple known programs recovered across internal/E1/external stress datasets | no known programs or random TFs explain equally          | acceptable | Known TF program recovery is the strongest regulatory use case, but not a novel mechanism.      | False               |
| ot_high_entropy_uncertainty_atlas      | Branch windows may be uncertainty-gated rather than condensation-only.                   | internal/E1 and one independent dataset show high uncertainty windows         | uncertainty is absent or indistinguishable from controls | fail       | Uncertainty is plausible but not cross-dataset strong.                                          | False               |
| agent_value_over_teacher_only          | Agent rollout adds order-parameter interpretability over teacher-only maps.              | M5 preserves fidelity while producing window taxonomy/order parameters        | teacher-only explains all outputs equally                | acceptable | Best method contribution: executable rollout exposes order-parameter audit, not OT superiority. | True                |
| failure_boundary_atlas                 | A rigorous boundary map is a publishable contribution.                                   | failures are categorized by data/method/biology limits                        | failures are untracked or selective                      | acceptable | Strong practical contribution for reviewer defense.                                             | True                |
| detector_robustness_calibration        | The detector can be calibrated across datasets and controls.                             | low false positives and baselines do not fully explain events                 | controls/baselines trigger similarly                     | fail       | Detector useful internally/E1 but not enough for broad independent claim.                       | False               |
| regulatory_perturbation_prioritization | GRN perturbation proxies can prioritize future validation candidates.                    | PathwayFinder-supported TFs overlap internal/E1 with random controls          | single-dataset or random TFs dominate                    | weak       | Candidate TF list is useful for future validation only.                                         | False               |
| spatial_validation_requirement         | Spatial validation can be specified despite current blockers.                            | clear assay/metadata/readout requirements are defined                         | requirements remain vague                                | acceptable | Future requirement is precise; current spatial evidence absent.                                 | False               |
| final_story_selection                  | The final story should focus on the most defensible contribution.                        | one coherent main line with unsupported mechanisms excluded                   | many weak claims compete                                 | acceptable | Select framework + taxonomy + failure-boundary story, not a rescued mechanism claim.            | True                |

## Retained Claim

SwarmLineage-OT is best framed as a native OT-guided finite-agent framework for developmental branch-window order-parameter taxonomy and strict evidence-boundary auditing. The transient condensation-before-divergence signature remains strongest in internal native and E1 MouseGastrulationData.

## Not Retained

- clone fate prediction
- topological-neighbor specificity
- swarm-required causality
- GRN causal mechanism
- birth/death, memory or CCI mechanisms
- diffusion as independent discovery

## Next Exact Experiment

A stage-resolved gastruloid or embryo time-series with matched scRNA-seq, spatial coordinates, curated cell-type/lineage labels and targeted TF perturbation around the predicted branch window. Primary readouts: branch-window taxonomy class, spatial condensation, OT transition entropy, regulon activity divergence and post-window lineage composition.
