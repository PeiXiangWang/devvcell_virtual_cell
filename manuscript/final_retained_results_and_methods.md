# Final Retained Results and Methods

## Central Claim

SwarmLineage-OT converts native OT-inferred developmental maps into finite-agent virtual-cell dynamics and currently retains a branch-window order-parameter hypothesis, transient condensation-before-divergence, rather than a clone-fate, CCI, memory, birth/death or topological-specific mechanism claim.

`M0b_ot_interpolation` is an oracle-like OT teacher/reference interpolation. The finite-agent model is evaluated by teacher fidelity, emergent-law robustness and mechanistic usefulness, not by beating the OT reference.

## Tier Summary

| teacher_fidelity_tier   | emergent_law_tier   | mechanistic_usefulness_tier   | mechanistic_gate_pass   | strong_gate   |   laws_at_least_acceptable |   laws_strong | native_or_external_teacher_validation   |
|:------------------------|:--------------------|:------------------------------|:------------------------|:--------------|---------------------------:|--------------:|:----------------------------------------|
| acceptable              | weak                | weak                          | False                   | False         |                          2 |             1 | True                                    |

## Retained Computational Hypotheses

| law               | tier       | interpretation_level              | rollout_based   | directly_supervised_or_encoded   |
|:------------------|:-----------|:----------------------------------|:----------------|:---------------------------------|
| diffusion         | acceptable | encoded_control_law_recovery      | False           | True                             |
| branch_nucleation | strong     | retained_computational_hypothesis | True            | False                            |

## Clone-Aware Evidence Boundary

- latest_audit_source: `outcome_preserving`
- final_clone_aware_status: `weinreb_sampling_specific_condensation_signal`
- final_clone_aware_tier: `weak`
- interpretation: Weinreb retains a sampling-specific primary condensation signal, but Jindal does not recover support after outcome-preserving native sampling; no general clone-aware support is established.
- boundary: Clone-level fate-diversification prediction is not supported under current tested datasets and native sampling strategies.

Retained: branch nucleation / transient condensation-before-divergence as a time-series order-parameter computational hypothesis, supported by internal native moscot and E1 MouseGastrulationData with M5_ot_swarm as the evidence-selected primary model.

Not retained: clone-level fate-diversification prediction from condensation, topological-neighbour-specific mechanism, swarm-required causality, birth/death, memory, CCI, or diffusion as an independent discovery.

Clone-aware analyses are stress tests for the time-series branch signature, not the retained main claim. They do not justify presenting condensation as a clone-level predictor unless the primary outcome is supported across datasets after covariate, matched and negative-control analyses.

## Developmental Time-Series Atlas

- atlas_tier: `weak`
- datasets_attempted: 13
- datasets_analyzed: 3
- acceptable_external_datasets: 0
- analyzed_datasets: E5_CellRank_Farrell_zebrafish_axial_mesoderm, E2_GSE212050_gastruloid_native_atlas, E3_MouseGastrulationData_wt_chimera_full_stage_mapped
- downloaded_new_dataset: E5_CellRank_Farrell_zebrafish_axial_mesoderm
- independent_native_analyzed: E5_CellRank_Farrell_zebrafish_axial_mesoderm, E2_GSE212050_gastruloid_native_atlas, E3_MouseGastrulationData_wt_chimera_full_stage_mapped
- interpretation: generalization beyond internal/E1 remains unresolved after the final independent EB and spatial/time-series sprint
- v2_final_sprint_directions_attempted: 2
- v2_final_sprint_usable_datasets: 1
- v2_final_sprint_acceptable_datasets: 0
- v2_spatial_validation_status: unavailable_with_current_metadata
- final_manuscript_line: internal native plus E1 support retained; independent EB/spatial sprint did not upgrade cross-system support

The atlas is used to define the current external boundary of the branch-window order-parameter hypothesis. Weak or failed atlas rows must not be written as cross-dataset validation. A detected branch-like window without condensation-before-divergence, or with unclean controls, does not upgrade the retained claim.

## GRN / Regulon Evidence Boundary

- final_GRN_tier: `fail`
- known_TF_program_recovery_tier: `acceptable`
- prior_source: PathwayFinder local OmniPath CollecTRI/DoRothEA TF-target gene-symbol priors, with co-expression fallback only for missing prior targets.
- interpretation: GRN features do not currently strengthen the branch-window mechanism beyond expression/OT geometry.
- boundary: GRN/regulon analysis is a computational audit and candidate-generation layer. It does not establish causal GRN control, validated TF perturbation, experimental validation, or a proven regulatory mechanism.

## Breakthrough Sprint Final Story

- selected_story: `computational branch-window framework with taxonomy and failure-boundary audit`
- directions_at_least_acceptable: 6
- strong_directions: 0
- main_story_directions: branch_window_taxonomy, agent_value_over_teacher_only, failure_boundary_atlas, final_story_selection
- interpretation: No new biological mechanism reached strong evidence; the defensible upgrade is branch-window taxonomy, agent rollout interpretability and evidence-boundary auditing.
- boundary: this is a computational taxonomy and hypothesis-generation framework, not a causal biological mechanism claim.

## Communication-Niche Search

- communication_niche_tier: `acceptable`
- conclusion: `communication_niche_priming_is_cross_dataset_candidate`
- analyzed_datasets: 6
- acceptable_datasets: 2
- activation_support_datasets: 2
- receiver_priming_support_datasets: 1
- strongest_module_candidate: `morphogen_patterning`
- module_detail: morphogen_patterning: positive activation in internal_native;E1_MouseGastrulationData;GSE154572_EB_WT
- boundary: this is a candidate extracellular-niche annotation layer. It does not establish confirmed ligand-receptor signalling, communication-driven cause-effect, or experimental perturbation support.

## Strict Morphogen Communication-Niche v2

- final_tier: `weak`
- strongest_module: `FGF_niche`
- support_datasets: E1_MouseGastrulationData;E5_zebrafish_Farrell
- acceptable_datasets: E1_MouseGastrulationData
- strict_extracellular_edges: 28314
- removed_intracellular_or_generic_edges: 1851922
- interpretation: module has cross-dataset direction but limited internal or control support
- boundary: strict family-level morphogen niches are candidate validation targets. Broad communication-niche priming remains more robust than any single confirmed morphogen family.

## Exploratory / Demonstration Only

| law           | tier   | interpretation_level    | rollout_based   | directly_supervised_or_encoded   |
|:--------------|:-------|:------------------------|:----------------|:---------------------------------|
| phase_diagram | weak   | exploratory_sensitivity | True            | False                            |

## Unsupported Claims

| law               | tier   | status   |
|:------------------|:-------|:---------|
| birth_death       | fail   | executed |
| memory_hysteresis | fail   | executed |
| cci_branch_bias   | fail   | executed |

## Core Metrics

| model                                 |   sinkhorn |    mmd_rbf |   energy |   celltype_composition_rmse |
|:--------------------------------------|-----------:|-----------:|---------:|----------------------------:|
| M0_linear_label_interpolation         |   0.281661 | 0.00822039 | 0.292374 |                  0.00396321 |
| M0b_ot_interpolation                  |   0.255076 | 0.00878499 | 0.335883 |                  0.00348103 |
| M10_shuffled_time_ot                  |   0.304366 | 0.00526067 | 0.213518 |                  0.014373   |
| M11_random_lr_labels                  |   0.300051 | 0.0135742  | 0.450108 |                  0.0202327  |
| M1_intrinsic_neural                   |   0.309283 | 0.0277509  | 0.917    |                  0.00420263 |
| M2_ot_teacher_force                   |   0.29088  | 0.0120656  | 0.38817  |                  0.0039636  |
| M3_ot_birth_death                     |   0.307424 | 0.0121     | 0.40647  |                  0.0210994  |
| M4_ot_adaptive_diffusion              |   0.291023 | 0.0120433  | 0.387111 |                  0.00384322 |
| M5_ot_swarm                           |   0.288679 | 0.013716   | 0.439041 |                  0.00402532 |
| M6_ot_swarm_birth_death               |   0.306787 | 0.0132087  | 0.451684 |                  0.0224666  |
| M7_ot_swarm_birth_death_diffusion     |   0.307358 | 0.014038   | 0.46045  |                  0.0219325  |
| M8_ot_swarm_birth_death_diffusion_cci |   0.307362 | 0.0140396  | 0.460488 |                  0.0219325  |
| M9_full_memory                        |   0.307017 | 0.0142543  | 0.468093 |                  0.0217414  |

## Limitations

- Current results are computational hypotheses.
- Some laws are encoded control-law recoveries and must not be written as independent biological discoveries.
- Native moscot teacher extraction removes the toy-fallback blocker for teacher construction, but not the need for external validation.
- No experimental validation or causal mechanism is claimed.
