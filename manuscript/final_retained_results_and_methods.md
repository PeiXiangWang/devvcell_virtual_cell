# Final Retained Results and Methods

## Central Claim

OT gives the developmental map; SwarmLineage-OT learns microscopic finite-agent rules that realize the map and reveals a rollout-supported branch-nucleation order-parameter signature.

`M0b_ot_interpolation` is an oracle-like OT teacher/reference interpolation. The finite-agent model is evaluated by teacher fidelity, emergent-law robustness and mechanistic usefulness, not by beating the OT reference.

## Native Teacher

| teacher_backend   | native_teacher_available   | external_teacher_validation   | native_teacher_claims_allowed   | strong_biological_claims_allowed   | nature_level_claim_allowed   |   native_temporalproblem_pairs | native_requirements_file                       | status                                                       |
|:------------------|:---------------------------|:------------------------------|:--------------------------------|:-----------------------------------|:-----------------------------|-------------------------------:|:-----------------------------------------------|:-------------------------------------------------------------|
| native_moscot     | True                       | False                         | True                            | False                              | False                        |                              6 | reproducibility/native_moscot_requirements.txt | native teacher available; external validation still required |

## Primary Mechanistic Model

- primary_model: M5_ot_swarm
- full model is not automatically the primary model.
- unsupported modules are excluded from retained main claims.
- architectural controls can show related condensation signals, so module necessity is not claimed.

| model                                 | teacher_fidelity_tier   |   relative_sinkhorn |   relative_mmd |   composition_rmse | branch_nucleation_tier   |   branch_nucleation_effect | branch_nucleation_seed_stability   |   unsupported_module_burden |   complexity_penalty | uses_unsupported_modules   |   selection_score |   mean_sinkhorn | recommendation                     |
|:--------------------------------------|:------------------------|--------------------:|---------------:|-------------------:|:-------------------------|---------------------------:|:-----------------------------------|----------------------------:|---------------------:|:---------------------------|------------------:|----------------:|:-----------------------------------|
| M5_ot_swarm                           | acceptable              |             1.13174 |        1.56129 |         0.00402532 | strong                   |                  -0.265822 | True                               |                           0 |                    0 | False                      |           13.7326 |        0.288679 | primary_mechanistic_model          |
| M9_full_memory                        | acceptable              |             1.20363 |        1.62257 |         0.0217414  | strong                   |                  -0.207681 | True                               |                           1 |                    1 | True                       |           11.8403 |        0.307017 | secondary_exploratory_model        |
| M7_ot_swarm_birth_death_diffusion     | acceptable              |             1.20497 |        1.59795 |         0.0219325  | strong                   |                  -0.21572  | True                               |                           1 |                    1 | True                       |           11.8397 |        0.307358 | models_not_retained_for_main_claim |
| M8_ot_swarm_birth_death_diffusion_cci | acceptable              |             1.20498 |        1.59814 |         0.0219325  | strong                   |                  -0.215772 | True                               |                           2 |                    2 | True                       |           10.3396 |        0.307362 | models_not_retained_for_main_claim |

## Retained Computational Hypotheses

- branch_nucleation: strong tier; interpretation=retained_computational_hypothesis; rollout_based=True; best mechanistic reading is a transient condensation-before-divergence order-parameter signature.
- diffusion: acceptable but encoded_control_law_recovery; not an independent discovery.

## Unsupported Modules

- birth/death, memory hysteresis and CCI branch bias are unsupported under current evidence and excluded from main claims.

## External Validation

External validation has been initiated through a public dataset registry but remains pending.

| dataset                                       | accession_or_url                                                 | doi_or_reference                                              | public_availability   | expression_matrix_availability           | metadata_availability                                  | time_stage_column_availability   | cell_type_fate_lineage_label_availability                | attempt_status                            | usable_for_main_validation      | blocker                                                                                             |
|:----------------------------------------------|:-----------------------------------------------------------------|:--------------------------------------------------------------|:----------------------|:-----------------------------------------|:-------------------------------------------------------|:---------------------------------|:---------------------------------------------------------|:------------------------------------------|:--------------------------------|:----------------------------------------------------------------------------------------------------|
| Waddington-OT iPSC reprogramming              | https://broadinstitute.github.io/wot/tutorial/                   | Schiebinger et al., Cell 2019, DOI:10.1016/j.cell.2019.01.006 | True                  | tutorial input data link provided        | cell_days/time metadata described                      | True                             | cell sets/signatures; no direct lineage tracing evidence | registry_confirmed_download_not_completed | potential                       | Google Drive/Terra download requires manual or network-enabled retrieval; not ingested in this run. |
| scLTdb lineage tracing datasets               | https://scltdb.com/scLT/ and https://zenodo.org/records/12176634 | scLTdb public database; record URL verified                   | True                  | h5ad/rds reported by database            | time, celltype and barcode fields reported by database | True                             | lineage barcodes reported                                | registry_confirmed_download_not_completed | potential_lineage_validation    | Dataset-specific selection/download not automated yet; no matrix was fabricated.                    |
| Tempora sample time-course scRNA-seq datasets | https://baderlab.org/Software/Tempora                            | Tran and Bader, Nucleic Acids Research 2020                   | True                  | supplementary/sample data links reported | time-course metadata expected                          | True                             | cell-type annotations vary by sample                     | registry_confirmed_download_not_completed | feasibility_time_series_support | Not yet converted to project AnnData schema in this run.                                            |

## Limitations

- Native moscot teacher extraction removes the toy-fallback blocker, but native downsample sensitivity and external validation remain incomplete.
- Experimental lineage tracing, wet-lab validation, causal proof and high-impact readiness are not claimed.
