# External Data Integrity Audit

- Registry is not described as validation unless `external_validation_tier` is acceptable or strong.
- External teacher backend is read from the actual generated teacher summary and coupling index.
- Time-series support is not lineage validation because no lineage barcode was present in the selected component.
- Computational branch signature is not described as a biological mechanism.
- Failed or unattempted fallback downloads remain marked as not usable for this run.
- Birth/death, memory and CCI remain unsupported and are not restored to the main claim.
- Primary model remains evidence-selected M5, not M9 by default.

| dataset_id                               | external_validation_tier   | external_teacher_backend   | branch_event_detected   | condensation_before_divergence_reproduced   | lineage_validated   |   lineage_separation_effect |   alignment_effect |   entropy_effect | negative_control_pass   | seed_stability_pass   | interpretation                           | status   | blocker   |
|:-----------------------------------------|:---------------------------|:---------------------------|:------------------------|:--------------------------------------------|:--------------------|----------------------------:|-------------------:|-----------------:|:------------------------|:----------------------|:-----------------------------------------|:---------|:----------|
| E1_mouse_gastrulation_wt_chimera_sample1 | acceptable                 | native_moscot              | True                    | True                                        | False               |                    -1.12489 |          0.0195907 |                0 | True                    | True                  | transient_condensation_before_divergence | analyzed |           |