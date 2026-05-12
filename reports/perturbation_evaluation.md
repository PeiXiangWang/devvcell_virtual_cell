# Perturbation Evaluation

No matched perturb-seq or drug time-series validation is part of the default SwarmLineage-OT quick dataset. The available scPerturb files are catalogued in the data audit and can support future expression-response benchmarks, but they do not validate developmental lineage dynamics here.

## LR Perturbation Candidates

| ligand   | receptor   |   mean_cci_score |   predicted_entropy_change |   predicted_growth_change |   baseline_ot_entropy | status                |
|:---------|:-----------|-----------------:|---------------------------:|--------------------------:|----------------------:|:----------------------|
| Kitl     | Kit        |           0.0135 |                    0.00135 |                 -0.000675 |              0.619676 | exploratory_in_silico |

## Gene Perturbation Candidates

| gene_or_axis                     | expected_effect                    | evidence                                   |
|:---------------------------------|:-----------------------------------|:-------------------------------------------|
| Mki67/Top2a proliferation module | growth hazard calibration          | cell-cycle marker proxy in current AnnData |
| Fgf8-Fgfr1                       | mesoderm/neural fate-balance shift | curated LR axis if genes are present       |
| Cxcl12-Cxcr4                     | dispersion/niche migration shift   | curated LR axis if genes are present       |
