# Morphogen Communication-Niche v2

## Strict Prior Audit

| metric                                     |   count |
|:-------------------------------------------|--------:|
| total_omnipath_joined_intercell_edges      | 1880236 |
| strict_extracellular_edges                 |   28314 |
| removed_intracellular_or_generic_edges     | 1851922 |
| suspicious_intercell_genes_flagged         |     223 |
| FGF_niche_strict_edges                     |      93 |
| WNT_niche_strict_edges                     |     143 |
| BMP_niche_strict_edges                     |      95 |
| TGF_NODAL_ACTIVIN_niche_strict_edges       |      33 |
| SHH_niche_strict_edges                     |      13 |
| Notch_Delta_niche_strict_edges             |      20 |
| ECM_adhesion_guidance_niche_strict_edges   |     346 |
| chemokine_growth_factor_niche_strict_edges |     266 |

## Cross-Dataset Module Tiers

| module                        | tier   |   support_dataset_count | support_datasets                                   | acceptable_datasets      | has_internal   | has_e1   | has_independent   |   mean_activation_effect |   mean_post_event_divergence_effect | interpretation                                                             |
|:------------------------------|:-------|------------------------:|:---------------------------------------------------|:-------------------------|:---------------|:---------|:------------------|-------------------------:|------------------------------------:|:---------------------------------------------------------------------------|
| SHH_niche                     | fail   |                       0 |                                                    |                          | False          | False    | False             |              -0.00198988 |                         -0.00165209 | no reliable module-specific branch-window niche signal                     |
| chemokine_growth_factor_niche | fail   |                       0 |                                                    |                          | False          | False    | False             |              -0.00446206 |                          0.00273673 | no reliable module-specific branch-window niche signal                     |
| BMP_niche                     | weak   |                       1 | E5_zebrafish_Farrell                               |                          | False          | False    | True              |              -0.0308252  |                          0.00375527 | dataset-specific candidate only                                            |
| ECM_adhesion_guidance_niche   | weak   |                       1 | E1_MouseGastrulationData                           |                          | False          | True     | False             |              -0.0246783  |                          0.00462032 | dataset-specific candidate only                                            |
| FGF_niche                     | weak   |                       2 | E1_MouseGastrulationData;E5_zebrafish_Farrell      | E1_MouseGastrulationData | False          | True     | True              |              -0.0139858  |                          0.00201757 | module has cross-dataset direction but limited internal or control support |
| Notch_Delta_niche             | weak   |                       1 | GSE154572_EB_WT                                    |                          | False          | False    | True              |              -0.00300377 |                          0.00106104 | dataset-specific candidate only                                            |
| TGF_NODAL_ACTIVIN_niche       | weak   |                       2 | E3_MouseGastrulationData_full;E5_zebrafish_Farrell |                          | False          | False    | True              |               0.00109626 |                          0.0027895  | dataset-specific candidate only                                            |
| WNT_niche                     | weak   |                       1 | GSE154572_EB_WT                                    | GSE154572_EB_WT          | False          | False    | True              |              -0.00413386 |                          0.0031672  | dataset-specific candidate only                                            |

## Final Interpretation

- final_communication_niche_tier: `weak`
- strongest_module: `FGF_niche`
- allowed_claim: strict extracellular morphogen/communication-niche priming is a candidate branch-window annotation.
- forbidden_claim: do not describe this as confirmed signalling, communication-driven fate control, or experimental validation.

## Future Spatial/Perturbation Design

Primary hypothesis: morphogen communication-niche priming occurs at the branch window and can be spatially observed as local ligand-producing sender neighborhoods around receptor-competent receiver cells.
Required assay: spatial transcriptomics, MERFISH, seqFISH or smFISH in a gastruloid/embryoid-body time course with branch-window stages, ligand/receptor readouts, cell type labels and optional FGF/WNT/BMP/TGF/NODAL perturbations.
Primary readouts: sender ligand density, receiver receptor competence, local niche field, branch-window order parameter and post-event lineage divergence.
