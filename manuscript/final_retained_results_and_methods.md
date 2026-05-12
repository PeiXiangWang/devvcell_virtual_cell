# Final Retained Results and Methods

This document keeps only the current retained prototype results and methods. Installation failures and exploratory noise are kept in `logs/` and `reports/negative_results.md` rather than used as claims.

## Data and Preprocessing

Default input: `data/processed/cell_level_subset_v1.h5ad`, subsampled and processed into `processed/swarmlineage_input.h5ad`. The dataset contains ordered developmental stages and lineage/cell-type labels. Stage is used as ordered developmental time.

## OT Teacher

Adjacent-stage entropic OT couplings are computed through the `run_moscot` entry point. In the current quick prototype, native moscot availability is recorded but an auditable POT/SciPy fallback generates the couplings. Outputs are `processed/ot_teacher.h5ad`, `processed/ot_couplings/*.npz`, and `processed/ot_fate_probabilities.parquet`.

## SwarmLineage-OT Model

A minimal PyTorch velocity model is trained to match OT barycentric descendant vectors. Agent simulations add optional birth/death resampling, adaptive diffusion, local swarm rules, CCI modulation and phenomenological memory.

## Main Quantitative Results

Best mean-rank model: `M0_ot_interpolation`.

| model                                 |   ('sinkhorn', 'mean') |   ('sinkhorn', 'std') |   ('mmd_rbf', 'mean') |   ('mmd_rbf', 'std') |   ('energy', 'mean') |   ('energy', 'std') |   ('knn_two_sample_accuracy', 'mean') |   ('knn_two_sample_accuracy', 'std') |   ('celltype_composition_rmse', 'mean') |   ('celltype_composition_rmse', 'std') |
|:--------------------------------------|-----------------------:|----------------------:|----------------------:|---------------------:|---------------------:|--------------------:|--------------------------------------:|-------------------------------------:|----------------------------------------:|---------------------------------------:|
| M0_ot_interpolation                   |               0.28892  |            0.00936219 |            0.0097388  |          0.00075215  |             0.345831 |           0.0115977 |                              0.571971 |                           0.0204607  |                              0.00381493 |                             0.00262295 |
| M10_shuffled_time_ot                  |               0.316176 |            0.00747013 |            0.00513459 |          0.000904152 |             0.218758 |           0.0329254 |                              0.645606 |                           0.00941169 |                              0.0108743  |                             0.00222133 |
| M11_random_lr_labels                  |               0.301723 |            0.0104156  |            0.0277966  |          0.00378125  |             0.94278  |           0.0765006 |                              0.559145 |                           0.0187483  |                              0.00983581 |                             0.00551719 |
| M1_intrinsic_neural                   |               0.29825  |            0.00876728 |            0.0243012  |          0.00128358  |             0.809875 |           0.0246339 |                              0.561995 |                           0.0144093  |                              0.00408543 |                             0.00133439 |
| M2_ot_teacher_force                   |               0.29825  |            0.00876728 |            0.0243012  |          0.00128358  |             0.809875 |           0.0246339 |                              0.561995 |                           0.0144093  |                              0.00408543 |                             0.00133439 |
| M3_ot_birth_death                     |               0.294263 |            0.00867486 |            0.0257157  |          0.00300539  |             0.849029 |           0.0612029 |                              0.55772  |                           0.018506   |                              0.00905709 |                             0.00460614 |
| M4_ot_adaptive_diffusion              |               0.298341 |            0.00868134 |            0.024264   |          0.00127684  |             0.807488 |           0.0248584 |                              0.562945 |                           0.015485   |                              0.00408543 |                             0.00133439 |
| M5_ot_swarm                           |               0.301489 |            0.00904874 |            0.0267753  |          0.00143666  |             0.909166 |           0.0269617 |                              0.55867  |                           0.0144874  |                              0.00390988 |                             0.00210723 |
| M6_ot_swarm_birth_death               |               0.297042 |            0.00876958 |            0.028327   |          0.00321477  |             0.951595 |           0.0637824 |                              0.556295 |                           0.016405   |                              0.00915502 |                             0.00496029 |
| M7_ot_swarm_birth_death_diffusion     |               0.297119 |            0.00882763 |            0.0282685  |          0.00321027  |             0.948511 |           0.0641559 |                              0.554869 |                           0.0190469  |                              0.00973161 |                             0.00388418 |
| M8_ot_swarm_birth_death_diffusion_cci |               0.297228 |            0.00881111 |            0.0283344  |          0.0032256   |             0.950664 |           0.0645807 |                              0.554869 |                           0.0190469  |                              0.00973161 |                             0.00388418 |
| M9_full_pheromone                     |               0.302977 |            0.0090176  |            0.0309501  |          0.00344266  |             1.05423  |           0.06745   |                              0.55962  |                           0.0161625  |                              0.0096355  |                             0.00391603 |

## Gate Status

Full model passes the quick superiority gate over OT interpolation on at least two core metrics: False.

## Limitations

- Current results are computational evidence only.
- Teacher lineage is OT-inferred pseudo-lineage, not true lineage tracing.
- Native moscot/WOT/TIGON/TrajectoryNet/MIOFlow/CellRank2 baselines remain required for high-impact claims.
- Perturbation predictions are exploratory until matched perturbation time-series or wet-lab validation is performed.

## Reproducibility

Run `bash reproducibility/run_all.sh` on Unix/WSL systems, or `powershell -ExecutionPolicy Bypass -File .\reproducibility\run_all.ps1` on this Windows environment. The PowerShell run was verified locally. Manifest is written to `reproducibility/manifest.json`.
