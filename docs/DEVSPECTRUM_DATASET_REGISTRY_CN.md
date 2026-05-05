# DevSpectrum 数据登记

版本：2026-05-06

| dataset | role | current status |
|---|---|---|
| GSE212050 strict sample | normal developmental time-course | main control spectrum source |
| MouseGastrulationData integrated chimera controls | E8.5 endpoint matched control | used for chimera endpoint expectation |
| MouseGastrulationData Tal1 chimera | perturbation endpoint | projected into endpoint residual space |
| MouseGastrulationData T chimera | perturbation endpoint | projected into endpoint residual space |
| GSE123187 tomo preview | future spatial/tomo validation | not used as main spectral training input |

## GSE212050 control spectrum

当前 main run 使用：

```text
data/processed/devguard/GSE212050_strict_sample_13285.h5ad
```

输出：

```text
results/devspectrum/timeseries/stage_lineage_module_timeseries.csv
results/devspectrum/spectral_fits/spectral_features.csv
results/devspectrum/missing_stage_reconstruction/reconstruction_summary.csv
```

## Chimera endpoint projection

输入：

```text
data/processed/devguard/MouseGastrulationData_integrated_chimera_controls.h5ad
data/processed/devguard/MouseGastrulationData_tal1_chimera_full.h5ad
data/processed/devguard/MouseGastrulationData_t_chimera_full.h5ad
```

DevGuard classification 输入：

```text
results/devguard_real/MouseGastrulationData_integrated_chimera_controls_e85_strict/heldout_control_classification/heldout_control_normality_classes.csv
results/devguard_real/MouseGastrulationData_tal1_chimera_full_integrated_e85_strict/perturbation_classification/cell_normality_classes.csv
results/devguard_real/MouseGastrulationData_t_chimera_full_integrated_e85_strict/perturbation_classification/cell_normality_classes.csv
```
