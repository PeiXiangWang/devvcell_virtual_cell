# DevSpectrum 结果报告

版本：2026-05-06

## 运行命令

```bash
python scripts/devspectrum/run_devspectrum_pipeline.py --mode quick
python scripts/devspectrum/run_devspectrum_pipeline.py --mode main
```

验证：

```text
python -m compileall src scripts/devspectrum
python -m pytest
```

## 正常发育频谱

GSE212050 strict sample main run 输出：

| item | value |
|---|---:|
| time-series rows | 2,232 |
| condition-level rows | 162 |
| lineages | 3 |
| time points | 5 |
| features/modules | 18 |
| gene modules | 13 |
| latent dimensions | 5 |

所有 13 个 gene modules 在 GSE212050 feature metadata 映射后 coverage 均为 1.0。核心输出：

```text
results/devspectrum/timeseries/stage_lineage_module_timeseries.csv
results/devspectrum/timeseries/module_gene_coverage.csv
results/devspectrum/spectral_fits/spectral_features.csv
```

## 缺失阶段重建

Leave-one-timepoint-out benchmark 显示：

| method | mean MSE | median MSE | mean MAE | shuffle-time mean MSE |
|---|---:|---:|---:|---:|
| linear | 1.3653 | 0.0034 | 0.5073 | 3.1952 |
| wavelet | 1.3653 | 0.0034 | 0.5073 | 3.1952 |
| DCT+wavelet | 1.6267 | 0.0034 | 0.5572 | 3.2061 |
| mean | 2.0761 | 0.0059 | 0.5806 | 2.0761 |
| spline | 2.1018 | 0.0026 | 0.6129 | 7.1443 |
| DCT | 2.3006 | 0.0051 | 0.6548 | 4.5633 |

解释：在这个 5-timepoint control set 中，linear/wavelet 是整体最强 baseline；DCT 在部分 lineage/module 上优于 linear，但不是整体最优。因此 DevSpectrum MVP 应被写成解释层和 benchmarked spectral-profile layer，而不是声称 DCT 是全局最优预测器。

## Tal1/T chimera endpoint projection

Chimera endpoint projection 生成 1,443 个 residual rows。

| perturbation | raw spectral distance | mean abs residual | cell-weighted mean abs residual | DevGuard failure fraction | DevGuard-linked failure burden | top residual module |
|---|---:|---:|---:|---:|---:|---|
| Tal1 chimera | 6.3200 | 0.0896 | 0.0758 | 0.4230 | 0.0321 | erythroid |
| T chimera | 2.8123 | 0.0518 | 0.0337 | 0.0506 | 0.0017 | erythroid |

Tal1 在 raw residual、cell-weighted residual 和 DevGuard-linked spectral failure burden 上均高于 T。Top residual module 为 erythroid，符合 Tal1/Scl 对血液发生必需性的已知生物学。

## DevGuard link

DTI 相关分析输出：

```text
results/devspectrum/devguard_link/devguard_spectral_correlation.csv
results/devspectrum/devguard_link/failure_mode_spectral_signature.csv
results/devspectrum/devguard_link/devguard_spectrum_link_summary.md
```

当前 Tal1 lineage-level DTI 与 endpoint residual 的 Spearman 相关不强，说明 endpoint module residual 与 DevGuard DTI 捕捉的是互补层面；T 的 DTI 与 spectral distance 呈负相关，提示其少量低 DTI lineages 对 endpoint residual 更敏感。

## Spectral rescue candidates

Tal1 in silico spectral rescue ranking 中，erythroid residual 占总 residual 的 67.6%，单独达到 50% correction coverage。后续模块包括 extraembryonic_mesoderm、cell_cycle、endothelial、apoptosis 和 WNT response。

输出：

```text
results/devspectrum/rescue_candidates/spectral_rescue_candidates.csv
results/devspectrum/rescue_candidates/tal1_spectral_rescue_report.md
```

该结果只能写成 in silico spectral rescue hypotheses。

## Figure outputs

```text
results/devspectrum/figures/figure2_low_frequency_energy.png
results/devspectrum/figures/figure2_high_frequency_energy.png
results/devspectrum/figures/figure2_spectral_entropy.png
results/devspectrum/figures/figure3_missing_stage_reconstruction.png
results/devspectrum/figures/figure4_perturbation_spectral_fingerprint.png
```

当前 MVP 已生成 Figure 2-4 草图。Figure 1 概念图、Figure 5/6 精修图可在后续 manuscript pass 中制作。

## 当前限制

- GSE212050 只有 5 个正常 time points，不能使用复杂高频 FFT。
- Tal1/T 是 E8.5 endpoint projection，不是完整 perturbation time spectrum。
- DCT 不应被表述为全局最优 reconstruction method。
- Rescue candidates 是频谱空间假设，不是机制验证。
