# DevSpectrum 方法说明

版本：2026-05-06

## Stage-lineage-module time series

输入单细胞 H5AD 后，DevSpectrum 先按：

```text
dataset_id × condition × sample_id × time_point × lineage × module_name
```

聚合模块分数。模块分数使用 log-normalized expression 的 module mean。GSE212050 H5AD 中基因名为 Ensembl ID，因此 DevSpectrum 通过 `GSE212050_feature_metadata_final.txt.gz` 临时补充 gene symbol 映射，不改写原 H5AD。

MVP 默认模块包括：

- mesoderm
- neural
- intermediate
- stress
- cell_cycle
- apoptosis
- haematoendothelial
- erythroid
- endothelial
- cardiac
- extraembryonic_mesoderm
- paraxial_mesoderm
- wnt_response

同时可选输出 SVD latent dimensions。

## Cosine/DCT basis

MVP 不做普通 FFT。对短时间序列使用连续 cosine basis：

```text
y(t) = c0 + sum_k c_k cos(pi k t_scaled)
```

输出 low/high-frequency energy、spectral entropy、dominant basis、phase proxy、amplitude 和 reconstruction error。

## Haar wavelet

Haar wavelet 只作为局部异常和 transient burst 指标，不解释为周期频谱。输出 wavelet low/high energy、transient burst score 和 dominant local index。

## Missing-stage reconstruction

对每条正常 trajectory 做 leave-one-timepoint-out。比较：

- mean baseline
- linear interpolation
- spline
- DCT
- wavelet
- DCT + wavelet

同时做 shuffle-time negative control。若 spectral method 不优于简单 baseline，结果报告必须诚实标注。

## Chimera endpoint projection

Tal1/T 当前只有 E8.5 endpoint，因此 DevSpectrum 只做 endpoint spectral projection。每个 perturbation lineage/module/class 的 observed endpoint module score 与 integrated heldout control endpoint mean 比较：

```text
spectral_residual = observed_perturbation - expected_matched_control
```

并输出 raw residual、cell-weighted residual、DevGuard failure fraction 和 DevGuard-linked spectral failure burden。
