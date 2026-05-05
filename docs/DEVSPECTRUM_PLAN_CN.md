# DevSpectrum 实施计划

版本：2026-05-06

DevSpectrum 是 DevGuard 的解释层，而不是 DevGuard 的替代模块。DevGuard 判断扰动细胞是否仍处于正常发育范围；DevSpectrum 分析发育程序随时间展开的趋势、相位和局部异常结构，用于解释 DevGuard-defined normality failure。

## MVP 范围

- 使用 GSE212050 strict sample-level mouse gastruloid control 构建正常 stage-lineage-module time series。
- 对每条 lineage/module trajectory 拟合 DCT/cosine、Haar wavelet 和 spline baseline。
- 做 leave-one-timepoint-out missing-stage reconstruction benchmark，并加入 shuffle-time negative control。
- 将 Tal1/T chimera E8.5 endpoint 投影到 matched-control endpoint module residual space。
- 将 residuals 与 DevGuard class、DTI 和 failure mode 关联。
- 输出 in silico spectral rescue hypotheses。

## 不做的事

- 不声称 Fourier/DCT 发现 GRN。
- 不把 Tal1/T 单时间点 endpoint 写成完整 perturbation time spectrum。
- 不把 rescue candidates 写成真实 rescue 机制。
