# DevGuard 结果报告

版本：2026-05-05

## 当前结论边界

本仓库已经完成 DevGuard 的代码重构、旧版 DevVCell response-recovery 归档、conformal normality MVP、stress fixture、sample-level split、分类歧义输出、DTI bootstrap CI，以及三个真实公开数据入口的本地 smoke test。

目前可以报告的是：

- 软件管线可运行；
- GSE212050 真实 mouse gastruloid control 数据已能构建 normality reference；
- MouseGastrulationData 的 WT/Tal1 chimera 样本已能完成 perturbation classification smoke test；
- GSE123187 RAW tar 已下载并可抽取 preview H5AD；
- 所有真实生物学结论仍需在完整 sample-level 数据接入后重跑。

## Quick Fixture

运行命令：

```bash
python scripts/devguard/run_devguard_pipeline.py --mode quick
```

输出：

```text
data/processed/devguard/devguard_quick_mouse.h5ad
results/devguard/normality_reference/
results/devguard/perturbation_classification/
results/devguard/tolerance_index/
results/devguard/figures/
```

规模：

| 项目 | 数值 |
|---|---:|
| control cells | 1,350 |
| perturbed cells | 120 |
| genes | 90 |
| reference groups | 9 |
| score-method records | 18 |

Heldout control FPR：

| 指标 | 数值 |
|---|---:|
| alpha | 0.05 |
| mean FPR | 0.0019 |
| max FPR | 0.0333 |

五类分类结果：

| class | cells |
|---|---:|
| within_stage_normal | 24 |
| developmental_delay | 24 |
| developmental_acceleration | 24 |
| fate_deviation | 24 |
| abnormal_off_normal | 24 |

Quick fixture 只证明代码路径和分类逻辑正确，不用于论文生物学结论。

## Stress Fixture

运行命令：

```bash
python scripts/devguard/run_devguard_pipeline.py --mode stress
```

Stress fixture 在 quick fixture 之外加入 Poisson/gamma count noise、batch shift、donor shift、cell-to-cell variation、dropout、ambiguous intermediate cells 和 label noise。

Heldout control FPR：

| 指标 | 数值 |
|---|---:|
| alpha | 0.05 |
| mean FPR | 0.0351 |
| max FPR | 0.2162 |

Perturbed-cell 分类：

| class | cells |
|---|---:|
| within_stage_normal | 259 |
| fate_deviation | 19 |
| developmental_delay | 9 |
| developmental_acceleration | 1 |

Stress fixture 的最大 FPR 明显高于 quick fixture，说明干净 synthetic fixture 会低估边界失败风险。后续论文主结果必须以真实 sample-level split 为准。

## GSE212050 真实 Control Reference

已完成：

- 下载 GEO metadata；
- 下载 `GSE212050_seurat_final.rds.gz`；
- 用 SeuratObject 导出 matrix/obs/var；
- 下采样 15,000 cells 转为 DevGuard H5AD；
- 以 `sample_id` 做 sample-level split，构建 control normality reference。

关键文件：

```text
data/processed/devguard/GSE212050_downsample_15000.h5ad
config/devguard/normality_model_gse212050_downsample.json
results/devguard_real/GSE212050_downsample/normality_reference/
```

规模与校准：

| 项目 | 数值 |
|---|---:|
| cells | 15,000 |
| genes | 30,406 |
| reference groups | 21 |
| score-method records | 42 |
| mean heldout FPR | 0.0472 |
| max heldout FPR | 0.3333 |
| low-heldout records | 22 / 42 |

解释：

GSE212050 下采样 control reference 的平均 FPR 接近 alpha=0.05，说明真实数据入口和 conformal calibration 可运行。最大 FPR 来自低 heldout 组，例如 `d4_d4.5 / Intermediate mesoderm` 只有 3 个 heldout cells，因此已在 quality table 中用 `low_heldout_flag` 标出。正式论文应使用完整数据并过滤低 sample-unit 组。

## MouseGastrulationData WT/Tal1 Chimera Smoke Test

已完成：

- 安装 Bioconductor `MouseGastrulationData`；
- 导出 embryo atlas sample 1；
- 导出 WT chimera sample 1；
- 导出 Tal1 chimera sample 1；
- 用 WT chimera sample 1 构建 control reference；
- 用 Tal1 chimera sample 1 做 perturbation classification；
- 输出 DTI 和 bootstrap CI。

关键文件：

```text
data/processed/devguard/MouseGastrulationData_wt_chimera_sample1.h5ad
data/processed/devguard/MouseGastrulationData_tal1_chimera_sample1.h5ad
config/devguard/normality_model_mouse_gastrulation_wt_chimera_sample1.json
config/devguard/perturbation_tests_tal1_chimera_sample1.json
results/devguard_real/MouseGastrulationData_wt_chimera_sample1/normality_reference/
results/devguard_real/MouseGastrulationData_tal1_chimera_sample1/perturbation_classification/
results/devguard_real/MouseGastrulationData_tal1_chimera_sample1/tolerance_index/
```

WT chimera reference：

| 项目 | 数值 |
|---|---:|
| reference groups | 12 |
| score-method records | 24 |
| mean heldout FPR | 0.0294 |
| max heldout FPR | 0.1389 |
| split strategy | cell_fallback |

Tal1 classification smoke test：

| class | cells | fraction |
|---|---:|---:|
| fate_deviation | 9,325 | 0.6322 |
| developmental_delay | 4,862 | 0.3296 |
| abnormal_off_normal | 564 | 0.0382 |

Tie/ambiguity：

| ambiguous_flag | cells | fraction |
|---|---:|---:|
| False | 11,721 | 0.7946 |
| True | 3,030 | 0.2054 |

解释：

这一步只能作为真实 perturbation ingestion smoke test。因为当前只导出一个 WT chimera 样本，`sample_split.strategy=sample` 无法形成独立 sample-level train/cal/test，所以 reference 自动退化为 `cell_fallback`。正式论文不能直接使用这组数值作为 Tal1 生物学发现。

## GSE123187 Spatial/Tomo Preview

已完成：

- 下载 `GSE123187_RAW.tar`；
- 新增 RAW tar 解析脚本；
- 抽取前 4 个 `coutb.tsv.gz` 文件生成 preview H5AD。

关键文件：

```text
scripts/devguard/build_gse123187_h5ad_from_raw_tar.py
data/processed/devguard/GSE123187_preview_4files.h5ad
results/devguard/dataset_metadata/GSE123187_h5ad_manifest.json
```

Preview 规模：

| 项目 | 数值 |
|---|---:|
| cells | 1,536 |
| genes | 24,086 |
| samples | 4 |
| time_point | 5dAA |
| mode | single_gastruloid |

解释：

GSE123187 已经进入 DevGuard schema，但 preview 目前没有 cell type/lineage 注释，因此只能作为 spatial/tomo ingestion 入口，不能直接用于 fate deviation validation。

## 本轮新增的关键工程改动

- `build_normality_groups` 支持 sample-level split，并修复 sample split 的全局 cell 索引问题；
- classification 输出 `assigned_class_by_priority`、`assigned_class_by_max_pvalue`、`pvalue_margin`、`ambiguous_flag`；
- classification 批量化按 reference group 打分，Tal1 14,751 cells 分类由超时变为约 9 秒；
- DTI 输出 bootstrap CI，并优先 sample-level bootstrap；
- quality table 输出 split unit counts 和 `low_heldout_flag`；
- 新增 `download_public_datasets.py`、Seurat RDS exporter、MouseGastrulationData exporter、MatrixMarket-to-H5AD converter、GSE123187 RAW tar parser。

## 验证

```text
python -m pytest
12 passed
```

## 下一步

正式结果需要：

- 导出完整 MouseGastrulationData embryo atlas 和 WT/Tal1/T chimera 数据；
- 对真实数据强制使用 embryo/organoid-level split，过滤低 sample-unit group；
- 对 GSE212050 使用完整 RDS 或更大 downsample 做 leave-one-organoid-out FPR；
- 为 GSE123187 补充 spatial/tomo axis 注释和 cell type/lineage 映射；
- 所有 perturbation 结论必须和 heldout control FPR、organoid heterogeneity FPR 对比后再写入 manuscript。
