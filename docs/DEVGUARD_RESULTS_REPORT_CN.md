# DevGuard 结果报告

版本：2026-05-05

## 当前结论边界

本仓库已经完成 DevGuard 的代码重构、旧版 DevVCell response-recovery 归档、conformal normality MVP、stress fixture、sample-level split、分类歧义输出、DTI bootstrap CI，以及真实公开小鼠数据接入。

本轮又完成了 GSE212050 的严格 sample-level rerun：只保留能够形成独立 sample-level train/cal/test 的 stage-lineage 组，不允许 cell-level fallback。

本轮还完成了 MouseGastrulationData chimera 全量干实验：WT chimera 10 个样本、Tal1 chimera 4 个样本、T chimera 14 个非 QC-fail 样本均已导出并转为 DevGuard H5AD。主 perturbation 分析使用 integrated chimera controls，即 WT control cells 加 Tal1/T chimera tomato-negative host controls，按 source-level latent centering 校正数据源偏移后建立 E8.5 strict sample-level normality boundary。

## Quick Fixture

运行命令：

```bash
python scripts/devguard/run_devguard_pipeline.py --mode quick
```

Heldout control FPR：

| 指标 | 数值 |
|---|---:|
| alpha | 0.05 |
| mean FPR | 0.0019 |
| max FPR | 0.0333 |

Quick fixture 只证明代码路径和分类逻辑正确，不用于论文生物学结论。

## Stress Fixture

运行命令：

```bash
python scripts/devguard/run_devguard_pipeline.py --mode stress
```

Stress fixture 加入 Poisson/gamma count noise、batch shift、donor shift、cell-to-cell variation、dropout、ambiguous intermediate cells 和 label noise。

Heldout control FPR：

| 指标 | 数值 |
|---|---:|
| alpha | 0.05 |
| mean FPR | 0.0351 |
| max FPR | 0.2162 |

Stress fixture 的最大 FPR 明显高于 quick fixture，说明干净 synthetic fixture 会低估 boundary failure risk。

## GSE212050 Strict Sample-Level Reference

目的：用真实 mouse gastruloid control 数据完成严格 sample-level split 后重跑，不让同一 organoid/sample 的 cells 同时进入 train、calibration 和 heldout test。

筛选规则：

```text
min_cells_per_group = 200
min_units_per_group = 8
max_cells_per_group = 1500
allow_cell_fallback = false
```

关键文件：

```text
scripts/devguard/select_gse212050_strict_sample_cells.py
data/processed/devguard/GSE212050_strict_sample_13285.h5ad
config/devguard/normality_model_gse212050_strict_sample.json
results/devguard_real/GSE212050_strict_sample/normality_reference/
```

保留的真实 reference groups：

| time_point | lineage | cells | sample units |
|---|---|---:|---:|
| d3 | mesodermal | 1,500 | 11 |
| d3 | neural | 1,500 | 13 |
| d3.5 | neural | 1,373 | 10 |
| d4 | mesodermal | 1,500 | 22 |
| d4.5 | mesodermal | 1,500 | 16 |
| d4.5 | neural | 1,412 | 12 |
| d5 | intermediate | 1,500 | 13 |
| d5 | mesodermal | 1,500 | 15 |
| d5 | neural | 1,500 | 12 |

严格 sample-level split 后的校准结果：

| 项目 | 数值 |
|---|---:|
| cells | 13,285 |
| genes | 30,406 |
| reference groups | 9 |
| score-method records | 18 |
| split_strategy=sample groups | 9 / 9 |
| cell_fallback groups | 0 |
| low-heldout records | 0 / 18 |
| mean heldout FPR | 0.0487 |
| max heldout FPR | 0.0723 |
| min heldout FPR | 0.0238 |

各组 split unit 情况：

| reference_group | train cells | cal cells | test cells | train units | cal units | test units |
|---|---:|---:|---:|---:|---:|---:|
| d3__mesodermal | 818 | 410 | 272 | 6 | 3 | 2 |
| d3__neural | 690 | 347 | 463 | 6 | 3 | 4 |
| d3.5__neural | 623 | 300 | 450 | 5 | 2 | 3 |
| d4__mesodermal | 756 | 412 | 332 | 11 | 6 | 5 |
| d4.5__mesodermal | 751 | 375 | 374 | 8 | 4 | 4 |
| d4.5__neural | 694 | 372 | 346 | 6 | 3 | 3 |
| d5__intermediate | 694 | 343 | 463 | 6 | 3 | 4 |
| d5__mesodermal | 800 | 400 | 300 | 8 | 4 | 3 |
| d5__neural | 750 | 375 | 375 | 6 | 3 | 3 |

解释：

严格 sample-level rerun 的 mean heldout FPR 为 0.0487，接近 alpha=0.05。最大 FPR 从先前 downsample 低样本组的 0.3333 降到 0.0723，且所有记录 `low_heldout_flag=False`。这组结果比早期随机 downsample 更适合写入“真实 control calibration”结果，但仍应在正式论文中用完整数据或更大规模采样复核。

## GSE212050 Earlier Random Downsample

早期 15,000-cell random downsample 已保留为对照：

| 项目 | 数值 |
|---|---:|
| reference groups | 21 |
| score-method records | 42 |
| mean heldout FPR | 0.0472 |
| max heldout FPR | 0.3333 |
| low-heldout records | 22 / 42 |

该结果说明随机 downsample 会产生低 sample-unit 组，不能作为主结果。

## MouseGastrulationData Integrated Chimera Main Analysis

已完成：

- 安装 Bioconductor `MouseGastrulationData`；
- 导出 WT chimera full：30,703 cells，10 samples；
- 导出 Tal1 chimera full：56,122 cells，4 samples；
- 导出 T chimera full：36,931 cells，14 non-QC-fail samples；
- 合并 WT controls、Tal1 tomato-negative controls 和 T tomato-negative controls，得到 76,514 control cells；
- 用 source-level latent centering 校正 WT/Tal1/T 数据源差异；
- 在 E8.5 组构建 strict sample-level reference；
- 分类 Tal1/T tomato-positive E8.5 perturbation cells；
- 输出 DTI、bootstrap CI 和 perturbation-vs-heldout-control Fisher/FDR。

关键文件：

```text
data/processed/devguard/MouseGastrulationData_wt_chimera_full.h5ad
data/processed/devguard/MouseGastrulationData_tal1_chimera_full.h5ad
data/processed/devguard/MouseGastrulationData_t_chimera_full.h5ad
data/processed/devguard/MouseGastrulationData_integrated_chimera_controls.h5ad
config/devguard/normality_model_mouse_gastrulation_integrated_chimera_controls_e85_strict.json
config/devguard/perturbation_tests_tal1_chimera_full_integrated_e85_strict.json
config/devguard/perturbation_tests_t_chimera_full_integrated_e85_strict.json
```

Integrated chimera control reference：

| 项目 | 数值 |
|---|---:|
| control cells | 76,514 |
| E8.5 reference groups | 28 |
| score method | KNN distance |
| k | 30 |
| split strategy | sample-level only |
| cell fallback groups | 0 |
| mean heldout FPR | 0.0547 |
| max heldout FPR | 0.1946 |

解释：

WT-only strict reference 的平均 FPR 为 0.149，不适合作主结果。Integrated matched-control reference 通过加入 Tal1/T host controls 和 source-level centering，将 KNN mean FPR 降到 0.0547，接近 alpha=0.05。少数 lineages 仍存在高 heldout FPR，应在解释 perturbation lineage-level 结果时作为 calibration-risk flag。

Heldout control baseline：

| class | cells | fraction |
|---|---:|---:|
| within_stage_normal | 19,061 | 0.9548 |
| fate_deviation | 675 | 0.0338 |
| abnormal_off_normal | 228 | 0.0114 |

Tal1 E8.5 tomato-positive classification：

| class | cells | fraction |
|---|---:|---:|
| within_stage_normal | 16,297 | 0.5758 |
| fate_deviation | 7,415 | 0.2620 |
| abnormal_off_normal | 4,593 | 0.1623 |

Tal1 vs heldout control enrichment：

| class | perturbation fraction | control fraction | odds ratio | FDR |
|---|---:|---:|---:|---:|
| abnormal_off_normal | 0.1623 | 0.0114 | 16.77 | 0 |
| fate_deviation | 0.2620 | 0.0338 | 10.14 | 0 |
| within_stage_normal | 0.5758 | 0.9548 | 0.064 | 0 |

T E8.5 tomato-positive classification：

| class | cells | fraction |
|---|---:|---:|
| within_stage_normal | 14,922 | 0.9434 |
| fate_deviation | 674 | 0.0426 |
| abnormal_off_normal | 221 | 0.0140 |

T vs heldout control enrichment：

| class | perturbation fraction | control fraction | odds ratio | FDR |
|---|---:|---:|---:|---:|
| abnormal_off_normal | 0.0140 | 0.0114 | 1.23 | 0.0353 |
| fate_deviation | 0.0426 | 0.0338 | 1.27 | 0.0000248 |
| within_stage_normal | 0.9434 | 0.9548 | 0.790 | 0.00000365 |

Interpretation：

Tal1 chimera 显示强烈的 off-normal 和 fate-deviation enrichment，远高于 heldout control baseline。T chimera 在 E8.5 以 normal retention 为主，仅有小幅但统计显著的 fate-deviation/abnormal enrichment。这提供了一个更强的论文叙事：DevGuard 不只是检测“是否异常”，还可区分强扰动 Tal1 与相对耐受的 T/Brachyury chimera 细胞状态。

Tal1 DTI 中最脆弱的大细胞数 lineages 包括 Allantois、ExE mesoderm、Cardiomyocytes 和 Paraxial mesoderm；高保留 lineages 包括 Spinal cord、Mesenchyme、Surface ectoderm、NMP 和 Haematoendothelial progenitors。T 的大多数 E8.5 lineages DTI 为正，整体更接近 control normality。

## MouseGastrulationData Early Smoke Test

早期 sample 1 smoke test 仍保留为开发记录：

| 项目 | 数值 |
|---|---:|
| WT sample 1 reference groups | 12 |
| WT sample 1 mean heldout FPR | 0.0294 |
| Tal1 sample 1 perturbed cells | 14,751 |

这组早期结果已不再作为主分析使用。

Tal1 classification smoke test：

| class | cells | fraction |
|---|---:|---:|
| fate_deviation | 9,325 | 0.6322 |
| developmental_delay | 4,862 | 0.3296 |
| abnormal_off_normal | 564 | 0.0382 |

这一步只能作为真实 perturbation ingestion smoke test。因为当前只导出一个 WT chimera 样本，无法形成独立 embryo-level split，不能作为 Tal1 生物学结论。

## GSE123187 Spatial/Tomo Preview

已完成：

- 下载 `GSE123187_RAW.tar`；
- 新增 RAW tar 解析脚本；
- 抽取前 4 个 `coutb.tsv.gz` 文件生成 preview H5AD。

Preview 规模：

| 项目 | 数值 |
|---|---:|
| cells | 1,536 |
| genes | 24,086 |
| samples | 4 |
| time_point | 5dAA |
| mode | single_gastruloid |

GSE123187 已经进入 DevGuard schema，但 preview 目前没有 cell type/lineage 注释，因此只能作为 spatial/tomo ingestion 入口。

## 本轮新增工程改动

- `build_normality_groups` 支持 `min_units_per_group`，可强制过滤低 sample-unit 组；
- 新增 `select_gse212050_strict_sample_cells.py`，把严格 sample-level 选择流程脚本化；
- `export_seurat_rds_components.R` 支持 cell-list subset；
- 新增 full chimera export/convert/classification configs；
- 新增 `combine_control_h5ads.py` 合并 matched chimera controls；
- 新增 source-level latent centering，用 Tal1/T host controls 校正数据源偏移；
- 新增 `classify_reference_test_cells.py` 输出 heldout control baseline；
- 新增 `compare_perturbation_to_control_baseline.py` 输出 Fisher exact test 和 BH FDR；
- quality table 输出 split unit counts、`low_heldout_flag`、`calibration_excess_fpr` 和 `high_fpr_flag`。

## 验证

```text
python -m pytest
14 passed
python -m compileall src scripts/devguard
```

LaTeX 草稿：

```text
xelatex -interaction=nonstopmode -halt-on-error main_cn.tex
```

已生成 `manuscript/main_cn.pdf`。

## 下一步

仍可继续增强的方向：

- 对 GSE212050 做 leave-one-organoid-out false-positive calibration；
- 为 GSE123187 补充 spatial/tomo axis 注释和 lineage mapping；
- 对 high-FPR lineages 做 sensitivity analysis 或从主 lineage-level 解释中降权；
- 追加更多公开小鼠 perturbation 数据集，扩大扰动类型覆盖。
