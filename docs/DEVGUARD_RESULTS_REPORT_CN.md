# DevGuard 结果报告

版本：2026-05-05

## 当前结论边界

本仓库已经完成 DevGuard 的代码重构、旧版 DevVCell response-recovery 归档、conformal normality MVP、stress fixture、sample-level split、分类歧义输出、DTI bootstrap CI，以及三个真实公开数据入口的本地 smoke test。

本轮又完成了 GSE212050 的严格 sample-level rerun：只保留能够形成独立 sample-level train/cal/test 的 stage-lineage 组，不允许 cell-level fallback。

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

## MouseGastrulationData WT/Tal1 Chimera Smoke Test

已完成：

- 安装 Bioconductor `MouseGastrulationData`；
- 导出 embryo atlas sample 1；
- 导出 WT chimera sample 1；
- 导出 Tal1 chimera sample 1；
- 用 WT chimera sample 1 构建 control reference；
- 用 Tal1 chimera sample 1 做 perturbation classification；
- 输出 DTI 和 bootstrap CI。

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
- 新增 `normality_model_gse212050_strict_sample.json`；
- 重跑 GSE212050 strict sample-level reference，并确认无 cell fallback；
- 保留 quality table 中的 split unit counts 和 `low_heldout_flag`。

## 验证

```text
python -m pytest
13 passed
python -m compileall src scripts/devguard
```

LaTeX 草稿：

```text
xelatex -interaction=nonstopmode -halt-on-error main_cn.tex
```

已生成 `manuscript/main_cn.pdf`。

## 下一步

正式结果仍需要：

- 导出完整 MouseGastrulationData embryo atlas 和 WT/Tal1/T chimera 数据；
- 对真实 perturbation 数据做 embryo-level split 或 matched control reference；
- 对 GSE212050 做 leave-one-organoid-out false-positive calibration；
- 为 GSE123187 补充 spatial/tomo axis 注释和 lineage mapping；
- 所有 perturbation 结论必须和 heldout control FPR、organoid heterogeneity FPR 对比后再写入 manuscript。
