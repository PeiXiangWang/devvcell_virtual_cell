# DevGuard 数据注册表

版本：2026-05-05

## 原则

DevGuard 主线只使用独立公开小鼠 embryo/gastruloid 数据。`scLine_pro.h5ad`、RDEG nodes、RDEG OT、GRN rollout 和 matrix writeback 不进入新版主图证据链。

配置文件：

```text
config/devguard/datasets_mouse.json
```

## 当前注册数据

| dataset_id | 角色 | 当前状态 |
|---|---|---|
| `E-MTAB-6967` | main normality reference | Bioconductor 包已安装，embryo atlas sample 1 已导出 |
| `MouseGastrulationData_WTChimera` | chimera control reference | WT chimera sample 1 已导出并完成 reference smoke test |
| `MouseGastrulationData_Tal1Chimera` | real mouse perturbation | Tal1 chimera sample 1 已导出并完成 classification smoke test |
| `MouseGastrulationData_TChimera` | real mouse perturbation | 已注册，尚未导出 |
| `GSE212050` | organoid heterogeneity negative control | Seurat RDS 已下载，strict sample-level H5AD 已构建，control reference 已重跑 |
| `GSE123187` | spatial/tomo validation | RAW tar 已下载，4-file preview H5AD 已构建 |
| `devguard_quick_mouse` | software fixture | quick mode 可生成 |
| `devguard_stress_mouse` | software stress fixture | stress mode 可生成 |

## 已验证的本地入口

```text
data/processed/devguard/GSE212050_strict_sample_13285.h5ad
data/processed/devguard/GSE212050_downsample_15000.h5ad
data/processed/devguard/MouseGastrulationData_integrated_chimera_controls.h5ad
data/processed/devguard/MouseGastrulationData_wt_chimera_full.h5ad
data/processed/devguard/MouseGastrulationData_tal1_chimera_full.h5ad
data/processed/devguard/MouseGastrulationData_t_chimera_full.h5ad
data/processed/devguard/MouseGastrulationData_embryo_atlas_sample1.h5ad
data/processed/devguard/MouseGastrulationData_wt_chimera_sample1.h5ad
data/processed/devguard/MouseGastrulationData_tal1_chimera_sample1.h5ad
data/processed/devguard/GSE123187_preview_4files.h5ad
```

这些文件均属于生成数据或下载数据，受 `.gitignore` 控制，不提交 Git。

## GSE212050 Strict Sample-Level Selection

用于当前真实 control calibration 主结果的入口是：

```text
scripts/devguard/select_gse212050_strict_sample_cells.py
config/devguard/normality_model_gse212050_strict_sample.json
```

筛选要求：

| 参数 | 数值 |
|---|---:|
| min cells per stage-lineage group | 200 |
| min sample units per stage-lineage group | 8 |
| max cells per stage-lineage group | 1,500 |
| allow cell fallback | false |

重跑结果中 9 / 9 个 reference groups 使用 `split_strategy=sample`，0 个 group 使用 cell fallback。

## 统一 obs schema

所有 processed H5AD 必须包含：

| 字段 | 含义 |
|---|---|
| `cell_id` | cell ID |
| `dataset_id` | dataset ID |
| `species` | 主线必须为 `Mus musculus` |
| `system` | embryo / gastruloid / organoid / embryo_chimera |
| `time_point` | E day / gastruloid day / hour |
| `time_numeric` | 数值化发育时间 |
| `condition` | control / perturbation |
| `perturbation_name` | KO gene / treatment |
| `perturbation_type` | genetic / chemical / pathway / environmental / none |
| `dose` | dose or NA |
| `duration` | treatment duration or NA |
| `sample_id` | embryo / organoid / donor ID |
| `batch` | batch |
| `cell_type` | cell type annotation |
| `lineage` | broad lineage |
| `is_control` | bool |
| `is_perturbed` | bool |

## 当前限制

- `GSE212050_strict_sample_13285` 是目前最可靠的真实 control calibration 入口，但仍是从 full Seurat object 中筛出的 balanced subset；正式论文应继续在完整数据或更大采样上复核。
- `MouseGastrulationData` full chimera 数据已经导出。当前主分析使用 integrated matched controls 和 E8.5 strict sample-level split；E7.5 WT chimera 只有 4 个 control sample units，不作为同等强度主校准。
- `GSE123187_preview_4files` 已进入 H5AD schema，但尚缺 cell type/lineage 注释和 axis mapping，不能直接作为 fate deviation validation。
