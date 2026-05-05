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
| `GSE212050` | organoid heterogeneity negative control | Seurat RDS 已下载，15k downsample H5AD 已构建，control reference 已跑通 |
| `GSE123187` | spatial/tomo validation | RAW tar 已下载，4-file preview H5AD 已构建 |
| `devguard_quick_mouse` | software fixture | quick mode 可生成 |
| `devguard_stress_mouse` | software stress fixture | stress mode 可生成 |

## 已验证的本地入口

```text
data/processed/devguard/GSE212050_downsample_15000.h5ad
data/processed/devguard/MouseGastrulationData_embryo_atlas_sample1.h5ad
data/processed/devguard/MouseGastrulationData_wt_chimera_sample1.h5ad
data/processed/devguard/MouseGastrulationData_tal1_chimera_sample1.h5ad
data/processed/devguard/GSE123187_preview_4files.h5ad
```

这些文件均属于生成数据或下载数据，受 `.gitignore` 控制，不提交 Git。

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

- `GSE212050_downsample_15000` 是可运行的真实 control reference，但一部分 stage-lineage group 的 sample-unit/heldout cell 数较低，正式结果需要完整数据或更大 downsample。
- `MouseGastrulationData` chimera smoke test 目前只导出 sample 1，不能形成真正 embryo-level split，reference 会退化为 cell-level fallback。
- `GSE123187_preview_4files` 已进入 H5AD schema，但尚缺 cell type/lineage 注释和 axis mapping，不能直接作为 fate deviation validation。
