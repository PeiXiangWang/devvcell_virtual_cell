# DevVCell response-recovery 正式运行报告

版本：2026-05-05

## 1. 执行结论

本轮已经完成从公开扰动数据下载、RDS 转 H5AD、外部 response dictionary、体内胚胎 manifold、response transfer、response-recovery classification、TS14-TS19 统计检验、minimal rescue control 和 niche-aware Tbx4-Glis3 case 的端到端正式入口运行。

统一入口已通过：

```bash
python scripts/run_response_recovery_pipeline.py --mode full
```

## 2. 环境和依赖

Python：

- `anndata 0.12.11`
- `pyarrow 21.0.0`
- `pandas/numpy/scipy/scikit-learn`

R：

- R 4.4.3
- 项目本地包库：`.r-lib/4.4`
- 已安装 `Seurat`、`SeuratObject`、`hdf5r`、`SeuratDisk`、`remotes`

注意：`SeuratDisk` 与当前 `SeuratObject 5.x` 在 `SaveH5Seurat()` 路径上存在 `slot=` 兼容性问题，因此正式转换改用自定义两步法：

```text
Seurat RDS -> Matrix Market/obs/var -> H5AD
```

对应脚本：

```text
scripts/export_seurat_rds_components.R
scripts/build_h5ad_from_seurat_components.py
```

## 3. 外部 perturbation 数据

主数据集：

```text
GSE208369
```

本地文件：

```text
data/external/GSE208369/GSE208369_KO_merged.RDS.gz
```

SHA256：

```text
DAFAD6C0567873FE449B876C010916409165A0E81AE164DE1CFE1A871E52907C
```

转换后 H5AD：

```text
data/external/response_transfer/primary_perturbation.h5ad
```

规模：

```text
20000 cells x 33577 genes
```

条件：

```text
NTC
KO_PAX3
KO_TBX6
RA_hGas_120h
```

## 4. 胚胎 developmental manifold

输入：

```text
data/scLine_pro.h5ad
```

输出：

```text
results/response_recovery/tables/embryo_manifold_cells.parquet
results/response_recovery/tables/stage_celltype_centroids.parquet
results/response_recovery/tables/manifold_quality_metrics.csv
```

本轮结果：

| 指标 | 值 |
|---|---:|
| sampled cells | 80,000 |
| stage/cell-type centroids | 632 |
| SVD explained variance ratio sum | 0.8277 |
| stage classifier accuracy | 0.2192 |
| cell type classifier accuracy | 0.7857 |

解释：cell type latent separation 较好；stage classification 仍弱，正式稿需要加入 time-aware/trajectory-aware encoder 或 donor/stage covariate control。

## 5. 外部 response dictionary 和 transfer

输出：

```text
results/perturbation_transfer/tables/external_response_dictionary.parquet
results/perturbation_transfer/tables/transferred_response_by_stage_celltype.parquet
```

规模：

| 输出 | 行数 |
|---|---:|
| external response dictionary | 33 |
| transferred responses | 20,856 |

外部响应来自 3 个 perturbation、11 个外部 cell type。转移到 632 个 embryo stage/cell-type centroids。

## 6. Response-recovery classification

输出：

```text
results/response_recovery/tables/response_recovery_classes.csv
results/response_recovery/tables/stage_vulnerability_response_recovery.csv
```

总分类数：

```text
20,856
```

类别分布：

| class | count |
|---|---:|
| reversible_response | 19,505 |
| fate_deflection | 622 |
| developmental_delay | 583 |
| off_manifold_collapse | 146 |

## 7. TS14-TS19 统计检验

输出：

```text
results/response_recovery/tables/window_enrichment_statistics.csv
```

关键结果：

| target | window mean | outside mean | p value | odds ratio |
|---|---:|---:|---:|---:|
| response_amplitude | 1.5108 | 1.2776 | 1.81e-92 | NA |
| recovery_cost | 1.5048 | 1.2749 | 1.85e-91 | NA |
| off_manifold_score | 1.4911 | 1.2665 | 1.78e-90 | NA |
| developmental_delay | 0.0338 | 0.0254 | 7.60e-04 | 1.34 |
| fate_deflection | 0.0487 | 0.0216 | 5.03e-26 | 2.31 |
| off_manifold_collapse | 0.0114 | 0.0051 | 4.62e-07 | 2.26 |

解释：TS14-TS19 在当前模型中确实呈现更高 response amplitude、recovery cost、fate deflection 和 off-manifold collapse 富集。但该结果仍需 no-OT、shuffled perturbation、donor/sex control 进一步确认。

## 8. Minimal rescue control

输出：

```text
results/response_recovery/tables/minimal_rescue_control_matrix.csv
```

规模：

```text
1,896 rows = 632 embryo centroids x 3 perturbations
```

方法：

对同一 stage/cell type 下的候选 perturbation response，求：

```text
min_beta || r_perturbation + beta * r_rescuer ||, 0 <= beta <= 1
```

输出字段包括：

- perturbation
- rescuer_perturbation
- optimal_rescuer_dose
- minimal_rescue_cost
- rescue_fraction
- response_recovery_class

当前结果中 rescue fraction 普遍偏低，说明 KO_PAX3、KO_TBX6 和 RA_hGas_120h 的响应方向在该 latent transfer 下并不形成强互相抵消。正式稿需要加入更丰富 perturbation dictionary 和 gene/module-level rescuer。

## 9. Tbx4-Glis3 niche-aware case

输出：

```text
results/niche_context/tables/niche_signature_by_stage_donor.csv
results/niche_context/tables/tbx4_glis3_niche_case.csv
results/niche_context/tables/cell_autonomous_vs_niche_scores.csv
```

本轮 sampled cells：

```text
120,000
```

结果：

| case | cell-autonomous Spearman | niche-mediated Spearman | interpretation |
|---|---:|---:|---|
| Tbx4-Glis3 | 0.1727 | NA | cell_autonomous_or_unresolved |

解释：Tbx4 和 Glis3 均在当前 3,000 gene embryo matrix 中被找到，但 donor/stage niche-mediated correlation 当前不足以支持强结论。正式稿只能写成 unresolved/candidate，不应写直接调控或强 niche-mediated。

## 10. 验证

已通过：

```bash
python scripts/validate_external_perturbation.py
python -m pytest
python -m compileall scripts src tests
```

测试：

```text
3 passed
```

验证报告：

```text
results/response_recovery/tables/external_validation_report.csv
```

## 11. 下一步必须做的正式消融

1. no-OT transfer：用 nearest/mean response 替代 transfer。
2. shuffled perturbation：打乱 perturbation labels 后重跑 classification。
3. no external response：与 quick RDEG proxy 对照。
4. no niche：去掉 niche features 后比较 Tbx4-Glis3 case。
5. donor/sex covariate control：检查 TS14-TS19 是否仍显著。
6. module/gene-level shared encoder：修正当前 external latent 与 embryo latent 不是同一训练空间的限制。
