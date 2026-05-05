# DevVCell Response-Recovery 方法说明

版本：2026-05-05

## 1. 模型对象

DevVCell 的基本对象是胚胎发育流形上的 stage/cell-type 状态：

```text
M = {z(stage, cell_type, donor, sex, niche)}
```

扰动响应不是单独预测一个表达矩阵，而是预测扰动状态与正常发育流形之间的关系：

```text
z_perturbed = z_normal + transferred_response + feedback_term
```

## 2. 发育流形构建

输入：

```text
data/scLine_pro.h5ad
```

脚本：

```text
scripts/build_developmental_manifold.py
```

最小正式实现使用 sparse matrix + TruncatedSVD latent。正式投稿版可以替换为 scVI/autoencoder，但必须保留同一套输出 schema：

- `embryo_manifold_cells`
- `stage_celltype_centroids`
- `manifold_quality_metrics`

质量指标：

- stage classifier accuracy
- cell type classifier accuracy
- SVD/scVI explained variance or reconstruction loss
- donor/sex mixing diagnostics
- heldout stage reconstruction error

## 3. 外部扰动响应字典

输入：

```text
data/external/response_transfer/primary_perturbation.h5ad
```

脚本：

```text
scripts/export_external_perturbation_response.py
```

输出：

- `external_response_dictionary`
- `external_response_metadata`

每个 perturbation response 至少包括：

- perturbation label
- control condition
- external cell type
- control cell count
- perturbed cell count
- response vector
- response norm
- encoding quality

跨物种数据必须先做 human-mouse ortholog mapping，或者在 module-level response 上工作。

## 4. 响应转移

脚本：

```text
scripts/transfer_perturbation_response_ot.py
```

正式模型应实现：

1. external control states 与 embryo centroid states 的 cost matrix；
2. entropic OT / Sinkhorn transfer plan；
3. barycentric response transfer；
4. transfer confidence；
5. no-OT nearest/mean-response 消融。

当前代码提供 soft confidence transfer 接口，正式数据接入后要把 response dictionary 与 embryo manifold 放到共享 latent 或 module space。

## 5. Response-recovery classification

脚本：

```text
scripts/classify_response_recovery.py
```

核心距离：

```text
response_amplitude = ||z_perturbed - z_current||
developmental_delay_score = d(z_perturbed, early_same_fate) - d(z_perturbed, z_current)
fate_deflection_index = d(z_perturbed, future_same_fate) - d(z_perturbed, future_other_fate)
off_manifold_score = min_z_in_M d(z_perturbed, z)
recovery_cost = min(d(z_perturbed, z_current), d(z_perturbed, future_same_fate))
```

分类规则：

1. `off_manifold_score` 超过正常 centroid nearest-neighbor 分布阈值：`off_manifold_collapse`。
2. `developmental_delay_score < -delay_margin`：`developmental_delay`。
3. `fate_deflection_index > fate_margin`：`fate_deflection`。
4. 否则：`reversible_response`。

阈值写在：

```text
config/response_recovery.json
```

## 6. TS14-TS19 脆弱性统计

正式稿不能只展示 TS14-TS19 图形富集。必须输出统计比较：

- TS14-TS19 vs outside 的 recovery cost 分布检验；
- delay/deflection/off-manifold class proportion enrichment；
- donor/sex/batch covariate 控制；
- shuffled perturbation negative control；
- downsampling sensitivity。

当前统计入口：

```bash
python scripts/analyze_response_recovery_statistics.py
```

输出：

```text
results/response_recovery/tables/window_enrichment_statistics.csv
```

## 7. Minimal rescue control

当前正式实现从 transferred response vectors 出发，对同一 stage/cell type 下的候选 perturbation response 计算最优抵消剂量：

```text
min_beta || r_perturbation + beta * r_rescuer ||, 0 <= beta <= 1
```

输出每个 perturbation 的 minimal residual cost、rescue fraction 和最佳 rescuer perturbation：

```bash
python scripts/compute_minimal_rescue_control.py
```

```text
results/response_recovery/tables/minimal_rescue_control_matrix.csv
```

## 8. 消融实验

必须包含：

| 消融 | 实现 |
|---|---|
| no external response | 只运行 `--mode quick` 或使用 RDEG proxy response |
| no OT transfer | 用 nearest centroid 或 mean response |
| no feedback | 去掉 rescue/recovery term |
| no niche | 不加入 `build_niche_context.py` 输出 |
| no sex/donor covariate | centroid 不分 donor/sex 或在统计模型中去掉 covariate |
| shuffled perturbation | 打乱 perturbation labels 后重跑 classification |

当前基础消融入口：

```bash
python scripts/run_response_recovery_ablations.py
```

输出：

```text
results/response_recovery/ablations/ablation_summary.csv
```

已实现的基础消融：

- `shuffled_response_vectors`
- `no_stage_celltype_transfer`
- `global_mean_response`
- `no_external_rdeg_proxy`

## 9. 当前可运行入口

快速检查：

```bash
python scripts/run_response_recovery_pipeline.py --mode quick
python scripts/validate_external_perturbation.py
```

正式运行：

```bash
python scripts/build_developmental_manifold.py
python scripts/export_external_perturbation_response.py
python scripts/transfer_perturbation_response_ot.py
python scripts/classify_response_recovery.py
python scripts/analyze_response_recovery_statistics.py
python scripts/compute_minimal_rescue_control.py
python scripts/build_niche_context.py
python scripts/validate_external_perturbation.py
```
