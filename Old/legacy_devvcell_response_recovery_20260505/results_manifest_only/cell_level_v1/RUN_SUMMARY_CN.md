# Cell-level baseline v1 运行摘要

运行日期：2026-05-02

## 数据子集

- 输入：`data/scLine_pro.h5ad`
- 输出：`data/processed/cell_level_subset_v1.h5ad`
- 细胞数：19,156
- 基因数：3,000
- 系统：neural、mesoderm_muscle、erythroid
- 阶段：Theiler stage 12--19
- 每个 system-stage 最多采样：800 个细胞

## 模型

- 表达编码：sparse TruncatedSVD
- latent 维度：64
- SVD 方差解释率总和：0.8231
- 伪靶标：同系统相邻阶段 Sinkhorn OT barycentric pseudo-targets
- 新增方法层：stage/system-conditioned context residual MLP
- Sinkhorn epsilon：0.12
- Sinkhorn iterations：120
- heldout target stages：15、18
- 训练伪配对：11,899
- heldout 伪配对：4,794

## 平均 heldout 指标

| model | pair latent MSE | centroid latent MSE | RBF-MMD |
|---|---:|---:|---:|
| context_residual_mlp | 0.3139 | 0.0178 | 0.0296 |
| ridge | 0.3254 | 0.0207 | 0.0445 |
| mlp | 0.3302 | 0.0177 | 0.0263 |
| identity | 0.3813 | 0.0230 | 0.0378 |
| mean_shift | 0.3856 | 0.0272 | 0.0444 |

## 解释边界

这些结果证明当前 cell-level 流水线可以从真实 H5AD 中导出细胞子集、构建相邻阶段 Sinkhorn OT 伪靶标并训练非平凡 transition operator。context residual MLP 将 source stage、target stage 和 DevVCell system 作为条件输入，并以残差形式学习发育推进，在当前 heldout pair MSE 上优于 ridge 和普通 MLP。当前研究包已经加入 zero-shot TF/GRN stimulus head、多随机种子 transition 消融和 GRN 组件消融；它们仍不是最终 Nature 水准结论，因为还缺少外部 perturbation benchmark、真实扰动标签、独立数据集泛化验证和更大规模统计检验。

## transition 统计比较

运行命令：

```powershell
python scripts\summarize_transition_statistics.py --metrics results\cell_level_v1\tables\cell_level_transition_metrics.csv --output-dir results\cell_level_v1 --reference-model context_residual_mlp
```

输出：

- `results/cell_level_v1/tables/transition_statistical_summary.csv`
- `results/cell_level_v1/tables/transition_paired_bootstrap_differences.csv`
- `results/cell_level_v1/tables/transition_statistical_summary.json`
- `results/cell_level_v1/figures/transition_bootstrap_ci.png`

配对 bootstrap 差异：

| reference | competitor | mean competitor-reference MSE | 95% CI | relative improvement |
|---|---|---:|---:|---:|
| context_residual_mlp | ridge | 0.0115 | 0.0014--0.0211 | 4.95% |
| context_residual_mlp | mlp | 0.0163 | 0.0058--0.0302 | 4.82% |
| context_residual_mlp | identity | 0.0674 | 0.0341--0.0966 | 23.99% |
| context_residual_mlp | mean_shift | 0.0717 | 0.0380--0.0996 | 25.31% |

## 可复用模型工件

- `results/cell_level_v1/models/state_svd.joblib`
- `results/cell_level_v1/models/state_scaler.joblib`
- `results/cell_level_v1/models/transition_ridge.joblib`
- `results/cell_level_v1/models/transition_mlp.pt`
- `results/cell_level_v1/models/transition_context_residual_mlp.pt`
- `results/cell_level_v1/models/mean_shift_vectors.npz`
- `results/cell_level_v1/models/model_config_snapshot.json`

这些文件固定了当前细胞状态空间和 transition baseline，可供下一步 TF/GRN stimulus head 复用。

## TF/GRN stimulus response head

运行命令：

```powershell
python scripts\run_stimulus_response_head.py --config config\cell_level_baseline.json
```

输出：

- `results/cell_level_v1/tables/cell_level_tf_grn_stimulus_response.csv`
- `results/cell_level_v1/tables/cell_level_tf_grn_stimulus_summary.csv`
- `results/cell_level_v1/figures/cell_level_tf_grn_stimulus_heatmap.png`
- `results/cell_level_v1/figures/cell_level_tf_recovery_scatter.png`
- `results/cell_level_v1/stimulus_response_summary.json`

当前 zero-shot GRN projection 覆盖：

- TF 数量：12
- DevVCell 系统：3
- TF--系统--阶段记录：252

按平均细胞级刺激响应强度排序的 top TF：

| TF | mean response norm | max response norm | mean recovery probability |
|---|---:|---:|---:|
| Eya1 | 0.5958 | 2.4844 | 0.2096 |
| Lef1 | 0.5799 | 2.1753 | 0.2024 |
| Rfx4 | 0.5144 | 1.9091 | 0.3045 |
| Nrg1 | 0.4968 | 1.8146 | 0.2357 |
| Satb2 | 0.4788 | 1.3016 | 0.2363 |

解释边界：该模块是基于 GRN 边的 zero-shot latent perturbation projection，没有使用外部 perturb-seq 或湿实验标签训练，因此当前用途是生成候选扰动假设和设计后续验证，而不是直接作为已验证扰动预测结论。

## 多 seed 和消融套件

运行命令：

```powershell
python scripts\run_ablation_suite.py --config config\ablation_suite.json
```

输出：

- `results/ablation_v1/tables/transition_ablation_metrics.csv`
- `results/ablation_v1/tables/transition_ablation_summary.csv`
- `results/ablation_v1/tables/stimulus_ablation_tf_summary.csv`
- `results/ablation_v1/tables/stimulus_ablation_summary.csv`
- `results/ablation_v1/figures/transition_ablation_summary.png`
- `results/ablation_v1/figures/stimulus_ablation_heatmap.png`

transition 消融摘要：

| pairing | model | mean heldout latent MSE | sd |
|---|---|---:|---:|
| nearest | ridge | 0.1413 | 0.0914 |
| nearest | context_residual_mlp | 0.1485 | 0.1009 |
| nearest | identity | 0.2218 | 0.1349 |
| nearest | mean_shift | 0.2232 | 0.1338 |
| sinkhorn | context_residual_mlp | 0.2974 | 0.2052 |
| sinkhorn | ridge | 0.3094 | 0.2041 |
| sinkhorn | identity | 0.3465 | 0.1940 |
| sinkhorn | mean_shift | 0.3528 | 0.1936 |

stimulus 消融摘要：

| variant | mean response | max response | mean recovery probability | n TF |
|---|---:|---:|---:|---:|
| no_system_edges | 0.5200 | 2.4820 | 0.2446 | 12 |
| no_global_grn | 0.5029 | 2.6597 | 0.2534 | 8 |
| full_grn | 0.4631 | 2.4109 | 0.2797 | 12 |

解释边界：当前消融评估的是模型内部稳定性和证据组件贡献，不能替代外部扰动数据验证。

## 统一复现清单

快速复现命令：

```powershell
python scripts\run_reproducible_pipeline.py --mode quick
```

该命令已重新聚合现有消融、两遍编译中文论文，并生成：

- `results/reproducibility_manifest.json`
- `manuscript/main_cn.pdf`

清单显示关键输出无缺失。
