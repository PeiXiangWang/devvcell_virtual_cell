# DevVCell-lite 代码实施计划

本文档记录本项目的敲代码计划。当前目标不是一次性完成大型虚拟细胞模型，而是把平行 RDEG 项目的数据资产转成一个可复现、可扩展、可写进论文计划的最小原型。

## 0. 项目边界

新项目路径：

```text
C:\Users\14915\PycharmProjects\devvcell_virtual_cell
```

复用源项目：

```text
C:\Users\14915\PycharmProjects\embryonic_development
```

当前不复制原始 `scLine_pro.h5ad`，只复制 RDEG 派生结果表。后续若要训练真正 cell-level 模型，再单独导出小型 AnnData 子集或以只读路径引用原始 H5AD。

## 1. 已完成的初始代码任务

- [x] 新建平行项目目录。
- [x] 建立 `config/`、`data/`、`docs/`、`scripts/`、`results/` 目录。
- [x] 复制 RDEG neural MVP 的轻量结果数据。
- [x] 编写 `config/devvcell_lite.json`。
- [x] 编写 `scripts/devvcell_lite.py`。
- [x] 编写 `README.md` 和数据 manifest。
- [x] 生成本文档作为代码计划。

## 2. DevVCell-lite 原型脚本计划

脚本：

```text
scripts/devvcell_lite.py
```

### 2.1 数据读取

读取以下输入：

```text
data/rdeg_neural_cell_mvp/stage_module_means.csv
data/rdeg_neural_cell_mvp/temporal_sensitivity.csv
data/rdeg_neural_cell_mvp/tf_knockout_results.csv
data/rdeg_neural_cell_mvp/rescue_experiments.csv
data/rdeg_neural_cell_mvp/rollout_error_matrix.csv
data/rdeg_neural_cell_mvp/ot_pair_metrics.csv
data/rdeg_neural_cell_mvp/next_step_pair_metrics.csv
data/rdeg_neural_cell_mvp/nodes_summary.csv
```

### 2.2 阶段脆弱性指标

生成 `stage_vulnerability.csv`，字段包括：

- `vulnerability_score`
- `mean_sensitivity`
- `mean_rollout_error`
- `baseline_dev_step_norm`
- `mean_outlier_ratio`
- `next_step_mse`
- `ot_reg_gain`

解释：

`vulnerability_score` 是 prototype 指标，用于把“阶段敏感性、rollout 不稳定性、正常发育模块步长、节点 outlier 和下一步预测误差”合成一个阶段脆弱性分数。

### 2.3 扰动优先级指标

生成 `perturbation_priority.csv`，字段包括：

- `devvcell_priority_score`
- `response_amplitude_proxy`
- `fate_displacement_proxy`
- `feedback_cost`
- `recovery_probability_proxy`
- `peak_stage`
- `best_rescuer_tf`

解释：

这些是 DevVCell 论文语言的初始指标，不等同于真实湿实验结论。它们用于筛出后续值得训练、benchmark 或实验验证的 TF perturbation。

### 2.4 虚拟细胞 rollout proxy

生成 `virtual_cell_rollout_proxy.csv`。当前不做深度生成，只用阶段模块状态和 top TF 扰动分数构造一个“虚拟细胞响应强度”的代理表。后续正式版本会替换为：

```text
cell encoder -> OT pseudo target -> transition operator -> stimulus head -> GRN feedback -> decoder
```

### 2.5 图件输出

输出：

```text
results/figures/stage_vulnerability.png
results/figures/perturbation_priority_top12.png
results/figures/recovery_cost_vs_response.png
results/figures/module_step_distance.png
```

## 3. 运行命令

```powershell
cd C:\Users\14915\PycharmProjects\devvcell_virtual_cell
python scripts\devvcell_lite.py
```

## 4. 下一阶段代码计划

- [ ] 从原始 H5AD 导出 3 个系统的小型 AnnData 子集：neural、mesodermal、reproductive。
- [ ] 实现 cell-level feature loader：HVG + TF + marker/pathway genes。
- [ ] 实现轻量 autoencoder 或 PCA/scVI adapter 接口。
- [ ] 实现相邻 stage 的细胞级 mini-batch Sinkhorn coupling。
- [ ] 实现 barycentric pseudo target。
- [ ] 实现 transition operator 的训练和 1 到 5 步 rollout。
- [ ] 实现 stimulus embedding：TF knockdown、TF overexpression、pathway activation/inhibition。
- [ ] 接入 GRN feedback head。
- [ ] 实现 ablation：no OT、no GRN、no feedback、no stimulus。
- [ ] 加入外部 perturbation benchmark 数据接口。
- [ ] 将结果自动写入 LaTeX 表格和图件目录。

## 5. 验收标准

当前 MVP 验收：

- [x] 新项目能独立运行，不依赖修改原 RDEG 代码。
- [x] 轻量数据已复制到新项目。
- [x] 原型脚本能生成表和图。
- [x] 文档清楚标注哪些是 proxy，哪些是下一步真实模型任务。

正式模型验收：

- [ ] 可从真实细胞矩阵训练 encoder/transition。
- [ ] 可在 heldout stage 上生成表达矩阵。
- [ ] 可接受 TF/pathway/environment perturbation。
- [ ] 可输出 response amplitude、feedback cost、fate displacement 和 recovery probability。
- [ ] 至少一个外部 perturbation 数据集或小规模实验能验证模型预测。
