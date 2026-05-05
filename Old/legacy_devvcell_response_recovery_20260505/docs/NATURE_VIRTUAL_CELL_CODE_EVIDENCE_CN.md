# DevVCell 正式论文级 evidence pipeline 说明

本项目使用正式配置驱动的 evidence pipeline 来论证“发育虚拟细胞”主张。入口为：

```powershell
python scripts\run_nature_virtual_cell_evidence_pipeline.py --config config\nature_virtual_cell_evidence.json
```

该管线不是轻量汇总脚本。它执行输入存在性检查、schema 校验、必需模型检查、统计检验、bootstrap 置信区间、provenance 记录、claim gate、manifest 和论文 Methods 输出。

## 配置

```text
config/nature_virtual_cell_evidence.json
```

配置项包括：

- 输入表路径。
- primary model 和 baseline model。
- competence window。
- minimum rows 和 claim thresholds。
- nullable numeric columns：只允许在配置中显式声明的“非适用”数值列为空，且仍必须至少包含一个有限值。
- failure policy。
- bootstrap reps 和 seed。

## 失败条件

默认硬失败条件：

- 缺失必需输入文件。
- 输入表缺少必需列。
- 必需数值列缺失、无法解析、含无穷值，或未在 `nullable_numeric_columns` 中声明却出现空值。
- 表行数低于配置阈值。
- 缺失 primary 或 baseline model。

默认非硬失败、但会进入 `claim_gate_matrix.csv` 的条件：

- transition paired strata 不足。
- external Perturb-seq condition 不足。
- 相对改进未达到阈值。
- competence window 行数不足。

如需把 claim gate 失败也作为硬失败：

```powershell
python scripts\run_nature_virtual_cell_evidence_pipeline.py --config config\nature_virtual_cell_evidence.json --fail-on-claim-threshold-failure
```

正式契约测试：

```powershell
python -m unittest discover -s tests -v
```

## 输出

```text
results/nature_virtual_cell_evidence/
```

关键文件：

- `evidence_manifest.json`：完整 provenance、输入 hash、输出 hash、claim gate 和运行环境。
- `paper_methods_evidence_pipeline.md`：可直接写入论文 Methods 的管线说明。
- `tables/input_file_inventory.csv`：输入文件清单和 SHA-256。
- `tables/input_schema_validation.csv`：schema 校验结果。
- `tables/required_model_validation.csv`：primary/baseline 模型检查。
- `tables/transition_model_performance.csv`：transition 主模型性能。
- `tables/transition_paired_statistical_tests.csv`：transition 配对统计检验。
- `tables/competence_window_group_statistics.csv`：competence window 分组统计。
- `tables/competence_window_statistical_tests.csv`：窗口内外统计检验。
- `tables/vulnerability_response_correlation_tests.csv`：vulnerability-response 相关检验。
- `tables/fate_recovery_virtual_screen.csv`：TF-system-stage fate/recovery 虚拟筛选。
- `tables/external_perturbation_model_performance.csv`：外部 Perturb-seq 模型性能。
- `tables/external_perturbation_paired_statistical_tests.csv`：外部扰动配对统计检验。
- `tables/external_perturbation_bias_control_summary.csv`：condition-centered bias-control 摘要。
- `tables/transition_ablation_formal_summary.csv`：transition 消融摘要。
- `tables/stimulus_ablation_formal_summary.csv`：stimulus 消融摘要。
- `tables/claim_gate_matrix.csv`：论文主张 gate 和当前缺口。

关键图：

- `figures/figure_transition_formal_benchmark.png`
- `figures/figure_competence_window_response.png`
- `figures/figure_fate_recovery_virtual_screen.png`
- `figures/figure_external_perturbation_paired_tests.png`

## 主张与证据

| 主张 | 代码证据 |
|---|---|
| DevVCell 是细胞级虚拟细胞，不是 RDEG 下游图谱包装 | `tables/virtual_cell_component_contract.csv` |
| stage/system-conditioned residual dynamics 改善 heldout transition | `tables/transition_model_performance.csv`、`tables/transition_paired_statistical_tests.csv` |
| competence window 调制扰动响应 | `tables/competence_window_statistical_tests.csv`、`figures/figure_competence_window_response.png` |
| TF 扰动输出命运偏移和恢复成本 | `tables/fate_recovery_virtual_screen.csv`、`figures/figure_fate_recovery_virtual_screen.png` |
| 外部 Perturb-seq 校准优于 naive baseline | `tables/external_perturbation_paired_statistical_tests.csv` |
| 当前证据边界可审计 | `tables/claim_gate_matrix.csv`、`evidence_manifest.json` |

## 解释边界

当前 TF/GRN stimulus、fate displacement 和 recovery cost 仍是计算 proxy。正式 Nature 级投稿还需要：

- CellOT、GEARS、scGen/CPA、foundation model embedding 的直接基线。
- 至少 3 个外部 Perturb-seq 数据集。
- heldout donor、heldout system、heldout gene、heldout context 和 heldout cell state。
- 阶段依赖 CRISPRi/a 或 Perturb-seq 湿实验验证。
