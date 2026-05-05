# DevVCell 投稿前研究包完成审计

审计日期：2026-05-03

## 目标拆解

目标：将 DevVCell-lite 推进为可复现的发育虚拟细胞研究包，包括细胞级模型训练评估流水线、高质量中文框架图、中文 LaTeX 论文稿、实验路线、结果表图和复现代码，定位为 Nature 水准投稿前研究包。

## Prompt-to-artifact 清单

| 要求 | 证据文件或命令 | 当前状态 |
|---|---|---|
| 项目说明和数据基础 | `README.md`；`data/scLine_pro.h5ad`；`data/rdeg_neural_cell_mvp/` | 已覆盖。README 记录 H5AD 规模、obs/var 字段、Theiler stage 分布和推荐读取方式。 |
| 细胞级子集导出 | `scripts/export_cell_level_subset.py`；`config/cell_level_baseline.json`；`data/processed/cell_level_subset_v1.h5ad`；`data/processed/cell_level_subset_v1.manifest.csv` | 已覆盖。当前子集为 19,156 cells x 3,000 genes。 |
| 细胞级模型训练评估流水线 | `scripts/train_cell_transition_baseline.py`；`results/cell_level_v1/training_summary.json`；`results/cell_level_v1/models/`；`results/cell_level_v1/tables/cell_level_transition_metrics.csv` | 已覆盖。包含 SVD latent、Sinkhorn OT pseudo-target、ridge、MLP、identity 和 mean-shift baseline。 |
| 条件化残差发育转移算子 | `config/cell_level_baseline.json`；`scripts/train_cell_transition_baseline.py`；`results/cell_level_v1/models/transition_context_residual_mlp.pt` | 已覆盖。`context_residual_mlp` 使用 source/target stage 与 system context，并通过 validation split 和 early stopping 训练。 |
| transition 统计比较 | `scripts/summarize_transition_statistics.py`；`results/cell_level_v1/tables/transition_statistical_summary.csv`；`results/cell_level_v1/tables/transition_paired_bootstrap_differences.csv`；`results/cell_level_v1/figures/transition_bootstrap_ci.png` | 已覆盖。对 6 个 heldout system-stage strata 做 bootstrap CI 和配对差异。 |
| TF/GRN stimulus 响应层 | `scripts/run_stimulus_response_head.py`；`results/cell_level_v1/stimulus_response_summary.json`；`results/cell_level_v1/tables/cell_level_tf_grn_stimulus_response.csv` | 已覆盖。当前输出 252 条 TF--系统--阶段响应记录，12 个 TF 和 3 个系统。 |
| 外部扰动基准 | `scripts/prepare_external_perturbation_benchmark.py`；`scripts/run_external_perturbation_benchmark.py`；`config/external_perturbation_benchmark.json`；`results/external_perturbation_v1/` | 已覆盖首版外部 benchmark。Datlinger/Bock 2021 H5AD 为 39,194 cells x 25,904 genes；guide `_1` 训练、guide `_2` heldout 测试中，`gene_context_ridge_residual` MSE 为 0.2209，优于 identity/global shift 的约 0.789。 |
| 多 seed 和组件消融 | `scripts/run_ablation_suite.py`；`config/ablation_suite.json`；`results/ablation_v1/tables/transition_ablation_metrics.csv`；`results/ablation_v1/tables/stimulus_ablation_tf_summary.csv` | 已覆盖。transition 消融 144 行；stimulus TF 消融 32 行。 |
| 结果表图 | `results/figures/`；`results/cell_level_v1/figures/`；`results/ablation_v1/figures/`；`results/external_perturbation_v1/figures/` | 已覆盖。包含阶段脆弱性、扰动优先级、OT+GRN heatmap、cell-level transition、stimulus heatmap、recovery scatter、消融图和外部扰动 benchmark 图。 |
| 中文深度模型生成框架图 | `manuscript/figures/devvcell_framework_cn_ai.png`；`manuscript/figures/FIGURE_PROVENANCE_CN.md`；`manuscript/figure_prompts/framework_figure_cn.md` | 已覆盖。来源记录声明由 Codex `image_gen` 深度图像生成工具生成，并保留 prompt。 |
| 中文 LaTeX 论文稿 | `manuscript/main_cn.tex`；`manuscript/main_cn.pdf` | 已覆盖。PDF 包含中文正文、10 张图、方法、数据/代码可用性、LLM 使用声明和参考文献。 |
| 实验路线和 Nature 投稿前路线 | `docs/NATURE_DELIVERY_ROADMAP_CN.md`；`docs/REPRODUCIBILITY_CN.md`；`results/cell_level_v1/RUN_SUMMARY_CN.md` | 已覆盖。路线图明确当前证据、扩展路线、风险和必须补强的外部验证。 |
| 同行方法对比 | `docs/PEER_METHOD_COMPARISON_CN.md`；`manuscript/main_cn.tex` | 已覆盖。对 Waddington-OT、moscot、CellRank、CellOT、scGen、scVI、scGPT 和 Systema 进行定位比较。 |
| 复现代码和一键入口 | `scripts/run_reproducible_pipeline.py`；`requirements.txt`；`results/reproducibility_manifest.json` | 已覆盖。`quick` 模式已重新聚合消融并两遍编译论文，关键输出无缺失。 |

## 实际验证

已运行的关键验证：

```powershell
python -m py_compile scripts\run_reproducible_pipeline.py scripts\train_cell_transition_baseline.py scripts\run_stimulus_response_head.py scripts\run_ablation_suite.py
python -m py_compile scripts\prepare_external_perturbation_benchmark.py
python -m py_compile scripts\run_external_perturbation_benchmark.py
python scripts\prepare_external_perturbation_benchmark.py --record-id 10044268 --file-name DatlingerBock2021.h5ad --output-dir data\external\scperturb
python scripts\run_external_perturbation_benchmark.py --config config\external_perturbation_benchmark.json
python scripts\summarize_transition_statistics.py --metrics results\cell_level_v1\tables\cell_level_transition_metrics.csv --output-dir results\cell_level_v1 --reference-model context_residual_mlp
python scripts\run_reproducible_pipeline.py --mode quick
rg -n "undefined|Undefined|LaTeX Warning|Error|Fatal|Missing|Citation" manuscript\main_cn.log
pdftotext manuscript\main_cn.pdf manuscript\main_cn_text.txt
```

验证结果：

- Python 脚本语法检查通过。
- 外部 scPerturb benchmark 已生成 manifest 和模型评估表，记录 39,194 个细胞、25,904 个基因、20 个扰动基因、2 个 stimulation context 和 39 个 heldout guide 条件。
- `summarize_transition_statistics.py` 通过，并显示 `context_residual_mlp` 为当前主 benchmark pair MSE 最低模型。
- `run_reproducible_pipeline.py --mode quick` 通过。
- `results/reproducibility_manifest.json` 中 `missing_critical_outputs` 为空。
- `manuscript/main_cn.log` 未匹配未解析引用、LaTeX warning、fatal error 或 missing 图件。
- PDF 文本抽取显示参考文献、数据可用性、代码可用性和 LLM 使用声明均存在。

## 审计结论

当前状态满足“Nature 水准投稿前研究包”的工程交付定义：已有可复现代码、细胞级模型训练评估结果、中文论文、中文 AI 框架图、结果表图、实验路线和复现清单。

仍需在正式投稿前由真实研究团队补强的事项：

- 外部 perturbation benchmark 或独立数据集验证。
- 更正式的统计检验、置信区间和 effect size 报告。
- 真实作者、单位、作者贡献、数据访问号和代码 release。
- 由作者或专业绘图人员对 AI 框架图进行最终标签和排版校验。
- 真实同行评议和期刊接收不能由本地代码或自动化流程保证。
