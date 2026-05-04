# DevVCell 复现说明

## 环境

当前本机环境已验证：

- Python 3.11
- pandas、numpy、matplotlib
- scipy、scikit-learn
- anndata、scanpy、h5py
- torch

建议从项目根目录运行：

```powershell
python -m pip install -r requirements.txt
```

## 一键复现入口

当前项目新增统一复现脚本：

```powershell
python scripts\run_reproducible_pipeline.py --mode quick
```

`quick` 模式会检查和聚合已有结果、重新编译中文论文，并写出 `results/reproducibility_manifest.json`。如果要重新运行主分析而不重复完整消融：

完成审计记录在：

```text
docs/COMPLETION_AUDIT_CN.md
```

同行方法对比和方法定位记录在：

```text
docs/PEER_METHOD_COMPARISON_CN.md
```

外部扰动 benchmark 入口记录在：

```text
docs/EXTERNAL_PERTURBATION_BENCHMARK_CN.md
```

```powershell
python scripts\run_reproducible_pipeline.py --mode main
```

如果要从现有配置重新运行多 seed 消融：

```powershell
python scripts\run_reproducible_pipeline.py --mode full
```

可选参数：

- `--force-subset`：强制重新导出 `data/processed/cell_level_subset_v1.h5ad`。
- `--skip-ablation`：跳过消融套件。
- `--skip-external-benchmark`：跳过公开 scPerturb H5AD 下载、校验和外部 benchmark。
- `--external-config`：指定外部扰动 benchmark 配置。
- `--no-compile-paper`：不运行 `xelatex`。

## 原型分析

```powershell
python scripts\devvcell_lite.py
python scripts\ot_grn_developmental_impact.py
```

输出：

- `results/tables/stage_vulnerability.csv`
- `results/tables/perturbation_priority.csv`
- `results/tables/virtual_cell_rollout_proxy.csv`
- `results/tables/ot_grn_tf_system_developmental_impact.csv`
- `results/figures/*.png`

## 细胞级基线

先从完整 H5AD 导出小型平衡子集：

```powershell
python scripts\export_cell_level_subset.py --config config\cell_level_baseline.json --force
```

再训练和评估相邻阶段 transition baseline：

```powershell
python scripts\train_cell_transition_baseline.py --config config\cell_level_baseline.json
```

当前主配置会同时训练 `context_residual_mlp`，该模型使用 source stage、target stage 和 DevVCell system one-hot 作为 context，并通过内部 validation split 和 early stopping 选择最佳 epoch。

输出：

- `data/processed/cell_level_subset_v1.h5ad`
- `data/processed/cell_level_subset_v1.manifest.csv`
- `results/cell_level_v1/tables/cell_level_transition_metrics.csv`
- `results/cell_level_v1/tables/cell_level_pair_manifest.csv`
- `results/cell_level_v1/figures/cell_level_transition_baseline_mse.png`
- `results/cell_level_v1/training_summary.json`
- `results/cell_level_v1/models/state_svd.joblib`
- `results/cell_level_v1/models/state_scaler.joblib`
- `results/cell_level_v1/models/transition_ridge.joblib`
- `results/cell_level_v1/models/transition_mlp.pt`
- `results/cell_level_v1/models/transition_context_residual_mlp.pt`
- `results/cell_level_v1/models/mean_shift_vectors.npz`

当前默认配置使用 Sinkhorn OT barycentric pseudo-targets；可在 `config/cell_level_baseline.json` 的 `model.pairing.method` 中改为 `nearest` 以复现最近邻伪配对基线。

运行 transition 统计比较：

```powershell
python scripts\summarize_transition_statistics.py --metrics results\cell_level_v1\tables\cell_level_transition_metrics.csv --output-dir results\cell_level_v1 --reference-model context_residual_mlp
```

输出：

- `results/cell_level_v1/tables/transition_statistical_summary.csv`
- `results/cell_level_v1/tables/transition_paired_bootstrap_differences.csv`
- `results/cell_level_v1/tables/transition_statistical_summary.json`
- `results/cell_level_v1/figures/transition_bootstrap_ci.png`

运行 TF/GRN stimulus response head：

```powershell
python scripts\run_stimulus_response_head.py --config config\cell_level_baseline.json
```

输出：

- `results/cell_level_v1/tables/cell_level_tf_grn_stimulus_response.csv`
- `results/cell_level_v1/tables/cell_level_tf_grn_stimulus_summary.csv`
- `results/cell_level_v1/figures/cell_level_tf_grn_stimulus_heatmap.png`
- `results/cell_level_v1/figures/cell_level_tf_recovery_scatter.png`
- `results/cell_level_v1/stimulus_response_summary.json`

## 多 seed 和消融

运行消融套件：

```powershell
python scripts\run_ablation_suite.py --config config\ablation_suite.json
```

若 transition runs 已经存在，只重跑 stimulus 和重画图：

```powershell
python scripts\run_ablation_suite.py --config config\ablation_suite.json --skip-transition
```

若所有 run 已经存在，只重新聚合和画图：

```powershell
python scripts\run_ablation_suite.py --config config\ablation_suite.json --skip-transition --skip-stimulus
```

输出：

- `results/ablation_v1/tables/transition_ablation_metrics.csv`
- `results/ablation_v1/tables/transition_ablation_summary.csv`
- `results/ablation_v1/tables/stimulus_ablation_tf_summary.csv`
- `results/ablation_v1/tables/stimulus_ablation_summary.csv`
- `results/ablation_v1/figures/transition_ablation_summary.png`
- `results/ablation_v1/figures/stimulus_ablation_heatmap.png`
- `results/ablation_v1/ablation_summary.json`

## 外部扰动数据入口

下载并检查公开 scPerturb H5AD：

```powershell
python scripts\prepare_external_perturbation_benchmark.py --record-id 10044268 --file-name DatlingerBock2021.h5ad --output-dir data\external\scperturb
python scripts\run_external_perturbation_benchmark.py --config config\external_perturbation_benchmark.json
```

输出：

- `data/external/scperturb/DatlingerBock2021.h5ad`
- `data/external/scperturb/scperturb_benchmark_manifest.json`
- `results/external_perturbation_v1/tables/external_perturbation_metrics.csv`
- `results/external_perturbation_v1/tables/external_perturbation_condition_metrics.csv`
- `results/external_perturbation_v1/figures/external_perturbation_benchmark_mse.png`
- `results/external_perturbation_v1/external_perturbation_summary.json`

该步骤当前实现 guide `_1` 训练、guide `_2` 测试的真实扰动标签 benchmark。完整外部验证仍需加入 CellOT、scGen/CPA 或 Systema 风格 baseline。

## 投稿合规要点

Nature 官方说明强调初投稿可将正文与图合并到一个 Word 或 PDF，Methods 需包含足以复现实验的要素，图件需足够清晰供审稿人评估。最终稿阶段，LaTeX 需要同时提供 PDF，图件通常按单栏或双栏宽度准备，并提供数据与代码可用性声明。

本项目的论文稿必须包含：

- Data availability
- Code availability
- Methods 中的统计与复现细节
- LLM use statement
- Competing interests
- Author contributions
