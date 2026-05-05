# 外部扰动基准入口

## 目的

DevVCell 当前的 TF/GRN stimulus head 主要是机制投影和假设生成层。为了把方法有效性推进到可投稿证据链，项目新增一个公开 CRISPR perturbation benchmark：先从 scPerturb 下载和检查 H5AD，再在真实扰动标签上运行 guide-transfer 评估。

## 已接入数据

脚本：

```powershell
python scripts\prepare_external_perturbation_benchmark.py --record-id 10044268 --file-name DatlingerBock2021.h5ad --output-dir data/external/scperturb
```

已生成：

- `data/external/scperturb/DatlingerBock2021.h5ad`
- `data/external/scperturb/scperturb_benchmark_manifest.json`

当前清单记录的关键事实：

- 数据来源：scPerturb Zenodo record `10044268`。
- 文件：`DatlingerBock2021.h5ad`。
- 文件大小：33,600,780 bytes。
- AnnData 规模：39,194 cells x 25,904 genes。
- 主要候选扰动列：`perturbation` 和 `perturbation_2`。
- 控制标签：`perturbation == control`，共 4,497 个控制细胞。
- 评估上下文：`perturbation_2` 中的 `stimulated` 和 `unstimulated`。
- 文献：Peidli et al., Nature Methods 2024, doi:10.1038/s41592-023-02144-y。

## 已实现评估

运行：

```powershell
python scripts\run_external_perturbation_benchmark.py --config config\external_perturbation_benchmark.json
```

设计：

- 数据集：Datlinger/Bock 2021 Jurkat T cell CRISPR perturbation 数据。
- 切分：同一基因的 guide `_1` 训练，guide `_2` 作为 heldout guide 测试。
- 条件：20 个扰动基因、2 个 TCR stimulation 上下文。
- 表示：3,000 个表达特征基因，64 维 TruncatedSVD latent。
- 伪靶标：控制细胞到同一 context 下扰动细胞的 Sinkhorn barycentric pseudo-target。
- 规模：8,144 个训练伪配对，8,391 个 heldout guide 伪配对，39 个 heldout 条件。

当前结果：

| 模型 | heldout guide MSE | centroid MSE | effect cosine | RBF-MMD |
|---|---:|---:|---:|---:|
| gene/context ridge residual | 0.2209 | 0.0101 | 0.5339 | 0.0699 |
| gene/context neural residual | 0.2818 | 0.0101 | 0.4997 | 0.0341 |
| global mean shift | 0.7886 | 0.0118 | 0.1775 | 0.0816 |
| identity | 0.7889 | 0.0121 | NA | 0.0803 |
| guide-transfer mean shift | 0.7975 | 0.0207 | 0.1893 | 0.0999 |

结果文件：

- `results/external_perturbation_v1/tables/external_perturbation_metrics.csv`
- `results/external_perturbation_v1/tables/external_perturbation_condition_metrics.csv`
- `results/external_perturbation_v1/tables/external_perturbation_pair_manifest.csv`
- `results/external_perturbation_v1/figures/external_perturbation_benchmark_mse.png`
- `results/external_perturbation_v1/external_perturbation_summary.json`

## 当前边界

该 benchmark 已经提供真实扰动标签上的首个外部效果数字，但仍不是完整 Nature 投稿级外部验证。下一步必须加入 CellOT、scGen/CPA 或 Systema 风格 baseline，增加 heldout gene 和 heldout context 切分，并在原始基因空间报告 DE gene rank recovery 与校准指标。
