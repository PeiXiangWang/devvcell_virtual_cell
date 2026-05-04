# DevVCell-lite

`DevVCell-lite` 是与原 `embryonic_development` 项目平行的新项目，用于把小鼠胚胎发育单细胞图谱从“阶段级 RDEG 重建”推进到“发育虚拟细胞”的原型研究。

新项目路径：

```text
C:\Users\14915\PycharmProjects\devvcell_virtual_cell
```

原项目路径：

```text
C:\Users\14915\PycharmProjects\embryonic_development
```

## 数据迁移状态

完整原始数据已经迁移到新项目：

```text
data/scLine_pro.h5ad
```

源文件：

```text
C:\Users\14915\PycharmProjects\embryonic_development\scLine_pro.h5ad
```

迁移后目标文件：

```text
C:\Users\14915\PycharmProjects\devvcell_virtual_cell\data\scLine_pro.h5ad
```

文件大小校验：

| 项目 | 值 |
|---|---:|
| 源文件字节数 | 12,620,552,772 |
| 目标文件字节数 | 12,620,552,772 |
| H5AD 矩阵形状 | 11,441,407 cells x 3,000 genes |
| 非零表达值数量 | 1,422,703,167 |
| 表达矩阵密度 | 4.14% |
| 表达矩阵稀疏度 | 95.86% |

除了完整 H5AD，本项目还保留了从 RDEG neural MVP 复制过来的轻量派生结果表：

```text
data/rdeg_neural_cell_mvp/
```

这些派生表用于快速运行 `DevVCell-lite` 原型，不替代原始表达矩阵。

## 数据集概览

`scLine_pro.h5ad` 是一个 AnnData/H5AD 格式的单细胞转录组数据集，主体是多阶段小鼠胚胎发育图谱。它包含从 Theiler stage 12 到 Theiler stage 27 的细胞表达矩阵、细胞注释、基因注释和二维 UMAP 坐标。

核心规模：

| 内容 | 数量 |
|---|---:|
| 细胞数 | 11,441,407 |
| 基因数 | 3,000 |
| Theiler stages | 16 |
| author_day 时间点 | 43 |
| donor embryos | 74 |
| 标准 cell_type 类别 | 134 |
| 作者原始 author_cell_type 类别 | 190 |
| 作者 major cluster 类别 | 26 |
| 性别类别 | 2 |
| disease 类别 | 1，全部为 normal |

## H5AD 内部结构

数据文件的主要结构如下：

```text
X
  CSR sparse matrix, shape = 11,441,407 x 3,000
  X/data dtype = float32
  X/indices dtype = int32
  X/indptr length = 11,441,408

obs
  cell-level metadata, 11,441,407 rows

var
  gene-level metadata, 3,000 rows

obsm/X_umap
  precomputed UMAP embedding, shape = 11,441,407 x 2
```

表达矩阵 `X` 是 CSR 稀疏矩阵。行是细胞，列是基因，值是表达量。由于稀疏度约 95.86%，后续模型训练应尽量使用稀疏读取、分块读取或按系统/阶段导出子集，不建议一次性转成稠密矩阵。

## obs 字段说明

`obs` 是细胞级元数据表，包含以下字段：

| 字段 | 类型 | 说明 |
|---|---|---|
| `donor_id` | categorical | 胚胎或样本供体编号，共 74 类 |
| `author_experimental_id` | categorical | 作者原始实验批次或实验 ID |
| `author_day` | categorical | 作者给出的发育时间点，共 43 类，例如 E0925、E1000、P0000 |
| `author_somite_count` | categorical | 作者记录的体节数或相关阶段信息 |
| `author_major_cell_cluster` | categorical | 作者原始大类细胞群，共 26 类 |
| `author_cell_type` | categorical | 作者原始细胞类型注释，共 190 类 |
| `cell_type` | categorical | 标准化后的细胞类型注释，共 134 类 |
| `disease` | categorical | 疾病状态；本数据中全部为 normal |
| `sex` | categorical | 性别，female 或 male |
| `development_stage` | categorical | Theiler stage，共 16 类，覆盖 stage 12 到 stage 27 |

## var 字段说明

`var` 是基因级元数据表，包含以下字段：

| 字段 | 类型 | 说明 |
|---|---|---|
| `gene_short_name` | categorical | 基因短名，例如 Eya1、Lef1、Glis3 |
| `chr` | categorical | 染色体位置 |
| `start` | numeric | 基因起始坐标 |
| `end` | numeric | 基因结束坐标 |
| `strand` | categorical | 链方向，`+` 或 `-` |
| `gene_type` | categorical | Ensembl/注释体系中的基因类型 |
| `feature_is_filtered` | numeric/bool | feature 是否被过滤的标记 |
| `feature_name` | categorical | feature 名称 |
| `feature_reference` | categorical | feature reference |
| `feature_biotype` | categorical | feature 生物类型；本数据中为 gene |
| `feature_length` | categorical | feature 长度信息 |
| `feature_type` | categorical | feature 类型 |

前 30 个基因短名示例：

```text
Gm37381, Rp1, Rgs20, St18, Sntg1, 2610203C22Rik, Mybl1, Sgk3,
Cpa6, Prex2, A830018L16Rik, Sulf1, Eya1, Defb41, Tfap2d, Tfap2b,
Pkhd1, Gm28653, Gsta3, Kcnq5, Col9a1, Col19a1, Khdrbs2, Arhgef4,
Gm37068, Zap70, Npas2, Tbc1d8, Il1r2, Il1r1
```

## Theiler stage 分布

`development_stage` 覆盖 16 个 Theiler stages：

| Theiler stage | 细胞数 |
|---|---:|
| Theiler stage 27 | 2,092,528 |
| Theiler stage 22 | 1,016,140 |
| Theiler stage 15 | 970,745 |
| Theiler stage 20 | 889,367 |
| Theiler stage 21 | 860,115 |
| Theiler stage 23 | 827,534 |
| Theiler stage 24 | 800,306 |
| Theiler stage 26 | 766,708 |
| Theiler stage 16 | 759,720 |
| Theiler stage 25 | 537,192 |
| Theiler stage 17 | 515,342 |
| Theiler stage 18 | 414,430 |
| Theiler stage 19 | 396,957 |
| Theiler stage 14 | 394,473 |
| Theiler stage 12 | 153,597 |
| Theiler stage 13 | 46,253 |

对 DevVCell 来说，Theiler stage 是正常发育 transition、heldout stage 验证、阶段脆弱性图谱和刺激窗口分析的核心时间轴。

## author_day 时间点分布

数据包含 43 个 `author_day` 时间点。细胞数最多的前 20 个时间点如下：

| author_day | 细胞数 |
|---|---:|
| P0000 | 2,092,528 |
| E0975 | 633,746 |
| E1000 | 489,516 |
| E0950 | 336,999 |
| E1400 | 326,063 |
| E1225 | 321,031 |
| E1325 | 305,678 |
| E1850 | 289,835 |
| E0925 | 276,424 |
| E1425 | 274,509 |
| E1433.3 | 273,184 |
| E1050 | 271,152 |
| E1025 | 270,204 |
| E1300 | 269,386 |
| E1075 | 244,190 |
| E1575 | 239,821 |
| E1150 | 233,843 |
| E1650 | 228,892 |
| E1675 | 222,748 |
| E1125 | 222,125 |

`author_day` 比 Theiler stage 更细，但不一定等间隔。后续可以把它作为辅助时间变量或 donor/time covariate。

## 主要 cell_type 分布

标准化 `cell_type` 共 134 类。前 20 类如下：

| cell_type | 细胞数 |
|---|---:|
| neural cell | 1,500,763 |
| mesodermal cell | 1,463,106 |
| erythroid progenitor cell | 1,033,409 |
| lateral mesodermal cell | 745,494 |
| glutamatergic neuron | 725,228 |
| fibroblast | 648,059 |
| neuron | 630,208 |
| cerebral cortex GABAergic interneuron | 603,448 |
| GABAergic neuron | 376,415 |
| hepatocyte | 289,569 |
| chondrocyte | 274,786 |
| myoblast | 185,751 |
| keratinocyte | 175,178 |
| osteoblast | 135,893 |
| neural progenitor cell | 123,796 |
| epithelial cell | 111,975 |
| retinal progenitor cell | 110,254 |
| astrocyte | 109,750 |
| primitive erythroid progenitor | 108,613 |
| endo-epithelial cell | 105,092 |

这些细胞类型适合分成几个 DevVCell 训练子任务：

1. neural/neuron/glia 系统：用于 Wnt、FGF、Lef1、neural differentiation 相关响应。
2. mesoderm/lateral mesoderm 系统：用于胚层转移、组织分支和阶段敏感性分析。
3. erythroid/hematopoietic 系统：用于高细胞数、强分化轨迹的 benchmark。
4. reproductive/PGC-like/Leydig-like 相关系统：用于连接原 RDEG 的 PGC-Leydig 生殖系统案例。

## 作者 major cluster 分布

`author_major_cell_cluster` 共 26 类。前 20 类如下：

| major cluster | 细胞数 |
|---|---:|
| Mesoderm | 3,267,338 |
| CNS_neurons | 2,106,206 |
| Neuroectoderm_and_glia | 1,733,663 |
| Definitive_erythroid | 1,033,409 |
| Intermediate_neuronal_progenitors | 628,251 |
| Epithelial_cells | 524,960 |
| Endothelium | 312,029 |
| Muscle_cells | 305,003 |
| Hepatocytes | 289,569 |
| White_blood_cells | 262,022 |
| Neural_crest_PNS_glia | 126,743 |
| Adipocytes | 114,478 |
| Primitive_erythroid | 108,613 |
| Neural_crest_PNS_neurons | 103,999 |
| Eye_and_other | 93,695 |
| T_cells | 89,429 |
| Lung_and_airway | 76,732 |
| Intestine | 51,796 |
| B_cells | 44,352 |
| Olfactory_sensory_neurons | 42,826 |

这个字段比 `cell_type` 更粗，可以用于分层抽样、系统级训练、跨系统验证和论文图中的大类展示。

## 作者原始 author_cell_type 分布

`author_cell_type` 共 190 类。前 20 类如下：

| author_cell_type | 细胞数 |
|---|---:|
| Definitive early erythroblasts (CD36-) | 940,963 |
| Lateral plate and intermediate mesoderm | 745,494 |
| Glutamatergic neurons | 725,228 |
| Fibroblasts | 648,059 |
| GABAergic cortical interneurons | 603,448 |
| Sclerotome | 476,671 |
| Telencephalon | 430,409 |
| Spinal cord/r7/r8 | 395,636 |
| Facial mesenchyme | 382,752 |
| GABAergic neurons | 376,415 |
| Hepatocytes | 289,569 |
| Upper-layer neurons | 238,829 |
| Deep-layer neurons | 236,402 |
| Limb mesenchyme progenitors | 223,422 |
| Dermatome | 198,888 |
| Hindbrain | 197,367 |
| Myoblasts | 185,751 |
| Midbrain | 145,124 |
| Early chondrocytes | 144,311 |
| Pre-osteoblasts (Sp7+) | 135,893 |

该字段更接近原始注释，可以在需要更细粒度命运标签时使用。

## donor 与 sex 分布

数据包含 74 个 donor。细胞数最多的前 10 个 donor：

| donor_id | 细胞数 |
|---|---:|
| embryo_74 | 1,613,834 |
| embryo_73 | 478,694 |
| embryo_53 | 326,063 |
| embryo_41 | 321,031 |
| embryo_71 | 289,835 |
| embryo_54 | 274,509 |
| embryo_55 | 273,184 |
| embryo_34 | 271,152 |
| embryo_33 | 270,204 |
| embryo_35 | 244,190 |

性别分布：

| sex | 细胞数 |
|---|---:|
| female | 6,134,796 |
| male | 5,306,611 |

疾病状态：

| disease | 细胞数 |
|---|---:|
| normal | 11,441,407 |

后续建模时，`donor_id` 和 `sex` 应作为 covariate 或 batch/context 信息，避免把 donor composition 当作发育信号。

## 基因注释概览

`gene_type` 前 20 类：

| gene_type | 基因数 |
|---|---:|
| protein_coding | 2,541 |
| lincRNA | 167 |
| antisense | 100 |
| processed_transcript | 47 |
| TEC | 41 |
| processed_pseudogene | 37 |
| unprocessed_pseudogene | 23 |
| transcribed_unprocessed_pseudogene | 10 |
| miRNA | 6 |
| TR_C_gene | 5 |
| sense_intronic | 5 |
| IG_C_gene | 3 |
| transcribed_processed_pseudogene | 3 |
| snoRNA | 2 |
| snRNA | 2 |
| TR_J_gene | 2 |
| unitary_pseudogene | 2 |
| scaRNA | 1 |
| pseudogene | 1 |
| polymorphic_pseudogene | 1 |

`feature_type` 主要分布：

| feature_type | 基因数 |
|---|---:|
| protein_coding | 2,537 |
| lncRNA | 326 |
| TEC | 41 |
| processed_pseudogene | 37 |
| transcribed_unprocessed_pseudogene | 21 |
| unprocessed_pseudogene | 8 |
| miRNA | 6 |

染色体分布前 10：

| chr | 基因数 |
|---|---:|
| chr1 | 250 |
| chr7 | 238 |
| chr6 | 229 |
| chr2 | 220 |
| chr11 | 206 |
| chr3 | 204 |
| chr5 | 201 |
| chr9 | 151 |
| chr10 | 147 |
| chr4 | 142 |

链方向：

| strand | 基因数 |
|---|---:|
| `-` | 1,588 |
| `+` | 1,412 |

## DevVCell 研究中的使用方式

完整 H5AD 在新项目中的用途如下：

1. 训练正常发育虚拟细胞：从相邻 Theiler stages 中学习 cell-level transition。
2. 构造 OT pseudo target：由于单细胞测序是破坏性的，不能直接追踪同一细胞，因此用相邻阶段分布匹配形成软监督。
3. 训练 stimulus-conditioned transition：把 TF knockdown、TF overexpression、pathway activation 或环境刺激作为外部输入。
4. 评估 feedback recovery：比较扰动后状态与正常发育流形的距离、恢复成本和替代命运概率。
5. 制作发育脆弱性图谱：按 stage、cell type、major cluster 和 perturbation 汇总 response amplitude 与 recovery cost。

## 推荐读取方式

不要把全矩阵一次性转成 dense array。推荐使用 backed 模式或按阶段/细胞类型分块导出。

示例：

```python
import scanpy as sc

adata = sc.read_h5ad("data/scLine_pro.h5ad", backed="r")
print(adata.shape)
print(adata.obs.columns.tolist())
print(adata.var.columns.tolist())
```

抽取一个阶段：

```python
stage = "Theiler stage 15"
mask = adata.obs["development_stage"] == stage
stage15 = adata[mask, :].to_memory()
```

抽取一个系统：

```python
mask = adata.obs["cell_type"].str.contains("neural|neuron|glia", case=False, regex=True)
neural = adata[mask, :].to_memory()
```

如果只做快速统计，优先读取 `obs` 和 `var`，避免触碰完整 `X`。

## 当前项目结构

```text
devvcell_virtual_cell/
  config/
    devvcell_lite.json
    cell_level_baseline.json
    ablation_suite.json
  data/
    scLine_pro.h5ad
    processed/
    rdeg_neural_cell_mvp/
  docs/
    REPRODUCIBILITY_CN.md
    NATURE_DELIVERY_ROADMAP_CN.md
    devvcell_work_plan_cn.tex
    devvcell_work_plan_cn.pdf
    CODE_PLAN.md
  manuscript/
    main_cn.tex
    main_cn.pdf
    figures/
  scripts/
    devvcell_lite.py
    ot_grn_developmental_impact.py
    export_cell_level_subset.py
    train_cell_transition_baseline.py
    run_stimulus_response_head.py
    run_ablation_suite.py
    run_reproducible_pipeline.py
  results/
    tables/
    figures/
    cell_level_v1/
    ablation_v1/
```

## DevVCell-lite 原型运行

当前原型脚本仍然使用轻量派生表，不直接读取 12GB H5AD：

```powershell
cd C:\Users\14915\PycharmProjects\devvcell_virtual_cell
python scripts\devvcell_lite.py
```

输出：

```text
results/tables/stage_vulnerability.csv
results/tables/perturbation_priority.csv
results/tables/virtual_cell_rollout_proxy.csv
results/tables/execution_summary.json
results/figures/stage_vulnerability.png
results/figures/perturbation_priority_top12.png
results/figures/recovery_cost_vs_response.png
results/figures/module_step_distance.png
```

## 统一复现入口

当前研究包可用统一入口检查已有结果并重新编译中文论文：

```powershell
python scripts\run_reproducible_pipeline.py --mode quick
```

重新运行主分析：

```powershell
python scripts\run_reproducible_pipeline.py --mode main
```

重新运行主分析和多 seed 消融：

```powershell
python scripts\run_reproducible_pipeline.py --mode full
```
