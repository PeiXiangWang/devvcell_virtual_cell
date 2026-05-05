# DevVCell niche-aware case study：Tbx4-Glis3

版本：2026-05-05

## 1. 核心问题

Tbx4-Glis3 不应被过早写成直接 TF-target。正式版问题应改为：

> Tbx4-Glis3 动态关联更像 cell-autonomous signal，还是 niche-mediated developmental response？

这样可以把第一篇 RDEG 的候选边升级为第二篇 DevVCell 的 niche-aware virtual cell 示例。

## 2. 生物学解释框架

可能机制：

1. Tbx4-high mesoderm/gonadal niche 在特定 stage 出现；
2. Glis3-high germ-cell/PGC-like state 在同一 donor/stage 增强；
3. 二者未必在同一细胞共表达；
4. 如果 donor/stage co-occurrence 或空间邻近强于 cell-autonomous 共表达，则写成 candidate niche-mediated response。

## 3. 当前实现

脚本：

```text
scripts/build_niche_context.py
```

配置：

```text
config/niche_context.json
```

输出：

```text
results/niche_context/tables/niche_signature_by_stage_donor.csv
results/niche_context/tables/tbx4_glis3_niche_case.csv
results/niche_context/tables/cell_autonomous_vs_niche_scores.csv
```

## 4. 当前指标

### 4.1 Niche signature

按 stage/donor 计算：

- cell type composition；
- candidate niche cell fraction；
- candidate target cell fraction；
- Tbx4-high niche fraction；
- Glis3-high target fraction；
- Tbx4-Glis3 co-occurrence score。

### 4.2 Cell-autonomous score

在 candidate target cells 内计算：

```text
Spearman(Tbx4 expression, Glis3 expression)
```

### 4.3 Niche-mediated score

在 donor/stage 层面计算：

```text
Spearman(Tbx4-high niche fraction, Glis3-high target fraction)
```

若 niche-mediated score 高于 cell-autonomous score，初步写成：

```text
candidate niche-mediated developmental response
```

## 5. 正式投稿所需增强

当前 donor/stage/cell-type composition 是 proxy。正式稿应至少补强一个层级：

1. 加入 ligand-receptor module；
2. 加入 spatial/tomo-seq validation；
3. 加入 fetal gonad/PGC external validation；
4. 对 donor、sex、stage 做 covariate control；
5. 做 shuffled donor/stage negative control。

## 6. 推荐外部验证

### GSE123187

mouse gastruloid scRNA-seq + tomo/spatial transcriptomics，可用于检查 niche proxy 是否有空间或轴向支持。

### GSE136441

prenatal/neonatal mouse gonad/ovary scRNA-seq，可用于 PGC/gonadal niche validation。

### GSE288206

murine fetal gonad single-nucleus multiomics，可用于 PGC/supporting-lineage regulatory validation。

## 7. 写作边界

可以写：

> DevVCell nominates Tbx4-Glis3 as a candidate niche-mediated developmental response linking mesodermal/gonadal context with Glis3-high germ-cell states.

不能写，除非有空间或实验验证：

> Tbx4 directly regulates Glis3 in PGCs.

## 8. 下一步

1. 在 `scLine_pro.h5ad` 上运行 `build_niche_context.py`。
2. 查看 `cell_autonomous_vs_niche_scores.csv`。
3. 如果 gene coverage 不足，转为 module-level Tbx4 niche signature 和 Glis3 target signature。
4. 下载并转换 GSE123187 或 GSE136441 做外部验证。
