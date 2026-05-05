# DevVCell 正式改动计划：发育响应-恢复虚拟细胞

版本：2026-05-05

## 1. 核心决策

第二篇文章不再定位为“把 RDEG 输出包装成 virtual cell”，而是正式转向：

> **A developmental virtual cell maps response-recovery landscapes across embryogenesis.**

中文主线：

> **发育虚拟细胞揭示胚胎扰动后的响应-恢复景观。**

文章要回答的问题是：胚胎细胞在特定发育阶段受到扰动后，能否回到正常发育流形，是否只是发育延迟，是否偏航到替代命运，还是进入正常图谱外的异常状态。

## 2. 和第一篇 RDEG 的边界

第一篇 RDEG 解决正常发育如何推进：宏观 OT、微观 GRN、rollout、矩阵回写和候选调控边。

第二篇 DevVCell 解决扰动后能否恢复：外部真实扰动响应、体内胚胎流形、response-recovery classification、minimal rescue cost 和 niche-mediated case。

因此，第一篇结果只作为 normal developmental manifold 和候选调控先验；第二篇主结果必须放在扰动后恢复、延迟、偏航和 off-manifold 风险上。

## 3. 正式版创新点

### 3.1 Response-recovery classification

每个 stage/cell type/perturbation 的虚拟扰动状态分为四类：

| 类别 | 含义 | 投稿叙事 |
|---|---|---|
| `reversible_response` | 短暂偏离后可回到正常流形 | 胚胎恢复力 |
| `developmental_delay` | 更接近早期正常状态 | 发育迟滞 |
| `fate_deflection` | 更接近替代命运分支 | 命运偏航 |
| `off_manifold_collapse` | 离所有正常状态过远 | 图谱外异常 |

核心指标：

- `response_amplitude`：扰动状态离同阶段正常状态的距离。
- `recovery_cost`：回到原命运正常流形的最小代价。
- `developmental_delay_score`：扰动状态相对当前状态是否更像早期状态。
- `fate_deflection_index`：扰动状态相对原命运是否更接近替代命运。
- `off_manifold_score`：扰动状态到所有正常状态的最近距离。

### 3.2 Perturbation-response transfer

正式版不只使用 RDEG-derived proxy。外部扰动数据必须进入模型：

1. 从公开 perturbation scRNA-seq 中学习 control-to-perturbed response。
2. 通过共享特征或 ortholog mapping 把外部状态编码到可对齐空间。
3. 用 OT 或 soft barycentric alignment 把响应转移到体内胚胎 stage/cell type。
4. 在胚胎流形上判定恢复、延迟、偏航和 off-manifold。

主张句：

> Existing methods learn either normal developmental couplings or perturbation maps. DevVCell transfers real perturbation responses into an in vivo embryonic developmental manifold and predicts whether cells recover, delay, or deflect.

### 3.3 Niche-aware virtual cell

虚拟细胞定义为：

```text
VirtualCell_i = (z_i, n_i)
```

其中 `z_i` 是细胞自身表达状态，`n_i` 是 donor/stage/cell-type composition、major cluster composition、ligand-receptor module 或空间邻近构成的 niche context。

Tbx4-Glis3 不再硬写成直接 TF-target，而是作为 candidate niche-mediated developmental response：

> RDEG proposed a Tbx4-Glis3 dynamic association; DevVCell resolves it as a candidate niche-mediated developmental response.

## 4. 已落地仓库改动

新增配置：

- `config/response_recovery.json`
- `config/perturbation_transfer.json`
- `config/niche_context.json`
- `config/external_datasets.json`

新增脚本：

- `scripts/build_developmental_manifold.py`
- `scripts/export_external_perturbation_response.py`
- `scripts/transfer_perturbation_response_ot.py`
- `scripts/classify_response_recovery.py`
- `scripts/analyze_response_recovery_statistics.py`
- `scripts/compute_minimal_rescue_control.py`
- `scripts/build_niche_context.py`
- `scripts/run_response_recovery_pipeline.py`
- `scripts/run_response_recovery_ablations.py`
- `scripts/validate_external_perturbation.py`
- `scripts/download_public_datasets.py`
- `scripts/convert_seurat_rds_to_h5ad.R`

新增文档：

- `docs/DEVVCELL_RESPONSE_RECOVERY_CHANGE_PLAN_CN.md`
- `docs/RESPONSE_RECOVERY_METHOD_CN.md`
- `docs/EXTERNAL_DATASET_PLAN_CN.md`
- `docs/NICHE_AWARE_CASE_STUDY_CN.md`

## 5. 正式实验路线

### Result 1：体内胚胎发育流形

输入 `data/scLine_pro.h5ad`，输出 sampled cell latent、stage/cell-type centroids、stage/cell type classifier accuracy、heldout stage baseline。

对应脚本：

```bash
python scripts/build_developmental_manifold.py
```

### Result 2：外部扰动响应字典

优先数据集是 GSE208369，因为它包含 conventional/RA gastruloids、LDN-treated、TBX6-KO 和 PAX3-KO RA-gastruloid samples，能支撑 RA/BMP/TF response transfer。

对应脚本：

```bash
python scripts/export_external_perturbation_response.py --allow-missing
```

`--allow-missing` 只用于记录数据获取计划；正式运行必须放入转换后的 H5AD。

### Result 3：OT response transfer

把外部 response dictionary 转移到 embryo stage/cell-type centroid：

```bash
python scripts/transfer_perturbation_response_ot.py
```

正式版要求补齐 human-mouse ortholog mapping 或统一 gene/module encoder，避免跨物种 latent basis 不一致。

### Result 4：response-recovery landscape

输出四类 response-recovery class，并统计 TS14-TS19 是否富集 delay/deflection/off-manifold：

```bash
python scripts/classify_response_recovery.py
```

### Result 5：minimal rescue control

基于 response amplitude、recovery cost、RDEG rescue fraction 和外部 perturbation dictionary，形成 perturbation-rescuer matrix。当前保留在 quick proxy 中，正式版要接入真实 external response。

### Result 6：Tbx4-Glis3 niche-aware case

计算 donor/stage cell-type composition、Tbx4-high niche 与 Glis3-high target co-occurrence，并区分 cell-autonomous 与 niche-mediated signal：

```bash
python scripts/build_niche_context.py
```

## 6. 数据集决策

优先级：

1. GSE208369：RA/BMP/TF perturbation response dictionary。
2. GSE212050：mouse gastruloid time course，用于 organoid heterogeneity 和 false-positive control。
3. GSE123187：mouse gastruloid scRNA-seq + tomo/spatial validation。
4. MouseGastrulationData / E-MTAB-6967：独立 early embryo manifold validation。
5. MOCA：organogenesis stage validation。
6. GSE136441 或 GSE288206：PGC/gonad niche validation。

原始大文件不提交 Git；只提交 registry、转换脚本、校验信息和分析结果摘要。

## 7. 不可降低的投稿标准

正式版必须满足：

- 至少一个真实外部 perturbation scRNA-seq 数据集进入 response dictionary。
- 必须有 no external response、no OT transfer、no niche、shuffled perturbation 消融。
- TS14-TS19 高风险窗口必须有统计检验，而不是只看图。
- Tbx4-Glis3 只能写成 candidate niche-mediated response，除非有空间或实验验证。
- README 和 manuscript 中必须明确 DevVCell-lite 只是 feasibility proxy，不是正式主模型。

## 8. 当前执行入口

快速回归：

```bash
python scripts/run_response_recovery_pipeline.py --mode quick
```

正式主线：

```bash
python scripts/run_response_recovery_pipeline.py --mode main
```

全量分析：

```bash
python scripts/run_response_recovery_pipeline.py --mode full
```

正式运行前需要把 `config/external_datasets.json` 中选定的数据下载并转换为 `data/external/response_transfer/primary_perturbation.h5ad`。
