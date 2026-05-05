# DevGuard 新版第二篇文章与代码实施方案

版本：2026-05-05

建议题目：**DevGuard: stage-specific developmental normality boundaries for perturbed mouse embryonic cells**

中文题目：**DevGuard：小鼠胚胎扰动细胞的阶段特异发育正常性判定框架**

## 核心问题

DevGuard 只回答一个问题：在小鼠胚胎或小鼠 gastruloid 发育系统中，扰动后的细胞是否仍处于对应发育阶段和谱系的正常范围内。如果不在正常范围内，它是发育延迟、发育提前、命运偏航，还是进入正常发育图谱外的异常状态。

## 与旧 DevVCell/RDEG 路线切割

新版主线不使用：

- RDEG pipeline；
- RDEG-derived node graph；
- OT transition 或 response transfer；
- GRN dynamics；
- matrix writeback；
- `data/scLine_pro.h5ad` 作为主证据链。

旧版 DevVCell response-recovery 原型已归档到：

```text
Old/legacy_devvcell_response_recovery_20260505/
```

## 方法路线

DevGuard 从 control mouse embryo/gastruloid time-course 中按 `time_point + lineage` 建立正常参考集：

```text
N(t, l) = control cells at time t and lineage l
```

对每个参考集分出 train、calibration 和 heldout test。训练集用于学习 embedding 和 reference distribution；calibration 集用于 nonconformity score 的 conformal 校准；heldout test 用于检查 control false-positive rate。

MVP 实现两种 nonconformity score：

- KNN distance；
- Mahalanobis distance。

Conformal p-value：

```text
p(x; t,l) = (1 + #{s_cal >= s(x; t,l)}) / (1 + n_cal)
```

当 `p(x; t,l) >= alpha` 时，认为细胞仍在该 stage-lineage 正常范围内。

## 五类扰动状态

分类优先级固定为：

1. `within_stage_normal`
2. `developmental_delay`
3. `developmental_acceleration`
4. `fate_deviation`
5. `abnormal_off_normal`

判定逻辑：

- 当前时间、同谱系通过：当前阶段正常；
- 当前不通过，但早期同谱系通过：发育延迟；
- 当前不通过，但晚期同谱系通过：发育提前；
- 当前不通过，但其他谱系通过：命运偏航；
- 所有正常参考集均不通过：正常图谱外异常。

## 第一版交付边界

当前代码完成 conformal normality MVP 和 quick-mode 软件验证。真实公开数据下载、MouseGastrulationData 导入、GSE212050 organoid barcode 解析、GSE123187 spatial/tomo 解析将作为后续数据接入工作推进。

Quick mode 生成小型 mouse gastruloid fixture，仅用于验证 schema、normality reference、FPR、五类分类、DTI 和 Figure 1-3 出图流程，不作为生物学结论。
