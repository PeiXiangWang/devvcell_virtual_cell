# DevGuard 方法说明

## Normality Reference

DevGuard 只用 control cells 建立正常参考。每个 `time_point + lineage` 组合是一个 stage-lineage reference group。组内 control cells 随机划分为：

- training set；
- calibration set；
- heldout test set。

训练集用于计算 reference embedding 和分布参数；calibration set 用于 conformal p-value；heldout test set 用于估计 control false-positive rate。

## Embedding

MVP 使用 counts normalization、`log1p`、高变基因筛选和 `TruncatedSVD`。该选择可审计、可复现，并避免引入深度模型依赖。后续可增加 scVI 作为 sensitivity analysis。

## Nonconformity Scores

当前实现：

- `knn_distance`：query cell 到 reference training cells 的 KNN 平均距离；
- `mahalanobis`：query cell 到 reference training distribution 的正则化 Mahalanobis distance。

## Conformal Calibration

对 calibration cells 计算 score。query cell 的 conformal p-value 为：

```text
(1 + number of calibration scores >= query score) / (1 + number of calibration cells)
```

若当前 stage-lineage p-value 大于等于 `alpha`，则判定为 `within_stage_normal`。

## Perturbed-Cell Classification

每个 perturbation cell 会对所有 reference group 计算 p-value，再汇总：

- `p_current_same`
- `p_early_same`
- `p_late_same`
- `p_other_lineage`
- `p_any_normal`

分类优先级：

```text
within_stage_normal
> developmental_delay
> developmental_acceleration
> fate_deviation
> abnormal_off_normal
```

## Developmental Tolerance Index

对每个 time point / lineage / perturbation：

```text
DTI = R_normal - R_delay - R_accel - R_deviation - R_abnormal
```

DTI 越高表示扰动越被该发育状态耐受；DTI 越低表示该状态更脆弱。
