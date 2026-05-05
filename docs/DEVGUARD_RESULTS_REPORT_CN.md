# DevGuard 结果报告

## 当前状态

本仓库已完成 DevGuard quick-mode MVP：

- 旧 DevVCell response-recovery 原型归档；
- 新建 `src/devguard/` 和 `scripts/devguard/`；
- 建立 dataset registry 和统一 obs schema；
- 实现 SVD embedding、KNN/Mahalanobis nonconformity、conformal p-value；
- 实现 heldout control FPR calibration；
- 实现 perturbed cell 五类分类；
- 实现 Developmental Tolerance Index；
- 实现 Figure 1-3 草图生成。

## Quick Mode 输出

运行：

```bash
python scripts/devguard/run_devguard_pipeline.py --mode quick
```

会生成：

```text
data/processed/devguard/devguard_quick_mouse.h5ad
results/devguard/normality_reference/
results/devguard/perturbation_classification/
results/devguard/tolerance_index/
results/devguard/figures/
```

## Quick Mode 详细验证结果

当前 quick fixture 是 deterministic mouse gastruloid software fixture，仅用于验证代码路径、schema、conformal calibration 和五类分类逻辑。

规模：

| 项目 | 数值 |
|---|---:|
| control cells | 1,350 |
| perturbed cells | 120 |
| genes | 90 |
| control time points | 3 |
| lineages | 3 |
| reference groups | 9 |
| score-method records | 18 |

Heldout control calibration：

| 指标 | 数值 |
|---|---:|
| alpha | 0.05 |
| max heldout control FPR | 0.0333 |
| score methods | KNN distance, Mahalanobis |

Perturbed-cell classification：

| quick fixture pattern | DevGuard class | cell count |
|---|---|---:|
| `within_stage_normal_like` | `within_stage_normal` | 24 |
| `delay_like` | `developmental_delay` | 24 |
| `acceleration_like` | `developmental_acceleration` | 24 |
| `fate_deviation_like` | `fate_deviation` | 24 |
| `abnormal_like` | `abnormal_off_normal` | 24 |

Condition-level summary：

| class | fraction |
|---|---:|
| `within_stage_normal` | 0.20 |
| `developmental_delay` | 0.20 |
| `developmental_acceleration` | 0.20 |
| `fate_deviation` | 0.20 |
| `abnormal_off_normal` | 0.20 |

DTI：

```text
DTI(E7.5, mesoderm, FGF8_KO)
= 0.20 - 0.20 - 0.20 - 0.20 - 0.20
= -0.60
```

测试：

```text
python -m pytest
9 passed
```

## 解释边界

Quick-mode fixture 只用于软件验证，不用于论文中的生物学结论。真实结果必须在 E-MTAB-6967、GSE212050、GSE123187 或其他独立公开小鼠扰动数据接入后重新生成。
