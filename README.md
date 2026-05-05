# DevGuard

DevGuard is an independent mouse developmental perturbation-normality
project. It does not use the RDEG pipeline, RDEG-derived node graphs, OT
transitions, GRN dynamics, or matrix-writeback outputs. Instead, DevGuard
learns stage- and lineage-specific normality boundaries from public mouse
embryo and gastruloid control time-course datasets, then classifies
perturbed cells as within-stage normal, developmentally delayed,
accelerated, fate-deviated, or abnormal/off-normal.

The previous DevVCell response-recovery prototype has been archived under:

```text
Old/legacy_devvcell_response_recovery_20260505/
```

## Current MVP

The first DevGuard implementation provides a conformal normality baseline:

- unified mouse developmental metadata schema;
- control-only SVD embedding;
- stage-lineage KNN and Mahalanobis nonconformity scores;
- conformal p-values and heldout control false-positive calibration;
- five-class perturbed-cell assignment;
- Developmental Tolerance Index (DTI);
- quick-mode synthetic mouse fixture for code and schema validation.

Quick mode intentionally does not depend on `data/scLine_pro.h5ad` or any
RDEG-derived tables.

## Quick Start

```bash
python scripts/devguard/run_devguard_pipeline.py --mode quick
```

Expected quick-mode outputs:

```text
data/processed/devguard/devguard_quick_mouse.h5ad
results/devguard/normality_reference/
results/devguard/perturbation_classification/
results/devguard/tolerance_index/
results/devguard/figures/
```

## Main Commands

```bash
python scripts/devguard/prepare_mouse_datasets.py --config config/devguard/datasets_mouse.json
python scripts/devguard/build_control_reference.py --config config/devguard/normality_model.json
python scripts/devguard/calibrate_conformal_boundaries.py --config config/devguard/conformal_thresholds.json
python scripts/devguard/classify_perturbed_cells.py --config config/devguard/perturbation_tests.json
python scripts/devguard/compute_developmental_tolerance_index.py
python scripts/devguard/build_devguard_figures.py
```

## Chinese Summary

DevGuard 是一个独立的小鼠发育扰动正常性判定项目。它不使用 RDEG
管线、RDEG 状态节点、OT 正常发育转运、GRN 动力学或矩阵回写结果。
DevGuard 从公开小鼠胚胎和 gastruloid control time-course 数据中学习阶段
和谱系特异的正常性边界，然后将扰动细胞分类为：当前阶段正常、发育延迟、
发育提前、命运偏航或正常图谱外异常。
