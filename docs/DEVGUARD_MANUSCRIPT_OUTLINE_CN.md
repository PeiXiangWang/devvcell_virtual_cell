# DevGuard 论文大纲

## 标题

DevGuard: stage-specific developmental normality boundaries for perturbed mouse embryonic cells

## 摘要逻辑

1. 胚胎扰动研究需要判断扰动后细胞是否仍处于正常发育范围。
2. 现有方法常预测表达变化或轨迹，但较少给出统计校准的阶段正常性判定。
3. DevGuard 从 mouse embryo/gastruloid control time-course 中学习 time-lineage normality boundary。
4. 对 perturbation cells，DevGuard 分类为 within-stage normal、delay、acceleration、fate deviation 或 abnormal。
5. 公开小鼠 gastrulation/gastruloid 数据将用于揭示阶段特异的扰动耐受性。
6. Organoid heterogeneity 和 spatial/tomo-seq validation 用于控制自然异质性和支持 fate deviation。

## Results

1. Public mouse developmental datasets enable stage-lineage normality calibration.
2. Conformal boundaries control false-positive rates in heldout control cells.
3. Perturbed cells partition into normal retention, delay, acceleration, fate deviation and abnormal states.
4. Developmental tolerance varies across time, lineage and perturbation.
5. Organoid heterogeneity controls separate biological abnormality from sample variation.
6. Spatial/tomo-seq validates fate-deviation assignments.

## Discussion 边界

正式文章应强调 DevGuard 是独立的问题驱动研究：它不做正常发育推演，不继承 RDEG OT/GRN/rollout/matrix-writeback，而是做扰动后细胞是否仍处于阶段正常范围内的统计校准判定。
