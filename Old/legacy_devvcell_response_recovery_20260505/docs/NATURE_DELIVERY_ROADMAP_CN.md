# DevVCell Nature 级研究包路线图

本文档把用户目标转化为可执行、可验收的研究工程路线。需要明确的是：期刊接收和真实同行评议不能由代码本身保证；本项目的可交付物是投稿前研究包，包括可复现代码、细胞级训练评估、中文论文稿、中文框架图、结果表图、数据与代码可用性说明，以及可供作者进一步提交和答辩的完整材料。

## 2026-05-03 重定位

DevVCell 必须从“RDEG 下游原型”重定位为“发育虚拟细胞”。RDEG 的核心对象是发育图谱、状态节点、跨阶段 OT 边和 GRN-OT 图演化；DevVCell 的核心对象必须是一个可调用的细胞状态模拟器：

```text
input = cell state + developmental stage + system + perturbation
output = next state + perturbed state + fate displacement + recovery cost + rescue candidate
```

新的主张是：DevVCell 是一个 stage-conditioned and perturbation-calibrated developmental virtual cell。它从正常胚胎发育图谱学习细胞级状态转移，用真实 Perturb-seq 数据校准扰动响应，并把模型输出从表达误差扩展为命运偏移、恢复成本和发育能力窗口。

详细重定位方案见 `docs/DEVVCELL_NATURE_VIRTUAL_CELL_STRATEGY_CN.md`。

## 中心科学问题

胚胎发育单细胞图谱提供了不同 Theiler stage 的横截面状态，但不能追踪同一细胞的真实未来。本项目拟构建一个发育虚拟细胞框架 DevVCell：以正常发育轨迹为主轴，用最优传输构造相邻阶段弱监督，用细胞级 transition operator 学习发育状态推进，并用真实 Perturb-seq 与 TF/GRN/stimulus 模块估计扰动响应、命运偏移和反馈恢复能力。

## 预期创新点

1. **发育时间条件化虚拟细胞**：把 Theiler stage 和 system context 作为细胞状态动力学的一部分，而不是事后分组标签。
2. **正常发育自监督 + 真实扰动校准**：用相邻发育阶段 OT pseudo-target 学习自然转移，用外部 Perturb-seq 学习 perturbation-conditioned residual。
3. **命运偏移与恢复成本**：从单纯表达 MSE 扩展到 fate displacement、recovery cost 和 rescue candidate。
4. **发育能力窗口**：把 TS14--TS19 等窗口定义为 competence window，检验同一 TF 在不同阶段的可塑性和可恢复性。
5. **抗偏差评估**：正面对比 perturbed mean、matched mean、ridge residual、CellOT、GEARS、scGen/CPA、foundation model embedding，并采用 Systema-style perturbation-specific effect 评估。

## 主要实验路线

### 阶段 1：数据与原型结果冻结

- 固化 `data/scLine_pro.h5ad` 的结构、细胞数、基因数、obs/var 字段和稀疏矩阵信息。
- 固化 `data/rdeg_neural_cell_mvp/` 中 RDEG 派生表。
- 重新运行 `scripts/devvcell_lite.py` 和 `scripts/ot_grn_developmental_impact.py`，生成 prototype/proxy 表和图。

### 阶段 2：细胞级数据子集

- 使用 `scripts/export_cell_level_subset.py` 从完整 H5AD 导出平衡子集。
- 初始系统：neural、mesoderm/muscle、erythroid。
- 初始阶段：Theiler stage 12-19。
- 每个 system-stage 最多采样 800 个细胞，保留全部 3,000 个基因和 donor/stage/cell-type metadata。

### 阶段 3：细胞级 transition 基线

- 对表达矩阵进行 sparse TruncatedSVD 编码。
- 构建同系统相邻 Theiler stage 的 Sinkhorn OT barycentric pseudo-targets。
- 训练并比较 identity、mean-shift、ridge、MLP transition。
- 以 heldout target stages 评估 pair latent MSE、centroid MSE 和 RBF-MMD。

### 阶段 4：正式模型扩展

后续版本需要把当前基线扩展为以下模块：

- cell encoder：PCA/SVD baseline -> scVI or sparse autoencoder。
- OT coupling：扩大 mini-batch Sinkhorn coupling 的规模并加入跨 batch 稳定性评估。
- transition operator：stage-conditioned residual dynamics。
- stimulus head：已实现 zero-shot TF knockdown GRN projection；下一步加入可训练 stimulus-conditioned residual dynamics。
- GRN feedback head：把 learned GRN edge weight 注入 response/recovery prediction。
- perturbation calibration：把 scPerturb/Perturb-seq 真实扰动标签作为校准分支，而不是只做外部补充实验。
- recovery/rescue head：输出最小恢复成本和候选 rescue TF/pathway。
- decoder：从 latent state 重建表达或模块分数。
- ablation：no OT、no GRN、no stimulus、no feedback、identity dynamics。

### 阶段 5：论文与图件

- 中文 Nature 风格论文：`manuscript/main_cn.tex`。
- 中文 AI 生成框架图 prompt：`manuscript/figure_prompts/framework_figure_cn.md`。
- 主图建议：
  - Fig. 1：DevVCell virtual-cell input/output framework。
  - Fig. 2：heldout stage/donor/system cell-state transition benchmark。
  - Fig. 3：developmental competence window and vulnerability-response coupling。
  - Fig. 4：TF perturbation virtual screen with fate displacement and recovery cost。
  - Fig. 5：external Perturb-seq calibration and peer benchmark。
  - Fig. 6：stage-dependent CRISPRi/a or Perturb-seq validation。

## 当前验收状态

- 已有 prototype/proxy 表和图。
- 已新增 cell-level 子集导出配置和脚本。
- 已新增 cell-level transition baseline 训练评估入口，包含 Sinkhorn OT 伪靶标。
- 已新增 stage/system-conditioned context residual MLP，作为显式条件化残差发育转移算子；主 heldout benchmark 中 pair MSE 优于 ridge 和普通 MLP。
- 已新增 transition bootstrap 统计比较，输出模型均值、95% CI 和 context residual 相对 baselines 的配对改进。
- 已新增 TF/GRN stimulus response head，输出细胞级 latent response、fate displacement 和 recovery proxy。
- 已新增多 seed 和组件消融套件，覆盖 transition pairing、ridge/identity/mean-shift baseline，以及 TF/GRN stimulus 的全局边/系统边消融。
- 已新增统一复现入口 `scripts/run_reproducible_pipeline.py`，可生成 `results/reproducibility_manifest.json` 并编译中文论文。
- 已新增中文论文骨架、AI 生成中文框架图、参考文献骨架、LLM 使用声明和同行方法对比文档 `docs/PEER_METHOD_COMPARISON_CN.md`。
- 已新增外部扰动 benchmark：`scripts/prepare_external_perturbation_benchmark.py`、`scripts/run_external_perturbation_benchmark.py` 和 `docs/EXTERNAL_PERTURBATION_BENCHMARK_CN.md`，已在 scPerturb Datlinger/Bock 2021 数据上完成 guide-transfer 首版评估。

## 风险与必须补强的证据

- 当前 TF perturbation 仍是计算 proxy，不等同真实实验扰动。
- 当前 pseudo-targets 由横截面分布和 OT 构造，不是真实 lineage tracing。
- Nature 水准需要更强外部验证：至少 3 个外部 Perturb-seq 数据集、直接同行基线、独立发育图谱或湿实验验证。
- 需要完整统计报告：n 值、误差定义、置信区间、消融、随机种子重复。
- 当前已具备第一版多 seed、消融表和一个外部 perturbation benchmark，但仍缺少 CellOT/GEARS/scGen/CPA/Systema-style 直接对比。
- 必须补 heldout donor、heldout system、heldout gene、heldout context 和 heldout cell state。
- 若目标是 Nature 主刊或强子刊，至少需要一个阶段依赖扰动实验验证 competence window。
- 使用 LLM 生成图件或文本辅助时，必须在 Methods 或声明中透明记录；LLM 不能作为作者。
