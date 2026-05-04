# DevVCell：Nature 级发育虚拟细胞重定位方案

更新日期：2026-05-03

## 核心判断

DevVCell 不能再被写成 RDEG 的“下游分析”或“轻量版复用”。它必须被定义为一个发育虚拟细胞系统：给定一个细胞状态、发育时间、所属系统和外部扰动，模型预测这个细胞的下一阶段状态、扰动响应方向、命运偏移、恢复成本和可救援性。

这和 RDEG 的区别必须写在论文第一页。RDEG 回答“胚胎发育图谱如何跨阶段组织和重排”；DevVCell 回答“一个细胞在发育时间和扰动条件下会变成什么”。前者是图谱演化框架，后者是可扰动的细胞状态模拟器。

## 一句话论文主张

DevVCell turns cross-sectional embryonic atlases into a perturbable developmental virtual cell by learning stage-conditioned cell-state dynamics, calibrating perturbation responses with real Perturb-seq data, and quantifying fate displacement and recovery cost across developmental competence windows.

中文表述：

DevVCell 将横截面胚胎单细胞图谱转化为一个可扰动的发育虚拟细胞：它学习阶段条件化的细胞状态动力学，用真实 Perturb-seq 数据校准扰动响应，并在发育能力窗口内量化命运偏移与恢复成本。

## 必须建立的五个创新点

### 1. 发育时间条件化的虚拟细胞

同行很多方法预测扰动后的表达，但多数把细胞看成静态输入。DevVCell 的第一创新点是把 Theiler stage 作为细胞状态动力学的一部分，而不是一个事后分组标签。

当前已有基础：

- `context_residual_mlp` 使用 source stage、target stage、stage delta 和 system one-hot。
- 主 heldout 结果中，`context_residual_mlp` mean pair latent MSE 为 0.3139，优于 ridge 的 0.3254、普通 MLP 的 0.3302、identity 的 0.3813。
- 配对 bootstrap 显示相对 ridge 的平均 MSE 改进为 0.0115，95% CI 为 0.0014--0.0211。

Nature 级补强：

- 把 stage 从简单数字扩展为 developmental clock embedding。
- 加入 stage uncertainty：同一 author day、donor、Theiler stage 的时间不确定性需要显式建模。
- 做 heldout stage、heldout donor、heldout system 三类泛化，不只做 heldout target stage。

### 2. 正常发育自监督 + 真实扰动校准的双任务虚拟细胞

DevVCell 不能只靠 OT pseudo-target，否则会被认为是 Waddington-OT 或 moscot 的小变体。必须把正常发育和真实扰动变成两个互相校准的训练任务：

- 正常发育任务：用相邻 Theiler stage 的 Sinkhorn barycentric pseudo-target 学习自然转移。
- 扰动任务：用 scPerturb / Perturb-seq 中 control-to-perturbed pseudo-target 学习 perturbation-conditioned residual。
- 联合任务：要求同一个 latent space 同时解释“正常发育推进”和“扰动导致的偏离”。

当前已有基础：

- `scripts/train_cell_transition_baseline.py` 已完成正常发育 cell-level transition baseline。
- `scripts/run_external_perturbation_benchmark.py` 已完成 Datlinger/Bock 2021 guide-transfer benchmark。
- 外部 benchmark 中，`gene_context_ridge_residual` heldout guide mean latent MSE 为 0.2209，优于 identity/global shift 约 0.789。

Nature 级补强：

- 将外部扰动 benchmark 从独立脚本升级为正式训练/验证分支。
- 统一 development transition 和 perturbation transition 的接口。
- 至少接入 3 个外部 perturb-seq 数据集，并包含 heldout gene、heldout guide、heldout context 和 heldout cell state。

### 3. 命运偏移和恢复成本，而不仅是表达 MSE

虚拟细胞论文不能只报告 latent MSE。真正有生物价值的问题是：扰动是否把细胞推离正常发育轨道？推到哪个命运附近？能否被 rescue？何时最脆弱？

DevVCell 的核心输出应固定为五个量：

- `developmental_step_prediction`：正常下一阶段状态。
- `perturbed_state_prediction`：扰动后状态。
- `fate_displacement`：扰动状态相对正常轨迹的偏移。
- `recovery_cost`：回到正常轨迹所需的最小调节强度。
- `rescue_candidate`：最可能降低 recovery cost 的 TF/pathway。

当前已有基础：

- `devvcell_lite.py` 已有 stage vulnerability、perturbation priority、recovery probability proxy。
- `run_stimulus_response_head.py` 已有 TF/GRN stimulus response、fate displacement 和 recovery proxy。

Nature 级补强：

- 把 proxy 改为可检验指标：用真实 Perturb-seq 或体外分化扰动数据校准。
- 对每个 top TF 报告 response amplitude、fate shift、recovery cost、rescue target 和 uncertainty。
- 加入剂量曲线：CRISPRi/a knockdown strength 或药物 dose 不应只是一档。

### 4. 发育能力窗口 competence window

这是 DevVCell 最应该和同行拉开距离的地方。GEARS、CellOT、scGen、CPA、scGPT 多数关注“扰动基因到表达响应”；DevVCell 应该关注“同一个扰动在不同发育阶段为什么效果不同”。

论文应把 TS14--TS19 写成第一版 competence window：

- 同一个 TF 在 TS14、TS16、TS18 的响应强度不同。
- 同一个扰动在 neural、mesoderm/muscle、erythroid 中命运偏移不同。
- 高 vulnerability stage 应该对 perturbation 更敏感，且 recovery cost 更高。

当前已有基础：

- 配置中已有 `stage_window_of_interest = TS14-TS19`。
- stage vulnerability 已整合 temporal sensitivity、rollout error、module shift、outlier ratio 和 next-step MSE。
- stimulus head 已按 TF--system--stage 输出响应。

Nature 级补强：

- 对 TS14--TS19 做系统性 perturbation simulation，而不是只展示 top TF。
- 设计真实实验：同一 TF 在窗口内外分别扰动，验证响应强度和恢复概率差异。
- 明确 competence window 不是 RDEG 的“阶段脆弱性图”，而是 DevVCell 的“虚拟细胞可塑性窗口”。

### 5. 抗偏差评估，而不是追求漂亮深度模型

2025 年扰动预测领域的一个核心教训是：复杂深度模型经常没有明显超过简单线性或平均效应基线，常用指标还会被 systematic variation 高估。DevVCell 如果想冲 Nature 子刊，必须主动采用这种更严格的评估框架。

必须加入的基线：

- identity。
- global mean shift。
- perturbed mean。
- matched mean。
- ridge residual。
- CellOT。
- GEARS。
- scGen / CPA。
- scGPT 或 scFoundation embedding + linear head。

必须加入的评估：

- Pearson delta 和 RMSE 只是基础指标。
- top DE gene recovery。
- effect cosine。
- perturbation-specific residual after removing systematic effect。
- calibration curve。
- heldout gene family 和 heterogeneous gene panel。
- batch、guide efficiency、cell cycle、stress response 控制。

这部分会让 DevVCell 的论文姿态更强：不是声称“深度模型天然更好”，而是证明“发育时间、系统上下文和恢复目标提供了线性基线无法解释的增量”。

## 和同行的清晰差异

| 方法 | 同行核心能力 | DevVCell 必须打出的区别 |
|---|---|---|
| Waddington-OT | 从未配对时间点恢复群体发育轨迹 | DevVCell 不止输出 coupling，而是训练可调用的细胞状态转移函数和扰动响应函数。 |
| moscot | atlas-scale、多模态、时空 OT 映射 | DevVCell 不和 moscot 拼 OT 工程规模，而是把 OT 作为虚拟细胞的监督信号，并预测扰动后的 fate/recovery。 |
| CellRank | RNA velocity / transition kernel / fate probability | DevVCell 不依赖 velocity，强调 stage-conditioned counterfactual dynamics。 |
| CellOracle | 用 GRN 模拟 TF 扰动后的 cell identity vector | DevVCell 输出高维状态转移、发育窗口、恢复成本和真实 Perturb-seq 校准，而不是只输出 2D identity vector。 |
| CellOT | 用 neural OT 学习 control-to-perturbed response | DevVCell 同时学习 normal development 和 perturbation response，并把扰动效果放回发育时间轴。 |
| GEARS | 图先验驱动的单/多基因扰动表达预测 | DevVCell 的差异是 stage competence、fate displacement 和 recovery/rescue，而不是单纯拼 perturbation MSE。 |
| scGen / CPA | 条件生成式扰动响应 | DevVCell 应强调发育时间条件和恢复目标；同时必须直接 benchmark，不能只口头比较。 |
| scGPT / scFoundation / Nicheformer | 大规模单细胞 foundation model 表征 | DevVCell 不做通用 foundation model，而做发育和扰动任务专用的 mechanistic virtual cell；foundation model 可作为 encoder baseline。 |
| Systema | 去除 systematic variation 后评估 perturbation-specific effects | DevVCell 应采用 Systema 式严格评估，证明不是复制平均扰动效应。 |

## Nature 级主图设计

### Figure 1：DevVCell 总框架

核心信息：

- 输入：cell state、stage、system、donor、perturbation。
- 双任务：normal development transition + perturbation transition。
- 输出：next state、perturbed state、fate displacement、recovery cost、rescue candidate。
- 明确和 RDEG 的接口：RDEG 只作为 atlas prior / teacher evidence，不是 DevVCell 本体。

### Figure 2：发育虚拟细胞能预测 heldout stage dynamics

必须展示：

- heldout stage/donor/system 切分。
- context residual dynamics 优于 identity、mean shift、ridge、MLP。
- Sinkhorn pseudo-target 和 nearest-target 的消融。
- uncertainty 或 bootstrap CI。

当前可用结果：

- 主 heldout benchmark。
- transition bootstrap CI。
- 多 seed ablation。

缺口：

- donor split。
- system split。
- 更大 stage 窗口。

### Figure 3：competence window 和发育可塑性

必须展示：

- TS14--TS19 窗口的 vulnerability。
- 同一 TF 在不同 stage 的 response curve。
- window 内外 recovery cost 对比。
- 高 vulnerability 是否预测高 perturbation sensitivity。

当前可用结果：

- stage vulnerability。
- TF--system--stage stimulus response。

缺口：

- 真实 perturbation 或独立数据校准。

### Figure 4：TF/GRN perturbation virtual screen

必须展示：

- top TF 的 response amplitude、fate displacement、recovery cost。
- rescue candidate network。
- GRN component ablation：full GRN、no global GRN、no system edge。

当前可用结果：

- stimulus response head。
- stimulus ablation。

缺口：

- rescue head 目前仍是 proxy，需要可训练或实验校准。

### Figure 5：外部 Perturb-seq 校准

必须展示：

- scPerturb 至少 3 个数据集。
- heldout gene/guide/context。
- DevVCell 与 CellOT、GEARS、scGen/CPA、linear baselines 对比。
- Systema-style perturbation-specific effect 评估。

当前可用结果：

- Datlinger/Bock 2021 guide-transfer 首版。

缺口：

- 直接运行同行方法。
- 多数据集。
- Systema-style 去偏指标。

### Figure 6：真实实验验证

最低要求：

- 选择 2--3 个 TF。
- 在窗口内外做 CRISPRi/a 或 Perturb-seq。
- 验证 response amplitude、fate marker shift、recovery/rescue。

更强版本：

- 同一 TF 在 TS14、TS16、TS18 的扰动效果呈现阶段依赖。
- DevVCell 能预测哪个阶段最容易偏离正常命运，哪个 rescue target 可降低偏移。

## 当前结果如何写成效果

### 已经可以写的效果

- 细胞级 transition：`context_residual_mlp` 是当前 heldout pair MSE 最优模型。
- 统计稳定性：bootstrap 和多 seed 消融支持该趋势，但 strata 仍少。
- TF/GRN stimulus：能生成 TF--system--stage 级别的 response/recovery hypothesis。
- 外部扰动：guide-transfer 上 gene/context residual ridge 显著优于 identity 和 global shift。
- 工程复现：已有一键 pipeline、manifest、论文 PDF 和关键输出。

### 不能过度写的效果

- 不能说已经实现真实虚拟细胞。
- 不能说 TF perturbation 预测已被生物实验验证。
- 不能说优于 CellOT、GEARS、CPA、scGPT，因为当前还没有直接运行这些基线。
- 不能把 RDEG 的节点图谱结果写成 DevVCell 的核心结果。

### 应该改成的强表述

当前版本：

“DevVCell-lite 是一个发育虚拟细胞原型。”

投稿版本：

“DevVCell is a stage-conditioned and perturbation-calibrated developmental virtual cell that learns normal developmental dynamics from embryonic atlases and estimates counterfactual perturbation responses with explicit fate-displacement and recovery-cost readouts.”

## 必须补齐的计算实验

### A. 发育泛化

- heldout target stage：已有。
- heldout source-target transition：需要补。
- heldout donor：必须补。
- heldout system：必须补。
- fine-stage split：建议补。
- full TS12--TS27 而非 TS12--TS19：必须补。

### B. 扰动泛化

- heldout guide：已有首版。
- heldout gene：必须补。
- heldout gene family：建议补。
- heldout context：必须补。
- heldout cell type/state：必须补。
- unseen combination：如果要和 GEARS 正面对比，必须补。

### C. 同行基线

最低可投稿计算基线：

- identity。
- mean shift / perturbed mean / matched mean。
- ridge residual。
- MLP residual。
- CellOT。
- scGen 或 CPA。
- GEARS。
- scGPT/scFoundation embedding + ridge。

如果时间有限，优先级是：

1. Systema-style simple baselines。
2. CellOT。
3. GEARS。
4. scGen/CPA。
5. foundation model embedding。

### D. 评估指标

每个扰动任务必须报告：

- pair latent MSE。
- centroid MSE。
- effect cosine。
- Pearson delta。
- top DE gene recovery。
- RBF-MMD。
- perturbation-specific effect after systematic-effect removal。
- calibration / uncertainty。

每个发育任务必须报告：

- next-stage prediction MSE。
- stage identity rank。
- fate-neighborhood accuracy。
- system composition shift。
- donor robustness。
- uncertainty across seeds。

## 必须补齐的湿实验

Nature 主刊或强子刊通常不能只靠计算 proxy。建议把实验目标收缩到“验证发育能力窗口”，这是 DevVCell 最独特的生物问题。

### 实验 1：阶段依赖扰动

目的：验证同一 TF 在不同发育阶段的响应不同。

设计：

- 选择 2 个 TF：一个 neural high-response，一个 mesoderm/erythroid high-response。
- 在 TS14、TS16、TS18 对应体外分化时间点做 CRISPRi/a 或 siRNA。
- readout：scRNA-seq 或 targeted Perturb-seq。
- 验证：DevVCell predicted response amplitude 与真实 delta expression / marker shift 相关。

### 实验 2：恢复/救援

目的：验证 recovery cost 和 rescue candidate。

设计：

- 对最高 recovery cost 的 TF knockdown 做 rescue TF/pathway co-perturbation。
- 比较单扰动和 rescue 后的 fate displacement。
- readout：关键 marker、细胞组成、trajectory affinity。

### 实验 3：窗口内外对照

目的：证明 competence window 是真实生物现象，不是模型指标。

设计：

- 同一 perturbation 在窗口内 TS14--TS19 和窗口外 stage 做对照。
- 预注册指标：response amplitude、normal trajectory distance、fate marker displacement、recovery probability。

## 审稿人最可能攻击的问题

### 问题 1：这不就是 RDEG 的下游包装吗？

回答策略：

- RDEG 预测 atlas graph；DevVCell 预测 individual cell state。
- RDEG 的输出是 prior/teacher evidence；DevVCell 的核心是 cell-level callable operator。
- DevVCell 的外部 Perturb-seq benchmark 不依赖 RDEG。

### 问题 2：OT pseudo-target 不是真实 lineage。

回答策略：

- 明确称为 weak supervision / distributional pseudo-target。
- 加入 heldout donor、heldout stage、external perturbation 和实验验证。
- 不把 OT coupling 写成真实细胞命运，只写成训练监督。

### 问题 3：深度模型没有超过简单线性基线。

回答策略：

- 主动纳入 perturbed mean、matched mean、ridge residual。
- 强调 DevVCell 的增量不是“深度”，而是 stage competence + fate/recovery readout。
- 使用 Systema-style 去偏评估。

### 问题 4：虚拟细胞这个词过大。

回答策略：

- 不声称通用 virtual cell。
- 精确定义为 developmental virtual cell for embryonic transcriptomic state transitions。
- 输出范围限定为 transcriptomic state、fate displacement、recovery cost，不涉及完整生理细胞。

### 问题 5：缺少真实实验。

回答策略：

- 至少完成一个阶段依赖扰动实验。
- 把实验设计和预测预注册。
- 不把未验证预测写成发现，只写成 candidate hypotheses。

## 最终投稿定位

### Nature 主刊

需要同时满足：

- 一个非常清晰的新生物学发现：例如“发育能力窗口决定 TF 扰动可恢复性”。
- 多系统/多数据集/多物种或空间验证。
- 至少一个关键湿实验验证。
- 方法本身可复现并明显超越强基线。

当前距离：较远，需要实验验证和更强 benchmark。

### Nature Methods

需要同时满足：

- 方法学定义清楚。
- 相比 CellOT、GEARS、CellOracle、scGen/CPA、foundation model baseline 有明确增量。
- 软件可复现、可扩展、接口稳定。
- 多数据集 benchmark 严格。

当前距离：中等，最现实。

### Nature Biotechnology / Nature Computational Science

需要同时满足：

- perturbation prediction 和 cell engineering 应用价值明确。
- 有真实 Perturb-seq 或 CRISPR 验证。
- 能提出可执行实验设计并节省筛选成本。

当前距离：中等偏远，取决于实验。

### Nature Communications / Cell Systems

需要同时满足：

- 计算主线完整。
- 有外部验证。
- 创新点清楚，和 RDEG 明确区分。

当前距离：较现实。

## 下一步执行清单

### 立即修改

- 把论文题目改为“stage-conditioned perturbable developmental virtual cell”。
- 摘要删除过多 RDEG proxy 语言，突出 cell-level operator、perturbation calibration、fate/recovery。
- 同行对比加入 GEARS、CellOracle、Systema、scFoundation、Nicheformer 和 2025 negative baseline。
- 图 1 改成 virtual cell input-output，不再像 RDEG 流程图。

### 一周内

- 扩展外部 perturbation benchmark 至至少 2 个 scPerturb 数据集。
- 加入 perturbed mean、matched mean、ridge residual 的 Systema-style baseline。
- 增加 heldout donor split。
- 将 stimulus head 的 outputs 改成统一 `fate_displacement`、`recovery_cost`、`rescue_candidate`。

### 一个月内

- 跑 CellOT、GEARS、CPA/scGen baseline。
- 做 full TS12--TS27 cell-level subset 或分批训练。
- 增加 scVI/scGPT/scFoundation embedding 对比。
- 设计并启动 1 个 TF 阶段依赖扰动实验。

### 投稿前

- 完成至少一个真实扰动验证。
- 冻结公开代码和数据 manifest。
- 写 Methods 里的 bias-control 和 systematic variation section。
- 将所有 proxy 语言替换为 calibrated / validated / hypothesis，避免过度宣称。

## 结论

DevVCell 要冲 Nature 及其子刊，不能靠“RDEG + OT + GRN”的旧叙事。它必须变成一个真正的发育虚拟细胞论文：细胞级、阶段条件化、可扰动、可恢复、可用真实 Perturb-seq 校准，并且在强线性基线和同行方法面前仍然显示增量。

最强的差异化不是“我们也能预测扰动表达”，而是：

同一个扰动在不同发育时间和不同系统中具有不同的命运偏移与恢复成本；DevVCell 是第一个把这种发育能力窗口系统化建模、预测和验证的虚拟细胞框架。
