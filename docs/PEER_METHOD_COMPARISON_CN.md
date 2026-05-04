# DevVCell 与同行方法的定位

## 新定位

DevVCell 的同行坐标不能只放在 RDEG、OT 或 GRN 推断里。为了成为“虚拟细胞”方向的论文，它必须和三类方法正面对比：

1. 发育/时空映射方法：Waddington-OT、moscot、CellRank。
2. 扰动响应方法：CellOT、GEARS、CellOracle、scGen/CPA、GPerturb。
3. 单细胞 foundation model 与评估框架：scGPT、scFoundation、Nicheformer、Systema，以及 2025 年关于深度扰动模型不一定超过简单线性基线的负面结果。

DevVCell 的差异化主张应固定为：

```text
stage-conditioned + perturbation-calibrated + fate/recovery readout
```

也就是说，DevVCell 不和 moscot 拼 OT 工程规模，不和 scGPT 拼预训练规模，不和 GEARS 只拼 perturbation MSE；它要证明“发育时间、系统上下文和恢复成本”能带来同行模型没有覆盖的虚拟细胞能力。

## 对比对象

| 方法 | 代表文献 | 核心能力 | DevVCell 当前差异化 |
|---|---|---|---|
| Waddington-OT | Schiebinger et al., Cell 2019, doi:10.1016/j.cell.2019.01.006 | 从未配对时间点分布中用最优传输恢复发育轨迹 | DevVCell 继承 OT 伪靶标思想，但进一步训练显式细胞级 transition operator，并加入 TF/GRN stimulus head。 |
| moscot | Klein et al., Nature 2025, doi:10.1038/s41586-024-08453-2 | 面向大规模单细胞、多模态、空间和时间映射的 OT 框架 | DevVCell 当前规模较小，但更聚焦发育虚拟细胞、扰动响应和恢复 proxy；后续应吸收 moscot 的可扩展 OT 和多模态设计。 |
| CellRank | Lange et al., Nature Methods 2022, doi:10.1038/s41592-021-01346-6 | 构造 directed transition kernel 并估计 fate probabilities | DevVCell 当前不依赖 RNA velocity，而是用 Theiler stage 分布和 OT pseudo-target 训练跨阶段算子。 |
| CellOracle | Kamimoto et al., Nature 2023, doi:10.1038/s41586-022-05688-9 | 从 GRN 推断出发做 in silico TF perturbation，并输出 cell identity transition vector | DevVCell 应输出高维细胞状态转移、stage competence、fate displacement 和 recovery cost，而不是只在 2D identity map 上模拟方向。 |
| CellOT | Bunne et al., Nature Methods 2023, doi:10.1038/s41592-023-01969-x | 学习 control 与 perturbation 分布之间的 neural OT map | DevVCell 当前 zero-shot TF/GRN stimulus 未使用真实 perturbation labels；未来应在外部 perturb-seq 上与 CellOT 直接 benchmark。 |
| GEARS | Roohani et al., Nature Biotechnology 2024, doi:10.1038/s41587-023-01905-6 | 用 gene graph 和 GNN 预测单基因/多基因 perturbation expression | DevVCell 不应只和 GEARS 比 MSE，而应证明 stage competence、fate displacement、recovery/rescue 是 GEARS 未覆盖的发育虚拟细胞输出。 |
| GPerturb | Xing and Yau, Nature Communications 2025, doi:10.1038/s41467-025-61165-7 | 用 Gaussian process 建模单细胞扰动数据，强调不确定性 | DevVCell 后续应吸收 uncertainty/calibration 报告，避免只给点预测。 |
| scGen | Lotfollahi et al., Nature Methods 2019, doi:10.1038/s41592-019-0494-8 | 用生成模型预测跨 cell type/study/species perturbation response | DevVCell 的扰动层目前更强调发育阶段、GRN evidence 和 recovery proxy。 |
| CPA | Lotfollahi et al., bioRxiv 2021 / 后续扰动模型基线 | 组合式 perturbation autoencoder，常用于 drug/genetic perturbation response | DevVCell 必须直接比较 CPA 或 scGen/CPA 类模型，证明 stage/system context 的增量。 |
| scVI | Lopez et al., Nature Methods 2018, doi:10.1038/s41592-018-0229-2 | 单细胞深度生成 latent 表示和批次建模 | DevVCell 当前用 sparse SVD baseline；后续可把 scVI 作为 encoder baseline 或替代模块。 |
| scGPT | Cui et al., Nature Methods 2024, doi:10.1038/s41592-024-02201-0 | 单细胞 foundation model，用生成式预训练支持多任务 | DevVCell 目前不是 foundation model；可在后续版本比较 scGPT embedding 是否提升 transition/stimulus 泛化。 |
| scFoundation | Hao et al., Nature Methods 2024, doi:10.1038/s41592-024-02305-7 | 大规模单细胞 foundation model，报告多种下游任务包括 perturbation prediction | DevVCell 可把 scFoundation embedding 作为强 encoder baseline，而不是声称自建 foundation model。 |
| Nicheformer | Tejada-Lapuerta et al., Nature Methods 2025, doi:10.1038/s41592-025-02814-z | 同时学习 dissociated single-cell 和 spatial context 的 foundation model | DevVCell 若未来接入空间胚胎数据，可把 Nicheformer 作为空间上下文 baseline。 |
| Systema | Viñas Torné et al., Nature Biotechnology 2025, doi:10.1038/s41587-025-02777-8 | 扰动响应预测评估框架，强调去除系统性偏差后的 perturbation-specific effects | DevVCell 当前应借鉴其评估原则，避免仅用全局差异或混杂信号高估扰动预测能力。 |
| Linear / perturbed mean / matched mean baselines | Ahlmann-Eltze et al., Nature Methods 2025, doi:10.1038/s41592-025-02772-6；Systema | 显示深度扰动模型常常不明显超过简单线性或平均效应基线 | DevVCell 必须主动纳入这些强基线，证明增量来自发育时间和恢复目标，而不是模型复杂度。 |
| scPerturb 数据资源 | Peidli et al., Nature Methods 2024, doi:10.1038/s41592-023-02144-y | 统一整理多种单细胞扰动数据集，提供 perturb-seq/CRISPR H5AD 入口 | DevVCell 已接入 Datlinger/Bock 2021 scPerturb H5AD，并完成 guide `_1` 训练、guide `_2` heldout 的外部扰动 benchmark。 |

## 当前新增方法层

本轮新增 `context_residual_mlp`：一个 stage/system-conditioned residual transition operator。其形式为：

```text
z_{t+1} = z_t + f_theta(z_t, source_stage, target_stage, system)
```

其中 `z_t` 是细胞 latent state，context 包括 source stage、target stage、stage delta 和 system one-hot。该设计相较普通 MLP 的优势是：

- 明确建模非自治发育动力学，而不是假设所有阶段共享同一个无条件映射。
- 可以把 heldout stage 推断表述为条件化泛化任务。
- 保留 residual form，使模型学习发育推进方向，而不是重新生成整个状态。
- 与 CellOT 的 learned map 思想相邻，但当前监督来自 stage OT pseudo-target，而不是真实 perturbation labels。

更强的论文表述应当是：`context_residual_mlp` 不是为了展示“深度模型更复杂”，而是为了检验发育虚拟细胞的关键假设：同一细胞状态在不同发育时间和不同系统中具有不同可达方向，因此 transition operator 必须是 stage/system-conditioned。

## 当前证据

主模型 heldout 结果：

- `context_residual_mlp` mean pair latent MSE：0.3139
- `ridge` mean pair latent MSE：0.3254
- `mlp` mean pair latent MSE：0.3302
- `identity` mean pair latent MSE：0.3813

bootstrap 配对差异显示，`context_residual_mlp` 相对 ridge 的平均 MSE 改进为 0.0115，95% CI 为 0.0014--0.0211；相对普通 MLP 的平均 MSE 改进为 0.0163，95% CI 为 0.0058--0.0302。

多 seed 消融中，`context_residual_mlp` 在 Sinkhorn OT 伪靶标下平均 MSE 为 0.2974，优于 ridge 的 0.3094；在 nearest-neighbor 伪靶标下 ridge 的 pair MSE 更低，但 `context_residual_mlp` 的 centroid MSE 和 RBF-MMD 更优。

外部扰动首版结果：

- scPerturb Datlinger/Bock 2021 guide-transfer benchmark 中，`gene_context_ridge_residual` mean pair latent MSE 为 0.2209。
- identity 和 global mean shift 约为 0.789。
- 这说明显式 perturbation gene + stimulation context 有真实扰动数据上的增量，但当前还没有直接和 CellOT、GEARS、CPA/scGen、Systema baselines 比较。

## Nature 级补强方向

1. 外部 perturbation benchmark：扩展到至少 3 个 scPerturb 数据集；除了 heldout guide，还要做 heldout gene、heldout context 和 heldout cell state。
2. 强基线正面对比：至少加入 perturbed mean、matched mean、ridge residual、CellOT、GEARS、scGen/CPA；foundation model embedding 可作为 encoder baseline。
3. 可扩展 OT：吸收 moscot 的 mini-batch、多模态和时空映射设计，但 DevVCell 的重点仍是 virtual-cell operator，而不是 OT 软件框架。
4. 表示学习对比：SVD vs scVI vs scGPT/scFoundation embedding 对 transition 和 stimulus 泛化的贡献。
5. 偏差控制：借鉴 Systema，报告 perturbation-specific effects，避免系统性 batch/cell composition 差异被误读为扰动预测能力。
6. 发育能力窗口：把 TS14--TS19 写成 competence window，验证同一 TF 在窗口内外具有不同 response amplitude、fate displacement 和 recovery cost。
7. 生物验证：把 top TF 候选转化为 stage/system-specific CRISPRi/a 或 Perturb-seq 实验，并预注册 response、fate shift、recovery/rescue 指标。
