# 图件来源记录

## Fig. 1 DevVCell 中文框架图

- 文件：`manuscript/figures/devvcell_framework_cn_ai.png`
- 生成方式：Codex 内置 `image_gen` 深度图像生成工具
- 原始生成文件：`C:\Users\14915\.codex\generated_images\019de963-2ab1-7370-ad64-7bf97a02e14d\ig_0e46ae57636b0e650169f620e41f1481918513e6993c1a1f00.png`
- 项目内复制时间：2026-05-03
- 用途：中文 Nature 风格主图框架图，展示 DevVCell 从胚胎单细胞图谱到细胞状态编码、OT 伪靶标、发育转移算子、TF/GRN 刺激-反馈模块和扰动输出的整体研究路线。

生成 prompt 保存在：

```text
manuscript/figure_prompts/framework_figure_cn.md
```

说明：该图为 AI 辅助生成的科研示意图，正式投稿前应由作者逐项检查标签、箭头含义、模型符号和图注，必要时用矢量绘图软件进行最终排版与校正。

## 2026-05-04 生成的 DevVCell 方法示意图组

- 文件夹：`manuscript/figures/generated_methods/`
- 生成方式：Codex 内置 `image_gen` 图像生成工具
- 原始生成目录：`C:\Users\14915\.codex\generated_images\019dee20-174c-7510-b695-c4ee71305acf`
- 项目内复制时间：2026-05-04
- 用途：补充论文中“总览 + 组件内部机制”的方法示意图组，使读者能够分别理解总框架、transition operator、扰动/恢复读出和正式 evidence pipeline。

项目内文件：

- `generated_methods/devvcell_overview_evidence_system.png`：DevVCell 总览升级图。
- `generated_methods/devvcell_transition_operator_method.png`：细胞级发育转移算子内部机制图。
- `generated_methods/devvcell_perturbation_recovery_method.png`：扰动校准与命运/恢复读出图。
- `generated_methods/devvcell_evidence_pipeline_method.png`：论文级 evidence pipeline 图。

说明：上述图件为 AI 辅助生成的科研示意图，已按“一个总图 + 三个组件图”的结构插入论文。正式投稿前仍建议作者对所有中文标签、数值标注、公式符号和图内小字进行人工校对。
