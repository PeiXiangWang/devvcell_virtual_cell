# DevGuard Tal1/T 稳健性与外部控制扩展

版本：2026-05-05

## 本轮新增产物

新增脚本：

- `scripts/devguard/analyze_chimera_robustness.py`
- `scripts/devguard/analyze_tal1_marker_modules.py`
- `scripts/devguard/analyze_gse212050_organoid_heterogeneity.py`
- `scripts/devguard/map_gse123187_spatial_tomo_lineages.py`

增强脚本：

- `scripts/devguard/build_gse123187_h5ad_from_raw_tar.py` 现在支持 `--mode-filter`、`--member-pattern`，并为 tomo-seq segment 写入 `tomo_axis_fraction` 和 `tomo_axis_bin`。

主要输出目录：

- `results/devguard_real/chimera_robustness/`
- `results/devguard_real/MouseGastrulationData_tal1_chimera_full_integrated_e85_strict/marker_modules/`
- `results/devguard_real/GSE212050_strict_sample/heldout_control_classification/`
- `results/devguard_real/GSE212050_strict_sample/organoid_heterogeneity_control/`
- `data/processed/devguard/GSE123187_tomo_3files.h5ad`
- `results/devguard_real/GSE123187_tomo_spatial_lineage_mapping/`

## Tal1/T sample-level robustness

使用 integrated chimera E8.5 matched-control reference 的既有 cell-level classification，按 `sample_id` 做 sample-level bootstrap。

| cohort | E8.5 classified cells | sample units | within-stage normal | fate deviation | abnormal/off-normal |
|---|---:|---:|---:|---:|---:|
| Control heldout | 19,964 | 12 | 0.9548 | 0.0338 | 0.0114 |
| T chimera | 15,817 | 4 | 0.9434 | 0.0426 | 0.0140 |
| Tal1 chimera | 28,305 | 2 | 0.5758 | 0.2620 | 0.1623 |

Tal1 的 sample-level bootstrap 区间仍很窄：fate deviation 0.2569-0.2666，abnormal/off-normal 0.1546-0.1693。T chimera 接近 heldout control：fate deviation 0.0392-0.0459，abnormal/off-normal 0.0085-0.0217。注意 Tal1 E8.5 tomato-positive 这里只有 2 个 sample units，因此 sample-level 稳健性支持效应强，但正式论文仍应把 Tal1 sample 数作为限制写清楚。

## High-FPR lineage 过滤/降权

Integrated reference 中 high-FPR lineages 为 4/28：

| lineage | heldout FPR |
|---|---:|
| Endothelium | 0.1946 |
| NMP | 0.1885 |
| Neural crest | 0.1366 |
| Doublet | 0.1324 |

对当前 lineage reference 做过滤或按 `min(1, alpha / heldout_fpr)` 降权后，Tal1 的异常信号仍存在：

| sensitivity scenario | Tal1 within | Tal1 fate | Tal1 abnormal |
|---|---:|---:|---:|
| full | 0.5758 | 0.2620 | 0.1623 |
| downweight current-reference FPR | 0.5820 | 0.2753 | 0.1427 |
| drop current high/missing-FPR lineages | 0.5361 | 0.2854 | 0.1785 |
| drop assigned high-FPR reference | 0.5349 | 0.2777 | 0.1873 |

结论：Tal1 的 fate-deviation/off-normal enrichment 不是由 high-FPR lineage 单独驱动。Doublet 作为 target reference 时必须降权解释；Endothelium、NMP、Neural crest、Doublet 的 lineage-level claims 不应作为强主结论。

## Balanced cell-count sensitivity

做了两种 downsampling：

- Equal total cells：每个 cohort 抽到 15,817 cells。
- Matched observed-lineage counts：只保留 Control/Tal1/T 共有 observed lineages，每个 lineage 每个 cohort 最多抽 200 cells，200 次重复。

关键结果：

| scenario | cohort | within | fate | abnormal |
|---|---|---:|---:|---:|
| equal total cells | Control | 0.9547 | 0.0339 | 0.0114 |
| equal total cells | T | 0.9434 | 0.0426 | 0.0140 |
| equal total cells | Tal1 | 0.5759 | 0.2620 | 0.1621 |
| matched observed-lineage counts | Control | 0.9346 | 0.0433 | 0.0221 |
| matched observed-lineage counts | T | 0.9288 | 0.0578 | 0.0134 |
| matched observed-lineage counts | Tal1 | 0.6215 | 0.2096 | 0.1689 |

Lineage-balanced 后 Tal1 fate-deviation fraction 下降，但 fate+abnormal 仍约 0.3785，远高于 T/control。因此 Tal1 信号不是单纯 cell-count 或 lineage-composition imbalance。

## Tal1 fate-deviation 指向哪些 reference lineages

Tal1 fate-deviation cells 共 7,415 个。主要 target reference lineages：

| target reference lineage | cells | fraction of Tal1 fate-deviation | calibration note |
|---|---:|---:|---|
| Stripped | 2,987 | 0.4028 | FPR 0.0732，非 high-FPR，但生物学解释弱 |
| Mesenchyme | 2,026 | 0.2732 | FPR 0.0052，较可靠 |
| Doublet | 759 | 0.1024 | high-FPR，需降权 |
| Spinal cord | 348 | 0.0469 | FPR 0.0285 |
| Pharyngeal mesoderm | 346 | 0.0467 | FPR 0.0423 |
| Rostral neurectoderm | 262 | 0.0353 | FPR 0.0463 |
| Haematoendothelial progenitors | 237 | 0.0320 | FPR 0.0365 |
| Caudal Mesoderm | 213 | 0.0287 | FPR 0.0404 |

主要 observed-to-target routes：

- Allantois -> Mesenchyme：873 cells，占 Allantois fate-deviation 的 0.8151。
- Forebrain/Midbrain/Hindbrain -> Stripped/Spinal cord/Rostral neurectoderm：神经相关 observed lineage 多数转向 Stripped 或 neural references。
- ExE mesoderm -> Mesenchyme/Stripped/Doublet：其中 Doublet route 需降权。
- Paraxial mesoderm -> Stripped：558 cells，占该 observed lineage fate-deviation 的 0.6186。

解释上应把 Tal1 fate-deviation 写成“多方向偏离”：一支指向 mesenchymal/extraembryonic-mesoderm-like 状态，一支指向 poorly resolved/stripped states，一小支指向 neural/spinal references。Doublet 和 Stripped 不应被过度生物学化。

## Marker/module 对 Tal1 异常的解释

全部 marker genes 在 Tal1 和 control H5AD 中均可用。相对于 Control within-stage-normal baseline：

- Tal1 abnormal/off-normal 最高的正向 module shifts 是 cardiac mesoderm (+0.1440)、allantois/ExE mesoderm (+0.1404)、endoderm (+0.0974)。
- Tal1 fate-deviation 最高的正向 shifts 是 allantois/ExE mesoderm (+0.0936)、WNT response (+0.0875)、paraxial/somite (+0.0201)。
- Tal1 within-stage-normal cells 仍保留 haematoendothelial (+0.0467) 和 endothelial (+0.0463) module signal。
- Tal1 所有三类的 erythroid module 均接近 0.01，而 control baseline erythroid module 很高。这与 Tal1/Scl 对血液发生的已知必需性一致。

文献约束：

- Pijuan-Sala et al. 2019 的原始 mouse gastrulation atlas 使用单细胞数据展示 Tal1-/- chimeric embryos 发生 early mesoderm diversification defects；该 paper 是当前 Tal1 chimera 数据的原始出处之一：[Nature 2019](https://www.nature.com/articles/s41586-019-0933-9)。
- MouseGastrulationData 文档说明 chimera 设计是把 fluorescent ESCs 注入 wild-type E3.5 embryos；host cells 可补偿整体发育，KO injected cells 可被捕获并用于分析 aberrant behaviour：[MouseGastrulationData vignette](https://bioconductor.org/packages/release/data/experiment/vignettes/MouseGastrulationData/inst/doc/MouseGastrulationData.html)。
- Tal1ChimeraData 文档确认 Tal1 processed data 共 4 samples，并包含 tomato、stage.mapped、celltype.mapped 等 metadata：[Tal1ChimeraData](https://rdrr.io/github/MarioniLab/MouseGastrulationData/man/Tal1ChimeraData.html)。
- TChimeraData 文档确认 T/Brachyury chimera 使用 paired embryo-pool design，默认保留 14 个 QC-passing samples：[TChimeraData](https://rdrr.io/github/MarioniLab/MouseGastrulationData/man/TChimeraData.html)。
- 经典 SCL/Tal1 knockout/chimera 研究显示 SCL/Tal1 对所有 hematopoietic lineages 的发育必不可少：[Porcher et al. 1996](https://www.sciencedirect.com/science/article/pii/S0092867400800768)。

因此，DevGuard 的 Tal1 结果应写成：Tal1 loss 在 E8.5 chimeric injected cells 中造成强烈 normality loss，主要体现为 erythroid/haematopoietic output failure 叠加 mesodermal fate allocation instability；而 T/Brachyury chimera 在当前 E8.5 tomato-positive subset 中整体接近 control normality，只表现出小幅 fate/off-normal enrichment。

## GSE212050 organoid heterogeneity control

GSE212050 的 GEO 摘要强调 individual gastruloids 存在 mesodermal 或 neural bias，并指出 organoid perturbation study 如果不建模 inter-organoid heterogeneity 会有高 false-positive/false-negative 风险：[GSE212050 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212050)。

本轮新增：

- GSE212050 strict sample heldout control classification。
- Heldout false-positive by organoid/sample。
- Fixed-embedding leave-one-organoid-out conformal FPR。
- Organoid lineage-bias table。

结果：

- Heldout sample-level false-positive fraction：mean 0.0511，max 0.0933。
- Leave-one-organoid-out FPR：mean 0.0504，median 0.0435，max 0.1825。
- 12/124 leave-one-organoid tests 的 FPR > 0.10。

解释：总体 calibration 接近 alpha=0.05，但少数 organoid/sample unit 仍显示明显高 FPR。这支持主文中继续保留 GSE212050 作为 heterogeneity stress control，而不是把单个 organoid 的偏差误写成 perturbation effect。

## GSE123187 spatial/tomo lineage mapping

GSE123187 是 gastruloid scRNA-seq + tomo-seq spatial transcriptomics 数据，原研究用于比较 gastruloid 和 mouse embryo 的 cell types 与 spatial expression patterns：[GSE123187 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123187)，[Nature 2020](https://www.nature.com/articles/s41586-020-2024-3)。

本轮新增：

- 从 RAW tar 选择 3 个 tomo-seq `coutb.tsv.gz` 文件，构建 `GSE123187_tomo_3files.h5ad`。
- 为每个 tomo segment 添加 ordinal `tomo_axis_fraction` 和 5-bin `tomo_axis_bin`。
- 用同一套 marker modules 对 tomo segments 做 top-module mapping。

结果：

- H5AD：1,146 segments，37,581 genes，3 tomo sources。
- Axis bins 分布均衡：每 bin 约 227-231 segments。
- Source-level top modules 主要为 haematoendothelial 与 paraxial/somite；axis-level table 已输出到 `tomo_axis_lineage_mapping.csv`。

限制：当前 mapping 是 marker-module preview，不是正式 cell-type transfer；tomo axis 方向尚未定向为 anterior/posterior。因此它可以支持 spatial/tomo ingestion 和粗 lineage-map preview，不能单独作为 fate-deviation validation 的强证据。

## 当前可写入论文的强弱分层

强证据：

- Tal1 比 T 和 heldout control 有大幅 fate-deviation/off-normal enrichment。
- 该 enrichment 经 sample bootstrap、high-FPR filtering/downweighting、equal-count downsampling、lineage-balanced downsampling 后仍保留。
- Tal1 marker/module 结果与 Tal1/Scl 血液发生和 mesoderm diversification 文献一致。

中等证据：

- Tal1 fate-deviation 主要指向 Mesenchyme、Stripped、Doublet、Spinal cord、Pharyngeal mesoderm 等 reference states；其中 Mesenchyme 最稳健，Stripped/Doublet 需要解释降权。
- GSE212050 leave-one-organoid-out 证明自然 organoid heterogeneity 会产生局部高 FPR，是必要 negative/stress control。

Preview 证据：

- GSE123187 tomo-seq 已完成 axis annotation 和 marker-module lineage mapping，但还需要正式 orientation、cell-type transfer 或 anchor-based mapping 才能作为 fate-deviation spatial validation。
