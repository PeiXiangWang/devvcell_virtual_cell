# DevVCell 外部公开数据集计划

版本：2026-05-05

## 1. 数据策略

正式版 DevVCell 必须纳入至少一个真实 perturbation scRNA-seq 数据集。数据不提交 Git；只提交：

- accession registry；
- 下载地址；
- 转换脚本；
- checksum；
- 处理后的 metadata；
- 分析结果摘要。

原始文件放在：

```text
data/external/
```

该目录已被 `.gitignore` 排除。

## 2. 第一优先级：GSE208369

GEO：

```text
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208369
```

用途：

- primary perturbation-response dictionary；
- RA gastruloid response；
- BMP inhibition response，LDN-treated；
- TF knockout response，TBX6-KO/PAX3-KO；
- human gastruloid 到 mouse embryo atlas 的 cross-species response transfer。

GEO 页面显示该数据公开于 2024-06-20，包含 conventional/RA gastruloid、LDN-treated、NTC、TBX6-KO 和 PAX3-KO samples。处理文件包括：

```text
GSE208369_KO_merged.RDS.gz 733.8 Mb
GSE208369_RAW.tar 3.8 Gb
```

正式处理步骤：

1. 下载 RDS 或 RAW；
2. 转换为 H5AD；
3. 统一 obs 字段：`condition`、`perturbation`、`cell_type`、`batch`、`time`、`dose`；
4. human-mouse ortholog mapping；
5. 输出 `data/external/response_transfer/primary_perturbation.h5ad`；
6. 运行 `export_external_perturbation_response.py`。

## 3. 第二优先级：GSE212050

GEO：

```text
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212050
```

用途：

- mouse gastruloid time course；
- individual organoid heterogeneity；
- donor/organoid covariate negative control；
- 检查 response-recovery 是否被 organoid composition bias 误导。

处理文件：

```text
GSE212050_seurat_final.rds.gz 1.7 Gb
GSE212050_sample_metadata_final.txt.gz 18.1 Mb
```

正式处理重点：

- 保留 individual gastruloid barcode；
- 建立 day 3-5 latent trajectory；
- 用 WNT-biased vs neural-biased organoids 作为 heterogeneity-aware validation。

## 4. 第三优先级：GSE123187

GEO：

```text
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123187
```

用途：

- mouse gastruloid scRNA-seq；
- spatial/tomo-seq validation；
- mesoderm/neural/somite response transfer validation；
- niche proxy 的空间合理性检查。

处理文件：

```text
GSE123187_RAW.tar 351.2 Mb
GSE123187_README.txt
```

正式处理重点：

- 解析 tomo-seq 和 scRNA-seq TSV；
- 将 tomo/spatial axis 映射为 niche validation；
- 检查 Tbx4-high niche 和 Glis3-high target 是否存在空间或轴向邻近证据。

## 5. 正常发育流形外部验证

### 5.1 MouseGastrulationData / E-MTAB-6967

用途：

- early embryo normal manifold validation；
- E6.5-E8.5 trajectory；
- chimaera/knockout-related resource。

入口：

```text
https://bioconductor.org/packages/release/data/experiment/html/MouseGastrulationData.html
```

### 5.2 MOCA

用途：

- organogenesis-stage trajectory validation；
- E9.5-E13.5；
- 2 million scale external atlas。

入口：

```text
https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/downloads
```

MOCA 很大，应在 primary perturbation pipeline 稳定后再接入。

## 6. PGC/gonad niche validation

### 6.1 GSE136441

用途：

- prenatal/neonatal mouse gonad/ovary atlas；
- PGC and gonadal somatic niche validation；
- 支撑 Glis3-high germ-cell state 和 Tbx4-high niche 的阶段共现。

### 6.2 GSE288206

用途：

- murine fetal gonad single-nucleus multiomics；
- PGC/supporting lineage regulatory validation；
- 如果 gene coverage 和 annotation 合适，用于 Tbx4-Glis3 niche case 的补充证据。

## 7. 下载原则

1. 优先下载 `GSE208369_KO_merged.RDS.gz`，因为体量相对可控且直接含 TF KO。
2. 大于 1GB 的文件下载前确认磁盘空间，不进入 Git。
3. 所有 RDS 转 H5AD 的转换步骤必须记录版本、包版本和 checksum。
4. 每个外部数据集必须生成一个 `*_metadata.csv`，保留 condition、perturbation、sample、batch、organism、cell type。
5. 跨物种数据必须明确 ortholog mapping 版本。

## 8. 当前仓库状态

数据 registry 已写入：

```text
config/external_datasets.json
```

如果外部 H5AD 尚未下载，先运行：

```bash
python scripts/export_external_perturbation_response.py --allow-missing
```

它会生成公开数据获取 manifest，不会伪造 response dictionary。

也可以用统一下载器先抓取小体量元数据：

```bash
python scripts/download_public_datasets.py --metadata-only --max-mb 50
```

下载器会写入 `data/external/download_manifest.csv`，包含本地路径、字节数、SHA256 和 gzip 校验状态。大文件需要显式提高 `--max-mb` 或不使用 `--metadata-only`。
