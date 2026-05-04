# DevVCell 外部数据下载记录

版本：2026-05-05

本记录只追踪公开数据的获取状态、体量和校验信息。原始数据文件位于 `data/external/`，该目录不进入 Git。

## 已下载数据与元数据

| 数据集 | 文件 | 本地路径 | 大小 | SHA256 |
|---|---|---|---:|---|
| GSE208369 | `GSE208369_KO_merged.RDS.gz` | `data/external/GSE208369/GSE208369_KO_merged.RDS.gz` | 769,426,387 bytes | `DAFAD6C0567873FE449B876C010916409165A0E81AE164DE1CFE1A871E52907C` |
| GSE123187 | `GSE123187_README.txt` | `data/external/metadata/GSE123187_README.txt` | 8,044 bytes | `E7C414B7FF0EB9700AA06016720C4F86A2572B95B44C8E4EC1C7F5D91C639B4B` |
| GSE212050 | `GSE212050_sample_metadata_final.txt.gz` | `data/external/metadata/GSE212050_sample_metadata_final.txt.gz` | 18,965,533 bytes | `52816ED28DA09E57F5DF3E77ED9966DDBB00BC1A6F2878678A114D7B00B4A072` |

## 暂未下载的大文件

| 数据集 | 文件 | 体量 | 原因 |
|---|---|---:|---|
| GSE208369 | `GSE208369_RAW.tar` | 3.8 Gb | 正式复现用，优先级低于 processed RDS |
| GSE212050 | `GSE212050_seurat_final.rds.gz` | 1.7 Gb | 先下载 metadata，待 primary perturbation pipeline 跑通后下载 |
| GSE123187 | `GSE123187_RAW.tar` | 351.2 Mb | 先下载 README 解析格式，后续用于 spatial/tomo validation |

## 下一步

1. 用 R/Seurat 读取 `GSE208369_KO_merged.RDS.gz` 并导出 H5AD。
2. 输出到 `data/external/response_transfer/primary_perturbation.h5ad`。
3. 运行 `scripts/export_external_perturbation_response.py`。
4. 若需要完全复现，再下载 `GSE208369_RAW.tar`。

## 转换环境检查

已新增转换脚本：

```bash
Rscript scripts/convert_seurat_rds_to_h5ad.R data/external/GSE208369/GSE208369_KO_merged.RDS.gz data/external/response_transfer/primary_perturbation.h5ad
```

2026-05-05 更新：本机 R 版本为 4.4.3；已将 `Seurat`、`SeuratObject`、`hdf5r`、`SeuratDisk` 安装到项目本地 `.r-lib/4.4`，并完成 GSE208369 RDS 转 H5AD。由于 `SeuratDisk` 与 `SeuratObject 5.x` 的 `slot=` 调用不兼容，实际采用两步转换：

```bash
Rscript scripts/export_seurat_rds_components.R data/external/GSE208369/GSE208369_KO_merged.RDS.gz data/external/response_transfer/GSE208369_components
python scripts/build_h5ad_from_seurat_components.py --components data/external/response_transfer/GSE208369_components --output data/external/response_transfer/primary_perturbation.h5ad
```

转换后数据：

```text
data/external/response_transfer/primary_perturbation.h5ad
20000 cells x 33577 genes
```
