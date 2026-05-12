# SwarmLineage-OT Data Audit

This audit was generated from local files only. It records fields that are available for modelling; it does not certify biological validity.

## Summary

- Files scanned: 84
- Readable expression objects: 25
- Time/stage-capable files: 21
- Perturbation/condition-capable files: 23
- Spatial candidate files: 2
- Ligand-receptor candidate files: 11

## File-Level Schema

### `data\external\download_manifest.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\GSE208369\GSE208369_KO_merged.RDS.gz`

- format: gz
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\GSE212050\components\component_manifest.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\GSE212050\components\matrix.mtx`

- format: mtx
- readable: True
- shape: 30406 cells/rows × 77683 genes/features
- time/day/stage/pseudotime: timepoint, timepoint.base, timepoint.demultiplexed, stage.mapped.original, cellstage.score.original, stage.mapped.meso, cellstage.score.meso, stage.mapped.extended, cellstage.score.extended, stage.mapped.descendant.extended, cellstage.score.descendant, velocity_pseudotime, day
- cell type/cluster/fate/lineage: seurat_clusters, cluster, celltype.mapped.original, celltype.score.original, celltype.mapped.original.base, celltype.mapped.original.meso, celltype.score.meso, celltype.mapped.extended, celltype.score.extended, celltype.mapped.descendant.extended, celltype.score.descendant
- lineage/barcode: barcode, cell_barcode
- batch/donor/condition/treatment/perturbation: batch, sample, estimated_growth_rates
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\GSE212050\components\obs.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\GSE212050\components\var.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\GSE212050\components_downsample_15000\component_manifest.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\GSE212050\components_downsample_15000\matrix.mtx`

- format: mtx
- readable: True
- shape: 30406 cells/rows × 15000 genes/features
- time/day/stage/pseudotime: timepoint, timepoint.base, timepoint.demultiplexed, stage.mapped.original, cellstage.score.original, stage.mapped.meso, cellstage.score.meso, stage.mapped.extended, cellstage.score.extended, stage.mapped.descendant.extended, cellstage.score.descendant, velocity_pseudotime, day
- cell type/cluster/fate/lineage: seurat_clusters, cluster, celltype.mapped.original, celltype.score.original, celltype.mapped.original.base, celltype.mapped.original.meso, celltype.score.meso, celltype.mapped.extended, celltype.score.extended, celltype.mapped.descendant.extended, celltype.score.descendant
- lineage/barcode: barcode, cell_barcode
- batch/donor/condition/treatment/perturbation: batch, sample, estimated_growth_rates
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\GSE212050\components_downsample_15000\obs.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\GSE212050\components_downsample_15000\var.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\GSE212050\components_strict_sample\component_manifest.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\GSE212050\components_strict_sample\matrix.mtx`

- format: mtx
- readable: True
- shape: 30406 cells/rows × 13285 genes/features
- time/day/stage/pseudotime: timepoint, timepoint.base, timepoint.demultiplexed, stage.mapped.original, cellstage.score.original, stage.mapped.meso, cellstage.score.meso, stage.mapped.extended, cellstage.score.extended, stage.mapped.descendant.extended, cellstage.score.descendant, velocity_pseudotime, day
- cell type/cluster/fate/lineage: seurat_clusters, cluster, celltype.mapped.original, celltype.score.original, celltype.mapped.original.base, celltype.mapped.original.meso, celltype.score.meso, celltype.mapped.extended, celltype.score.extended, celltype.mapped.descendant.extended, celltype.score.descendant
- lineage/barcode: barcode, cell_barcode
- batch/donor/condition/treatment/perturbation: batch, sample, estimated_growth_rates
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\GSE212050\components_strict_sample\obs.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\GSE212050\components_strict_sample\var.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\GSE212050\GSE212050_feature_metadata_final.txt.gz`

- format: gz
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\GSE212050\GSE212050_sample_metadata_final.txt.gz`

- format: gz
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\GSE212050\GSE212050_seurat_final.rds.gz`

- format: gz
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\metadata\GSE212050_sample_metadata_final.txt.gz`

- format: gz
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\embryo_atlas_sample1\component_manifest.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\embryo_atlas_sample1\matrix.mtx`

- format: mtx
- readable: True
- shape: 29452 cells/rows × 360 genes/features
- time/day/stage/pseudotime: stage, theiler, cluster.stage, cluster.theiler
- cell type/cluster/fate/lineage: cluster, cluster.sub, cluster.stage, cluster.theiler, celltype
- lineage/barcode: barcode, cell_barcode
- batch/donor/condition/treatment/perturbation: sample, sequencing.batch
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\embryo_atlas_sample1\obs.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\embryo_atlas_sample1\var.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\t_chimera_full\component_manifest.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\t_chimera_full\matrix.mtx`

- format: mtx
- readable: True
- shape: 29453 cells/rows × 36931 genes/features
- time/day/stage/pseudotime: stage, stage.mapped, somite.subct.mapped
- cell type/cluster/fate/lineage: celltype.mapped
- lineage/barcode: barcode, cell_barcode
- batch/donor/condition/treatment/perturbation: sample
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\t_chimera_full\obs.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\t_chimera_full\var.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\tal1_chimera_full\component_manifest.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\tal1_chimera_full\matrix.mtx`

- format: mtx
- readable: True
- shape: 29453 cells/rows × 56122 genes/features
- time/day/stage/pseudotime: stage, stage.mapped
- cell type/cluster/fate/lineage: celltype.mapped
- lineage/barcode: barcode, cell_barcode
- batch/donor/condition/treatment/perturbation: sample
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\tal1_chimera_full\obs.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\tal1_chimera_full\var.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\tal1_chimera_sample1\component_manifest.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\tal1_chimera_sample1\matrix.mtx`

- format: mtx
- readable: True
- shape: 29453 cells/rows × 14751 genes/features
- time/day/stage/pseudotime: stage, stage.mapped
- cell type/cluster/fate/lineage: celltype.mapped
- lineage/barcode: barcode, cell_barcode
- batch/donor/condition/treatment/perturbation: sample
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\tal1_chimera_sample1\obs.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\tal1_chimera_sample1\var.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\wt_chimera_full\component_manifest.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\wt_chimera_full\matrix.mtx`

- format: mtx
- readable: True
- shape: 29453 cells/rows × 30703 genes/features
- time/day/stage/pseudotime: stage, stage.mapped
- cell type/cluster/fate/lineage: celltype.mapped
- lineage/barcode: barcode, cell_barcode
- batch/donor/condition/treatment/perturbation: sample
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\wt_chimera_full\obs.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\wt_chimera_full\var.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\wt_chimera_sample1\component_manifest.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\wt_chimera_sample1\matrix.mtx`

- format: mtx
- readable: True
- shape: 29453 cells/rows × 2882 genes/features
- time/day/stage/pseudotime: stage, stage.mapped
- cell type/cluster/fate/lineage: celltype.mapped
- lineage/barcode: barcode, cell_barcode
- batch/donor/condition/treatment/perturbation: sample
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\wt_chimera_sample1\obs.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\MouseGastrulationData\wt_chimera_sample1\var.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\response_transfer\GSE208369_components\component_manifest.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\response_transfer\GSE208369_components\matrix.mtx`

- format: mtx
- readable: True
- shape: 33577 cells/rows × 20000 genes/features
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: seurat_clusters, orig.clusters, cell_type
- lineage/barcode: cell_barcode
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\response_transfer\GSE208369_components\obs.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\response_transfer\GSE208369_components\var.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\external\response_transfer\primary_perturbation.h5ad`

- format: h5ad
- readable: True
- shape: 20000 cells/rows × 33577 genes/features
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: seurat_clusters, orig.clusters, cell_type
- lineage/barcode: cell_barcode
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: PCNA, MCM5, TYMS, Pcna, Mcm5, Tyms
- ligand-receptor genes: Fgf8, Fgfr1, Wnt3, Fzd1, Kit, Cxcl12, Cxcr4, Vegfa, Kdr, Dll1, Notch1, FGF8, FGFR1, WNT3, FZD1, KITLG, KIT, CXCL12, CXCR4, VEGFA
- true lineage tracing likely present: True
- perturb-seq/treatment labels likely present: False

### `data\external\scperturb\AdamsonWeissman2016_GSM2406675_10X001.h5ad`

- format: h5ad
- readable: True
- shape: 5768 cells/rows × 35635 genes/features
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: celltype
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: perturbation, disease, perturbation_type
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: MKI67, TOP2A, PCNA, MCM5, TYMS, UBE2C, Mki67, Top2a, Pcna, Mcm5, Tyms, Ube2c
- ligand-receptor genes: Fgf8, Fgfr1, Wnt3, Fzd1, Kit, Cxcl12, Cxcr4, Vegfa, Kdr, Dll1, Notch1, FGF8, FGFR1, WNT3, FZD1, KITLG, KIT, CXCL12, CXCR4, VEGFA
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: True

### `data\external\scperturb\DatlingerBock2021.h5ad`

- format: h5ad
- readable: True
- shape: 39194 cells/rows × 25904 genes/features
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: celltype
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: sample, perturbation, perturbation_2, disease, perturbation_type, perturbation_type_2
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: MKI67, TOP2A, PCNA, MCM5, TYMS, UBE2C, Mki67, Top2a, Pcna, Mcm5, Tyms, Ube2c
- ligand-receptor genes: Fgf8, Fgfr1, Wnt3, Fzd1, Kit, Cxcl12, Cxcr4, Vegfa, Kdr, Dll1, Notch1, FGF8, FGFR1, WNT3, FZD1, KITLG, KIT, CXCL12, CXCR4, VEGFA
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: True

### `data\processed\cell_level_subset_v1.h5ad`

- format: h5ad
- readable: True
- shape: 19156 cells/rows × 3000 genes/features
- time/day/stage/pseudotime: author_day, author_somite_count, development_stage, stage_num
- cell type/cluster/fate/lineage: author_major_cell_cluster, author_cell_type, cell_type
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: donor_id, disease
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: Kitl, Kit, Kdr, KIT, KDR
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\processed\cell_level_subset_v1.manifest.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\processed\devguard\devguard_quick_mouse.h5ad`

- format: h5ad
- readable: True
- shape: 1470 cells/rows × 90 genes/features
- time/day/stage/pseudotime: time_point, time_numeric
- cell type/cluster/fate/lineage: cell_type, lineage
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: condition, perturbation_name, perturbation_type, sample_id, batch, is_perturbed
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: True

### `data\processed\devguard\devguard_stress_mouse.h5ad`

- format: h5ad
- readable: True
- shape: 2178 cells/rows × 120 genes/features
- time/day/stage/pseudotime: time_point, time_numeric, stress_fixture_truth_time
- cell type/cluster/fate/lineage: cell_type, lineage, stress_fixture_truth_lineage
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: condition, perturbation_name, perturbation_type, sample_id, batch, is_perturbed
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: True

### `data\processed\devguard\GSE123187_preview_4files.h5ad`

- format: h5ad
- readable: True
- shape: 1536 cells/rows × 24086 genes/features
- time/day/stage/pseudotime: time_point, time_numeric
- cell type/cluster/fate/lineage: cell_type, lineage
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: condition, perturbation_name, perturbation_type, sample_id, batch, is_perturbed
- spatial coordinates: tomo_position
- multimodal fields: not detected
- cell-cycle/proliferation markers: MKI67, TOP2A, PCNA, MCM5, TYMS, UBE2C, Mki67, Top2a, Pcna, Mcm5, Tyms, Ube2c
- ligand-receptor genes: Fgf8, Fgfr1, Wnt3, Fzd1, Kitl, Kit, Cxcl12, Cxcr4, Vegfa, Kdr, Dll1, Notch1, FGF8, FGFR1, WNT3, FZD1, KIT, CXCL12, CXCR4, VEGFA
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: True

### `data\processed\devguard\GSE123187_tomo_3files.h5ad`

- format: h5ad
- readable: True
- shape: 1146 cells/rows × 37581 genes/features
- time/day/stage/pseudotime: time_point, time_numeric
- cell type/cluster/fate/lineage: cell_type, lineage
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: condition, perturbation_name, perturbation_type, sample_id, batch, is_perturbed
- spatial coordinates: tomo_position, tomo_axis_fraction, tomo_axis_bin
- multimodal fields: not detected
- cell-cycle/proliferation markers: MKI67, TOP2A, PCNA, MCM5, TYMS, UBE2C, Mki67, Top2a, Pcna, Mcm5, Tyms, Ube2c
- ligand-receptor genes: Fgf8, Fgfr1, Wnt3, Fzd1, Kitl, Kit, Cxcl12, Cxcr4, Vegfa, Kdr, Dll1, Notch1, FGF8, FGFR1, WNT3, FZD1, KIT, CXCL12, CXCR4, VEGFA
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: True

### `data\processed\devguard\GSE212050_downsample_15000.h5ad`

- format: h5ad
- readable: True
- shape: 15000 cells/rows × 30406 genes/features
- time/day/stage/pseudotime: timepoint, timepoint.base, timepoint.demultiplexed, stage.mapped.original, cellstage.score.original, stage.mapped.meso, cellstage.score.meso, stage.mapped.extended, cellstage.score.extended, stage.mapped.descendant.extended, cellstage.score.descendant, velocity_pseudotime, day, time_point, time_numeric
- cell type/cluster/fate/lineage: seurat_clusters, cluster, celltype.mapped.original, celltype.score.original, celltype.mapped.original.base, celltype.mapped.original.meso, celltype.score.meso, celltype.mapped.extended, celltype.score.extended, celltype.mapped.descendant.extended, celltype.score.descendant, cell_type, lineage
- lineage/barcode: barcode, cell_barcode
- batch/donor/condition/treatment/perturbation: batch, sample, estimated_growth_rates, condition, perturbation_name, perturbation_type, sample_id, is_perturbed
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: True
- perturb-seq/treatment labels likely present: True

### `data\processed\devguard\GSE212050_strict_sample_13285.h5ad`

- format: h5ad
- readable: True
- shape: 13285 cells/rows × 30406 genes/features
- time/day/stage/pseudotime: timepoint, timepoint.base, timepoint.demultiplexed, stage.mapped.original, cellstage.score.original, stage.mapped.meso, cellstage.score.meso, stage.mapped.extended, cellstage.score.extended, stage.mapped.descendant.extended, cellstage.score.descendant, velocity_pseudotime, day, time_point, time_numeric
- cell type/cluster/fate/lineage: seurat_clusters, cluster, celltype.mapped.original, celltype.score.original, celltype.mapped.original.base, celltype.mapped.original.meso, celltype.score.meso, celltype.mapped.extended, celltype.score.extended, celltype.mapped.descendant.extended, celltype.score.descendant, cell_type, lineage
- lineage/barcode: barcode, cell_barcode
- batch/donor/condition/treatment/perturbation: batch, sample, estimated_growth_rates, condition, perturbation_name, perturbation_type, sample_id, is_perturbed
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: True
- perturb-seq/treatment labels likely present: True

### `data\processed\devguard\MouseGastrulationData_embryo_atlas_sample1.h5ad`

- format: h5ad
- readable: True
- shape: 360 cells/rows × 29452 genes/features
- time/day/stage/pseudotime: stage, theiler, cluster.stage, cluster.theiler, time_point, time_numeric
- cell type/cluster/fate/lineage: cluster, cluster.sub, cluster.stage, cluster.theiler, celltype, cell_type, lineage
- lineage/barcode: barcode, cell_barcode
- batch/donor/condition/treatment/perturbation: sample, sequencing.batch, condition, perturbation_name, perturbation_type, sample_id, batch, is_perturbed
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: MKI67, TOP2A, PCNA, MCM5, TYMS, UBE2C, Mki67, Top2a, Pcna, Mcm5, Tyms, Ube2c
- ligand-receptor genes: Fgf8, Fgfr1, Wnt3, Fzd1, Kitl, Kit, Cxcl12, Cxcr4, Vegfa, Kdr, Dll1, Notch1, FGF8, FGFR1, WNT3, FZD1, KIT, CXCL12, CXCR4, VEGFA
- true lineage tracing likely present: True
- perturb-seq/treatment labels likely present: True

### `data\processed\devguard\MouseGastrulationData_integrated_chimera_controls.h5ad`

- format: h5ad
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\processed\devguard\MouseGastrulationData_t_chimera_full.h5ad`

- format: h5ad
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\processed\devguard\MouseGastrulationData_tal1_chimera_full.h5ad`

- format: h5ad
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\processed\devguard\MouseGastrulationData_tal1_chimera_sample1.h5ad`

- format: h5ad
- readable: True
- shape: 14751 cells/rows × 29453 genes/features
- time/day/stage/pseudotime: stage, stage.mapped, time_point, time_numeric
- cell type/cluster/fate/lineage: celltype.mapped, cell_type, lineage
- lineage/barcode: barcode, cell_barcode
- batch/donor/condition/treatment/perturbation: sample, condition, perturbation_name, perturbation_type, sample_id, batch, is_perturbed
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: MKI67, TOP2A, PCNA, MCM5, TYMS, UBE2C, Mki67, Top2a, Pcna, Mcm5, Tyms, Ube2c
- ligand-receptor genes: Fgf8, Fgfr1, Wnt3, Fzd1, Kitl, Kit, Cxcl12, Cxcr4, Vegfa, Kdr, Dll1, Notch1, FGF8, FGFR1, WNT3, FZD1, KIT, CXCL12, CXCR4, VEGFA
- true lineage tracing likely present: True
- perturb-seq/treatment labels likely present: True

### `data\processed\devguard\MouseGastrulationData_wt_chimera_full.h5ad`

- format: h5ad
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\processed\devguard\MouseGastrulationData_wt_chimera_sample1.h5ad`

- format: h5ad
- readable: True
- shape: 2882 cells/rows × 29453 genes/features
- time/day/stage/pseudotime: stage, stage.mapped, time_point, time_numeric
- cell type/cluster/fate/lineage: celltype.mapped, cell_type, lineage
- lineage/barcode: barcode, cell_barcode
- batch/donor/condition/treatment/perturbation: sample, condition, perturbation_name, perturbation_type, sample_id, batch, is_perturbed
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: MKI67, TOP2A, PCNA, MCM5, TYMS, UBE2C, Mki67, Top2a, Pcna, Mcm5, Tyms, Ube2c
- ligand-receptor genes: Fgf8, Fgfr1, Wnt3, Fzd1, Kitl, Kit, Cxcl12, Cxcr4, Vegfa, Kdr, Dll1, Notch1, FGF8, FGFR1, WNT3, FZD1, KIT, CXCL12, CXCR4, VEGFA
- true lineage tracing likely present: True
- perturb-seq/treatment labels likely present: True

### `data\rdeg_neural_cell_mvp\edge_importance_scores.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\rdeg_neural_cell_mvp\grn_learned_network.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\rdeg_neural_cell_mvp\hub_tf_ranking.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\rdeg_neural_cell_mvp\next_step_pair_metrics.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\rdeg_neural_cell_mvp\nodes_summary.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\rdeg_neural_cell_mvp\ot_pair_metrics.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\rdeg_neural_cell_mvp\perturbation_scores.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\rdeg_neural_cell_mvp\rescue_experiments.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\rdeg_neural_cell_mvp\rollout_error_matrix.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\rdeg_neural_cell_mvp\stage_module_means.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\rdeg_neural_cell_mvp\system_specific_edges.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\rdeg_neural_cell_mvp\temporal_sensitivity.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\rdeg_neural_cell_mvp\tf_fate_mutual_info.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\rdeg_neural_cell_mvp\tf_knockout_results.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `data\scLine_pro.h5ad`

- format: h5ad
- readable: False
- time/day/stage/pseudotime: author_day, development_stage
- cell type/cluster/fate/lineage: author_major_cell_cluster, author_cell_type, cell_type
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: donor_id, author_experimental_id, sex
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `processed\ot_couplings\moscot_coupling_index.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `processed\ot_couplings\wot_coupling_index.csv`

- format: csv
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `processed\ot_fate_probabilities.parquet`

- format: parquet
- readable: False
- time/day/stage/pseudotime: not detected
- cell type/cluster/fate/lineage: not detected
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: not detected
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: not detected
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: False

### `processed\ot_teacher.h5ad`

- format: h5ad
- readable: True
- shape: 8000 cells/rows × 2000 genes/features
- time/day/stage/pseudotime: author_day, author_somite_count, development_stage, stage_num, time_numeric, time_point, ot_target_time
- cell type/cluster/fate/lineage: author_major_cell_cluster, author_cell_type, cell_type, lineage, fate_prob_neural, fate_prob_erythroid, fate_prob_mesoderm_muscle, ot_fate_max
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: donor_id, disease, condition, batch, is_perturbed
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: Kitl, Kit, Kdr, KIT, KDR
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: True

### `processed\swarmlineage_input.h5ad`

- format: h5ad
- readable: True
- shape: 8000 cells/rows × 2000 genes/features
- time/day/stage/pseudotime: author_day, author_somite_count, development_stage, stage_num, time_numeric, time_point
- cell type/cluster/fate/lineage: author_major_cell_cluster, author_cell_type, cell_type, lineage
- lineage/barcode: not detected
- batch/donor/condition/treatment/perturbation: donor_id, disease, condition, batch, is_perturbed
- spatial coordinates: not detected
- multimodal fields: not detected
- cell-cycle/proliferation markers: not detected
- ligand-receptor genes: Kitl, Kit, Kdr, KIT, KDR
- true lineage tracing likely present: False
- perturb-seq/treatment labels likely present: True

## Main Modelling Input Decision

`data/processed/cell_level_subset_v1.h5ad` is used as the default quick real-data input because it has eight ordered developmental stages, lineage/cell-type labels, and a manageable 19,156 × 3,000 matrix. The 12.6GB `data/scLine_pro.h5ad` is retained as a larger extension target.

If no real time field exists, stage or pseudotime is treated as an exploratory fallback. Main claims in generated manuscripts are gated accordingly.
