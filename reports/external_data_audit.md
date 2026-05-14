# External Data Audit

- selected_dataset_id: E1_mouse_gastrulation_wt_chimera_sample1
- source: local matrix.mtx/obs.csv/var.csv components derived from public MouseGastrulationData.
- cells_after_filtering: 1800
- genes_after_hvg: 2000
- time/stage field: `stage.mapped` -> `time_numeric` / `time_point`
- lineage field: `celltype.mapped` -> `lineage` / `cell_type`
- lineage barcode: unavailable; no lineage-validated claim is made.
- embedding fit: external data only; no internal teacher or internal embedding was used.

## Stage Counts

| time_point   |   count |
|:-------------|--------:|
| E7.0         |     161 |
| E7.5         |     453 |
| E7.25        |      87 |
| E7.75        |     618 |
| E8.0         |     388 |
| E8.25        |      93 |
