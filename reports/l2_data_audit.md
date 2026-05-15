# L2 Data Audit

- selected_dataset: Biddy_2018_Nature / GSE99915 CellTag reprogramming
- full_loaded_shape: 18803 x 2000
- downsampled_shape: 3200 x 2000
- time_points: ['Day12', 'Day15', 'Day21', 'Day28', 'Day6', 'Day9']
- clone_id_field: barcode_all
- usable_clones_size_ge_5: 536
- usable_clones_size_ge_20: 123

Counts by time and lineage:

| time_point   | lineage      |   n_cells |
|:-------------|:-------------|----------:|
| Day6         | Failed       |         0 |
| Day6         | Others       |        82 |
| Day6         | Reprogrammed |         5 |
| Day9         | Failed       |         0 |
| Day9         | Others       |       893 |
| Day9         | Reprogrammed |        10 |
| Day12        | Failed       |         0 |
| Day12        | Others       |      2981 |
| Day12        | Reprogrammed |        49 |
| Day15        | Failed       |         0 |
| Day15        | Others       |      3249 |
| Day15        | Reprogrammed |        39 |
| Day21        | Failed       |        10 |
| Day21        | Others       |      4079 |
| Day21        | Reprogrammed |       839 |
| Day28        | Failed       |       591 |
| Day28        | Others       |      4101 |
| Day28        | Reprogrammed |      1875 |
