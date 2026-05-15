# L2 Leakage Audit

- L2 preprocessing uses only the downloaded Biddy_2018_Nature h5ad.
- No internal teacher, internal PCA, internal labels or E1 labels are used to fit L2 embeddings.
- Clone identifiers are copied from the source `barcode_all` field; no clone labels are invented.
- Processed h5ad files are written under ignored data/processed paths and are not staged for Git.
