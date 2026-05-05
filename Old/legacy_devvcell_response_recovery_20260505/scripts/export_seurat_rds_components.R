#!/usr/bin/env Rscript

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
local_lib <- file.path(root, ".r-lib", paste0(R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1]))
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages({
  required <- c("Seurat", "SeuratObject", "Matrix")
})

args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[[1]] else "data/external/GSE208369/GSE208369_KO_merged.RDS.gz"
out_dir <- if (length(args) >= 2) args[[2]] else "data/external/response_transfer/GSE208369_components"

missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

read_rds_auto <- function(path) {
  attempts <- list(
    function() readRDS(path),
    function() {
      con <- gzfile(path, open = "rb")
      on.exit(close(con), add = TRUE)
      readRDS(con)
    },
    function() {
      con <- gzcon(gzfile(path, open = "rb"))
      on.exit(close(con), add = TRUE)
      readRDS(con)
    }
  )
  errors <- character()
  for (attempt in attempts) {
    result <- tryCatch(attempt(), error = function(err) err)
    if (!inherits(result, "error")) {
      return(result)
    }
    errors <- c(errors, conditionMessage(result))
  }
  stop("Unable to read RDS. Attempts failed with: ", paste(errors, collapse = " | "), call. = FALSE)
}

extract_assay_matrix <- function(obj, assay) {
  for (layer in c("counts", "data", "scale.data")) {
    matrix <- tryCatch(
      SeuratObject::GetAssayData(obj, assay = assay, layer = layer),
      error = function(err) err
    )
    if (!inherits(matrix, "error") && nrow(matrix) > 0 && ncol(matrix) > 0) {
      return(list(matrix = matrix, layer = layer))
    }
  }

  assay_obj <- obj[[assay]]
  for (slot_name in c("counts", "data", "scale.data")) {
    matrix <- tryCatch(slot(assay_obj, slot_name), error = function(err) err)
    if (!inherits(matrix, "error") && nrow(matrix) > 0 && ncol(matrix) > 0) {
      return(list(matrix = matrix, layer = slot_name))
    }
  }
  stop("Could not extract a non-empty matrix from assay: ", assay, call. = FALSE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading RDS: ", input)
obj <- read_rds_auto(input)
message("Object class: ", paste(class(obj), collapse = ", "))
assays <- SeuratObject::Assays(obj)
message("Assays: ", paste(assays, collapse = ", "))
assay <- if ("RNA" %in% assays) "RNA" else SeuratObject::DefaultAssay(obj)
if (is.null(assay) || !assay %in% assays) {
  assay <- assays[[1]]
}
message("Using assay: ", assay)

extracted <- extract_assay_matrix(obj, assay)
mat <- extracted$matrix
if (!inherits(mat, "dgCMatrix")) {
  mat <- methods::as(mat, "dgCMatrix")
}
message("Using layer: ", extracted$layer)
message("Matrix dimensions genes x cells: ", paste(dim(mat), collapse = " x "))

obs <- obj@meta.data
obs$cell_barcode <- rownames(obs)
var <- data.frame(gene_id = rownames(mat), gene_name = rownames(mat), stringsAsFactors = FALSE)

Matrix::writeMM(mat, file.path(out_dir, "matrix.mtx"))
write.csv(obs, file.path(out_dir, "obs.csv"), row.names = FALSE, quote = TRUE)
write.csv(var, file.path(out_dir, "var.csv"), row.names = FALSE, quote = TRUE)

manifest <- data.frame(
  input = input,
  assay = assay,
  layer = extracted$layer,
  n_genes = nrow(mat),
  n_cells = ncol(mat),
  n_obs_columns = ncol(obs),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(out_dir, "component_manifest.csv"), row.names = FALSE, quote = TRUE)
message("Wrote components to: ", out_dir)
