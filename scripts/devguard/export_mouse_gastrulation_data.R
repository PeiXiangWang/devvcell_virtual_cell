#!/usr/bin/env Rscript

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
local_lib <- file.path(root, ".r-lib", paste0(R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1]))
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1) args[[1]] else "embryo_atlas"
out_dir <- if (length(args) >= 2) args[[2]] else file.path("data", "external", "MouseGastrulationData", dataset)
sample_arg <- if (length(args) >= 3) args[[3]] else ""

required <- c("MouseGastrulationData", "SingleCellExperiment", "SummarizedExperiment", "Matrix")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(MouseGastrulationData)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
})

parse_samples <- function(x) {
  if (!nzchar(x) || x == "all") {
    return(NULL)
  }
  as.integer(strsplit(x, ",")[[1]])
}

samples <- parse_samples(sample_arg)
message("Exporting MouseGastrulationData dataset: ", dataset)
sce <- switch(
  dataset,
  embryo_atlas = if (is.null(samples)) EmbryoAtlasData() else EmbryoAtlasData(samples = samples),
  wt_chimera = if (is.null(samples)) WTChimeraData() else WTChimeraData(samples = samples),
  tal1_chimera = if (is.null(samples)) Tal1ChimeraData() else Tal1ChimeraData(samples = samples),
  t_chimera = if (is.null(samples)) TChimeraData() else TChimeraData(samples = samples),
  stop("Unknown dataset: ", dataset, call. = FALSE)
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
mat <- SingleCellExperiment::counts(sce)
if (!inherits(mat, "dgCMatrix")) {
  mat <- methods::as(mat, "dgCMatrix")
}
obs <- as.data.frame(SummarizedExperiment::colData(sce))
obs$cell_barcode <- colnames(sce)
obs$mouse_gastrulation_dataset <- dataset
var <- as.data.frame(SummarizedExperiment::rowData(sce))
var$gene_id <- rownames(sce)
if (!"SYMBOL" %in% colnames(var)) {
  var$SYMBOL <- rownames(sce)
}

Matrix::writeMM(mat, file.path(out_dir, "matrix.mtx"))
write.csv(obs, file.path(out_dir, "obs.csv"), row.names = FALSE, quote = TRUE)
write.csv(var, file.path(out_dir, "var.csv"), row.names = FALSE, quote = TRUE)
manifest <- data.frame(
  dataset = dataset,
  samples = ifelse(is.null(samples), "all", paste(samples, collapse = ",")),
  n_genes = nrow(mat),
  n_cells = ncol(mat),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(out_dir, "component_manifest.csv"), row.names = FALSE, quote = TRUE)
message("Wrote components to: ", out_dir)
