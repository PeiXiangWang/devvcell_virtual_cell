#!/usr/bin/env Rscript

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
local_lib <- file.path(root, ".r-lib", paste0(R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1]))
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages({
  required <- c("Seurat", "SeuratDisk")
})

args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[[1]] else "data/external/GSE208369/GSE208369_KO_merged.RDS.gz"
output <- if (length(args) >= 2) args[[2]] else "data/external/response_transfer/primary_perturbation.h5ad"

missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing, collapse = ", "),
    "\nInstall them before conversion, then rerun:\n",
    "Rscript scripts/convert_seurat_rds_to_h5ad.R ",
    input,
    " ",
    output,
    call. = FALSE
  )
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

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
h5seurat <- sub("\\.h5ad$", ".h5seurat", output)

message("Reading RDS: ", input)
obj <- read_rds_auto(input)
message("Writing intermediate h5Seurat: ", h5seurat)
SeuratDisk::SaveH5Seurat(obj, filename = h5seurat, overwrite = TRUE)
message("Converting to H5AD: ", output)
SeuratDisk::Convert(h5seurat, dest = "h5ad", overwrite = TRUE)
message("Done: ", output)
