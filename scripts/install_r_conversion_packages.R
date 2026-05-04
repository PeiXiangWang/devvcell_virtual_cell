#!/usr/bin/env Rscript

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
local_lib <- file.path(root, ".r-lib", paste0(R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1]))
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))

options(repos = c(CRAN = "https://cloud.r-project.org"))

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing R package: ", pkg)
    install.packages(pkg, lib = local_lib, dependencies = TRUE, type = "binary")
  } else {
    message("R package already available: ", pkg)
  }
}

for (pkg in c("remotes", "SeuratObject", "Seurat", "hdf5r", "Matrix")) {
  install_if_missing(pkg)
}

if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  message("Installing R package from GitHub: mojaveazure/seurat-disk")
  remotes::install_github(
    "mojaveazure/seurat-disk",
    lib = local_lib,
    upgrade = "never",
    dependencies = TRUE
  )
} else {
  message("R package already available: SeuratDisk")
}

status <- vapply(
  c("remotes", "SeuratObject", "Seurat", "hdf5r", "SeuratDisk"),
  requireNamespace,
  logical(1),
  quietly = TRUE
)
print(status)
if (!all(status)) {
  stop("Some required R conversion packages are still unavailable.", call. = FALSE)
}
