# The slide table and the coordinate loader shared by the two Moran scripts,
# bench_moran_by_method.R and bench_B_moran.R. Definitions only: sourcing this file
# runs nothing. platform is what those scripts read to set the neighbourhood size.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(zellkonverter)
})

HEST_DIR <- Sys.getenv("HEST_ROOT", "data/HEST")

slides <- list(
  TENX115 = list(file = file.path(HEST_DIR, "skin/st/TENX115.h5ad"), platform = "Xenium"),
  TENX117 = list(file = file.path(HEST_DIR, "skin/st/TENX117.h5ad"), platform = "Xenium"),
  MEND40 = list(file = file.path(HEST_DIR, "skin/st/MEND40.h5ad"), platform = "Visium"),
  NCBI776 = list(file = file.path(HEST_DIR, "breast/st/NCBI776.h5ad"), platform = "Visium"),
  NCBI600 = list(file = file.path(HEST_DIR, "bladder/st/NCBI600.h5ad"), platform = "Visium"))

load_coords <- function(h5ad) {
  sce <- readH5AD(h5ad)
  cd <- as.data.frame(colData(sce))
  if ("spatial" %in% reducedDimNames(sce)) {
    xy <- as.data.frame(reducedDim(sce, "spatial"))[, 1:2]
  } else {
    xc <- intersect(c("pxl_col_in_fullres", "array_col", "imagecol"), colnames(cd))[1]
    yc <- intersect(c("pxl_row_in_fullres", "array_row", "imagerow"), colnames(cd))[1]
    xy <- cd[, c(xc, yc)]
  }
  colnames(xy) <- c("x", "y")
  rownames(xy) <- colnames(sce)
  xy
}
