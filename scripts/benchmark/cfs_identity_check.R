# Checks the statement in cfs_ICASpatial.R on generated data: the copy and CellsFromSpace's
# own ICASpatial return the same loadings, and the seed the copy exposes changes nothing.
# Reads no data and writes nothing. Run from the project root; needs CFS, Seurat, ica and
# purrr installed.

suppressPackageStartupMessages({ 
  library(CFS)
  library(Seurat)
  library(ica)
  library(purrr) })
source("scripts/benchmark/cfs_ICASpatial.R")

# a small object carrying the one slot both functions read. @scale.data is a v3 assay slot,
# as SCTransform's own SCTAssay is, so the object is built as v3.
options(Seurat.object.assay.version = "v3")
set.seed(2026)
ng <- 200
ns <- 120
k <- 5
counts <- matrix(rpois(ng * ns, 5), ng, ns,
                 dimnames = list(paste0("g", seq_len(ng)), paste0("c", seq_len(ns))))
obj <- CreateSeuratObject(counts, assay = "SCT")
obj <- SetAssayData(obj, assay = "SCT", layer = "scale.data",
                    new.data = t(scale(t(as.matrix(counts)))))

a <- CFS:::ICASpatial(obj, nics = k, maxit = 600) # upstream, seed fixed at 42
b <- custom_ICASpatial(obj, nics = k, maxit = 600, seed.use = 1)
d <- custom_ICASpatial(obj, nics = k, maxit = 600, seed.use = 2)

L <- function(obj1) obj1@reductions$ica@feature.loadings
cat("upstream vs copy      identical:", identical(L(a), L(b)),
    " max |difference|:", max(abs(L(a) - L(b))), "\n")
cat("copy at seed 1 vs 2  identical:", identical(L(b), L(d)),
    " max |difference|:", max(abs(L(b) - L(d))), "\n")
cat("icafast default Rmat  :", deparse(formals(ica::icafast)$Rmat), "\n")
