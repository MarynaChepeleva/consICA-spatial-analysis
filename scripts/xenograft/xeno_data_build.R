# Builds the three objects the other xenograft scripts read, from the public
# Visium deposit unpacked one folder per sample. Each sample was aligned twice,
# against GRCh38 and against mm10, and the two matrices are stacked. Seventeen
# symbols occur in both and the human copy takes an hg- prefix.
#
# The gene set comes from xeno_genes.csv, not from the filter. A slide with no
# human alignment is skipped, which on the public deposit is slide 1, so the
# control section cannot be rebuilt from it.
#
# tumor_spots.csv holds the units annotated as tumour, based on the Space Ranger
# clustering and confirmed by the histology annotation. The clustering is outside
# this deposit.
#
# Writes heap/counts_mat.rds, heap/merged_seurat.rds and heap/gene_species.rds.

library(Matrix)
library(Seurat)

XENO_DIR <- "data/zenodo/MBrM_Visium_FF"   # <- set to the unpacked archive
f_genes <- "data/xeno_genes.csv"   # genome per gene, written outside this deposit
out_dir <- "heap"

SLIDES <- c(1, 5, 6, 7, 8)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gene_table <- read.csv(f_genes, stringsAsFactors = FALSE)
gene_set <- gene_table$gene
cat("gene set:", length(gene_set), "genes,",
    sum(gene_table$genome == "human"), "human\n")

# the images keep the names the figure scripts ask for, slice1 for the first
# slide and the slide number after it
image_name <- function(slide) if (slide == 1) "slice1" else paste0("slice1.", slide)


# ---- one slide ----
# A gene the deposit does not carry enters as a zero row, so every slide ends up
# on the same gene set in the same order and the slides can be merged.

read_slide <- function(slide) {
  d <- file.path(XENO_DIR, sprintf("%02d", slide))
  f_hm <- file.path(d, "filtered_feature_bc_matrix_human.h5")
  f_mm <- file.path(d, "filtered_feature_bc_matrix_mouse.h5")
  if (!file.exists(f_hm)) {
    message("slide ", slide, ": the deposit carries no human alignment, skipped")
    return(NULL)
  }

  hm <- Read10X_h5(f_hm)
  mm <- Read10X_h5(f_mm)
  dups <- intersect(rownames(hm), rownames(mm))
  rownames(hm)[rownames(hm) %in% dups] <-
    paste0("hg-", rownames(hm)[rownames(hm) %in% dups])

  common <- intersect(colnames(hm), colnames(mm))
  x <- rbind(hm[, common, drop = FALSE], mm[, common, drop = FALSE])

  absent <- setdiff(gene_set, rownames(x))
  if (length(absent)) {
    x <- rbind(x, Matrix(0, nrow = length(absent), ncol = ncol(x), sparse = TRUE,
                         dimnames = list(absent, colnames(x))))
  }
  x <- x[gene_set, , drop = FALSE]

  obj <- CreateSeuratObject(counts = x, assay = "Spatial")

  # the deposit carries the high-resolution stain and no low-resolution one, so
  # the scale factor that goes with the stored image is the high-resolution one
  img <- Read10X_Image(d, image.name = "tissue_hires_image.png",
                       filter.matrix = FALSE)
  img@scale.factors$lowres <- img@scale.factors$hires
  img <- img[colnames(obj), ]
  DefaultAssay(img) <- "Spatial"
  obj[[image_name(slide)]] <- img

  obj$sample_id <- slide
  cat(sprintf("slide %d: %d genes x %d spots, %d symbols prefixed, %d zero rows\n",
              slide, nrow(obj), ncol(obj), length(dups), length(absent)))
  RenameCells(obj, add.cell.id = paste0("S", slide))
}


# ---- every slide, merged ----

objs <- lapply(SLIDES, read_slide)
names(objs) <- SLIDES
objs <- objs[!vapply(objs, is.null, logical(1))]
stopifnot(length(objs) > 1)

merged_seurat <- merge(objs[[1]], y = objs[-1])
cat("\nmerged:", nrow(merged_seurat), "genes x", ncol(merged_seurat), "spots,",
    "images", paste(Images(merged_seurat), collapse = " "), "\n")

counts_mat <- SeuratObject::LayerData(
  JoinLayers(merged_seurat[["Spatial"]]), layer = "counts")
counts_mat <- counts_mat[gene_set, colnames(merged_seurat), drop = FALSE]

gene_species <- setNames(gene_table$genome, gene_table$gene)

saveRDS(counts_mat, file.path(out_dir, "counts_mat.rds"))
saveRDS(merged_seurat, file.path(out_dir, "merged_seurat.rds"))
saveRDS(gene_species, file.path(out_dir, "gene_species.rds"))

cat("\nwritten to", out_dir, "\n")
cat("counts_mat", paste(dim(counts_mat), collapse = " x "),
    " spots per slide:",
    paste(table(merged_seurat$sample_id), collapse = " "), "\n")
