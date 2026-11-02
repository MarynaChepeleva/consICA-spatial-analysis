# consICA decompositions for the xenograft section of the ST paper
#
# counts_mat -> log2(x+1) -> consICA -> getGO against both genomes -> saveReport.
# Writes heap/cica_ntry20/ and output/cica_ntry20/. A run whose .rds already
# exists is skipped, so the script can be restarted.

library(dplyr)
library(Seurat)
library(consICA)

# Slide 6 is the reference, 5, 7 and 8 the replicates, 1 the tumour-free
# control. Slides 2, 3 and 4 also received no cells but are too shallow to
# decompose.
HAVE_SLIDE1_HM <- FALSE   # have human alignment for slide 1
SLIDES <- if (HAVE_SLIDE1_HM) c(6, 1, 5, 7, 8) else c(6, 5, 7, 8)

K_VALUES <- 12
# K_VALUES <- c(12, 20)

GO_SLIDES <- SLIDES
K20_SLIDES <- c(6, 5, 7, 8)

# 20 to match the rest of the paper.
NTRY <- 20
NCORES <- 1
SEED <- 2026
RUN_GO <- TRUE

# getGO on several cores dies serialising the object to the PSOCK workers.
GO_NCORES <- 1

cica_dir <- "heap/cica_ntry20"
rep_dir <- "output/cica_ntry20"
dir.create(cica_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rep_dir, showWarnings = FALSE, recursive = TRUE)

# Only the annotation columns are needed, and the object holds over 5 GB.
merged_seurat <- readRDS("heap/merged_seurat.rds")
meta <- merged_seurat@meta.data
rm(merged_seurat); gc()

counts_mat <- readRDS("heap/counts_mat.rds")

for (ncomp in K_VALUES) {
  for (slide_id in SLIDES) {

    if (ncomp == 20 && !slide_id %in% K20_SLIDES) next

    slide_spots <- rownames(meta)[meta$sample_id == slide_id]
    if (!length(slide_spots)) {
      message("slide ", slide_id, ": not in the merged object, skipped")
      next
    }
    counts_slide <- counts_mat[, slide_spots]

    pref <- paste0(cica_dir, "/cica", ncomp, "_slide_", slide_id, "_",
                   nrow(counts_slide), "x", ncol(counts_slide))

    # Each artefact is checked on its own, so a run that died during getGO
    # resumes there instead of skipping the slide.
    want_go <- RUN_GO && slide_id %in% GO_SLIDES
    have_cica <- file.exists(paste0(pref, ".rds"))
    have_hs <- file.exists(paste0(pref, "_GO.rds"))
    have_mm <- file.exists(paste0(pref, "_GO_mm.rds"))

    if (have_cica && (!want_go || (have_hs && have_mm))) {
      message("k=", ncomp, " slide ", slide_id, ": already done, skipped")
      next
    }

    message("k=", ncomp, " slide ", slide_id, ": ", ncol(counts_slide), " spots")

    if (have_cica) {
      message("  decomposition present, resuming at getGO")
      cica <- readRDS(paste0(pref, ".rds"))
    } else {
      mat_slide_log <- log2(counts_slide + 1) %>% as.matrix()
      gc()

      set.seed(SEED)
      cica <- consICA(mat_slide_log, ncomp = ncomp, ntry = NTRY,
                      ncores = NCORES, verbose = TRUE)

      saveRDS(cica, paste0(pref, ".rds"))
      rm(mat_slide_log); gc()
    }

    if (want_go) {

      # from the metadata frame: subsetting the Seurat object touches its images
      vars <- intersect(c("spot_type", "SR_hm_Graph", "SR_mm_6", "SR_hm_6"),
                        colnames(meta))
      var <- meta[colnames(cica$M), vars, drop = FALSE]

      rep_pref <- paste0(rep_dir, "/Report_cica", ncomp, "_slide_", slide_id,
                         "_", nrow(counts_slide), "x", ncol(counts_slide))

      if (have_hs) {
        cica <- readRDS(paste0(pref, "_GO.rds"))
      } else {
        cica <- getGO(cica, db = "BP", genome = "org.Hs.eg.db", ncores = GO_NCORES)
        saveRDS(cica, paste0(pref, "_GO.rds"))
        saveReport(cica, Var = var, file = paste0(rep_pref, ".pdf"))
      }

      if (!have_mm) {
        cica_mm <- getGO(cica, db = "BP", genome = "org.Mm.eg.db", ncores = GO_NCORES)
        saveRDS(cica_mm, paste0(pref, "_GO_mm.rds"))
        saveReport(cica_mm, Var = var, file = paste0(rep_pref, "_mm.pdf"))
        rm(cica_mm)
      }
    }

    rm(counts_slide, cica)
    gc()
  }
}

message("done")
