# Do IC12 and IC4 separate inside the tumour, or only outside it?
#
# Over the whole reference section IC4 accounts for 80.5% of the deviance and
# IC12 for 36.6%. This splits the section at the 47 annotated tumour units. The
# figures the manuscript quotes, 5% for IC12 and 48% for IC4, are the
# rest-of-section rows. Writes tumour_split.csv.

library(dplyr)
library(mgcv)

cica_dir <- "heap/cica_ntry20"
f_counts <- "heap/counts_mat.rds"
f_seurat <- "heap/merged_seurat.rds"
f_species <- "heap/gene_species.rds"
f_tumour <- "data/tumor_spots.csv"
out_dir <- "output/ST_paper"

ncomp <- 12
ref_slide <- 6
n_boot <- 2000
seed <- 2026
gam_k <- 5 # basis dimension for the fit on 47 units

# the spatial part needs the Seurat object, which is the slow read of the script
DO_SPATIAL <- TRUE

melanocytic <- c("PMEL", "TYR", "TYRP1", "MLANA", "SOX10", "TRPM1", "DCT")


# ---- data ----
roles <- read.csv(file.path(out_dir, "components.csv"), stringsAsFactors = FALSE)
ic_prog <- roles$ic[roles$role == "programme"]
ic_abun <- roles$ic[roles$role == "abundance"]
cat("programme", ic_prog, " abundance", ic_abun, "\n")

f_cica <- list.files(cica_dir, full.names = TRUE,
                     pattern = sprintf("^cica%d_slide_%d_[0-9]+x[0-9]+\\.rds$",
                                       ncomp, ref_slide))
M <- readRDS(f_cica[1])$M
spots <- colnames(M)

gene_species <- readRDS(f_species)
human_genes <- names(gene_species)[gene_species == "human"]

counts <- readRDS(f_counts)[, spots]
depth <- Matrix::colSums(counts)
human_frac <- Matrix::colSums(counts[intersect(human_genes, rownames(counts)), ]) /
  depth

tumour_barcodes <- as.character(read.csv(f_tumour)[[1]])
is_tumour <- spots %in% tumour_barcodes
cat("units", length(spots), " annotated tumour", sum(is_tumour), "\n")


# ---- correlation of a component with the human proportion, by subset ----
# The sign of a component is arbitrary, so only the size of a correlation and
# the agreement between subsets mean anything. The Spearman interval is a
# percentile bootstrap, since cor.test gives none.

set.seed(seed)

boot_rho <- function(x, y) {
  vals <- replicate(n_boot, {
    i <- sample(length(x), replace = TRUE)
    suppressWarnings(cor(x[i], y[i], method = "spearman"))
  })
  quantile(vals, c(0.025, 0.975), na.rm = TRUE)
}

cor_row <- function(ic, subset_name, keep) {
  w <- M[ic, spots][keep]
  f <- human_frac[keep]
  pe <- cor.test(w, f)
  sp <- suppressWarnings(cor.test(w, f, method = "spearman"))
  ci <- boot_rho(w, f)
  data.frame(component = ic, subset = subset_name, n = sum(keep),
             pearson_r = unname(pe$estimate),
             pearson_lo = pe$conf.int[1], pearson_hi = pe$conf.int[2],
             pearson_p = pe$p.value,
             spearman_rho = unname(sp$estimate),
             spearman_lo = unname(ci[1]), spearman_hi = unname(ci[2]),
             spearman_p = sp$p.value,
             row.names = NULL, stringsAsFactors = FALSE)
}

subsets <- list("whole section" = rep(TRUE, length(spots)),
                "annotated tumour" = is_tumour,
                "rest of the section" = !is_tumour)

cors <- lapply(c(ic_prog, ic_abun), function(ic) {
  lapply(names(subsets), function(s) cor_row(ic, s, subsets[[s]])) %>% bind_rows()
}) %>% bind_rows()

cat("\n=== correlation with the human proportion ===\n")
print(cors %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))) %>%
        as.data.frame(),
      row.names = FALSE)


# ---- how much of a component the human proportion accounts for ----
# The deviance of a smooth is not comparable between 3441 units and 47 at one
# basis dimension, since ten basis functions on 47 points take a share of the
# variance by their degrees of freedom alone. The linear R squared compares
# across the two.

fit_row <- function(ic, subset_name, keep) {
  w <- M[ic, spots][keep]
  f <- human_frac[keep]
  g_default <- gam(w ~ s(f, bs = "cs"), method = "REML")
  g_small <- gam(w ~ s(f, bs = "cs", k = gam_k), method = "REML")
  data.frame(component = ic, subset = subset_name, n = sum(keep),
             dev_gam_default = summary(g_default)$dev.expl,
             edf_default = sum(summary(g_default)$edf),
             dev_gam_k5 = summary(g_small)$dev.expl,
             edf_k5 = sum(summary(g_small)$edf),
             r2_linear = cor(w, f)^2,
             row.names = NULL, stringsAsFactors = FALSE)
}

fits <- lapply(c(ic_prog, ic_abun), function(ic) {
  lapply(names(subsets), function(s) fit_row(ic, s, subsets[[s]])) %>% bind_rows()
}) %>% bind_rows()

cat("\n=== deviance explained by the human proportion ===\n")
print(fits %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))) %>%
        as.data.frame(),
      row.names = FALSE)

write.csv(bind_rows(cors %>% mutate(block = "correlation"),
                    fits %>% mutate(block = "fit")),
          file.path(out_dir, "tumour_split.csv"), row.names = FALSE)
cat("\nwritten to", file.path(out_dir, "tumour_split.csv"), "\n")


# ---- what is left of the programme component once the proportion is removed ----
# Residuals of the smooth fitted on the tumour units alone. Spatial structure or
# a link to the melanocytic set means the component holds something the human
# proportion does not.

if (DO_SPATIAL) {
  w <- M[ic_prog, spots][is_tumour]
  f <- human_frac[is_tumour]
  res <- residuals(gam(w ~ s(f, bs = "cs", k = gam_k), method = "REML"))

  merged_seurat <- suppressWarnings(
    suppressMessages(Seurat::UpdateSeuratObject(readRDS(f_seurat))))
  xy <- Seurat::GetTissueCoordinates(merged_seurat, image = "slice1.6")
  cols <- if (all(c("imagecol", "imagerow") %in% colnames(xy))) {
    c("imagecol", "imagerow")
  } else {
    c("y", "x")
  }
  xy <- as.matrix(xy[spots[is_tumour], cols])
  rm(merged_seurat)
  gc()

  # the hexagonal Visium lattice has six neighbours, as in the benchmark scripts
  nb <- spdep::knn2nb(spdep::knearneigh(xy, k = 6))
  lw <- spdep::nb2listw(nb, style = "W")
  set.seed(seed)
  mc <- spdep::moran.mc(res, lw, nsim = 999)
  cat(sprintf("\n=== residuals of %s on the tumour units ===\nMoran's I %.3f, permutation p %.4f, %d units, k = 6\n",
              ic_prog, unname(mc$statistic), mc$p.value, length(res)))

  # markers as a share of the unit's counts, so depth does not drive the
  # association; the residual is already adjusted for the human proportion
  present <- intersect(melanocytic, rownames(counts))
  marker_cor <- lapply(present, function(g) {
    share <- as.vector(counts[g, is_tumour]) / depth[is_tumour]
    if (sd(share) == 0) return(NULL)
    ct <- suppressWarnings(cor.test(res, share, method = "spearman"))
    data.frame(gene = g, detected_in = sum(share > 0),
               spearman_rho = unname(ct$estimate), p = ct$p.value,
               row.names = NULL, stringsAsFactors = FALSE)
  }) %>% bind_rows()

  cat("\nresiduals against the melanocytic set, on the tumour units:\n")
  print(marker_cor %>%
          mutate(across(where(is.numeric), ~ signif(.x, 3))) %>%
          as.data.frame(),
        row.names = FALSE)
}
