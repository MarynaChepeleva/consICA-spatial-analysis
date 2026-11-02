# Moran's I by method for the consensus-step comparison (Results).
# Reads saved pure_stability runs (K=20), takes M_matrix (iteration 1 = one
# representative decomposition per method), and computes per-component Moran's I
# with the SAME parameters as the demo scripts, so numbers are comparable to
# §3.1 / §3.2:
#   kNN graph on unit coordinates, row-standardised (nb2listw style="W")
#   k = 8 on Xenium, k = 6 on Visium (hex lattice)
#   permutation null = 500 shuffles of the component weights; analytic p from moran.test
#
# The slide table and load_coords() are in slide_coords.R, shared with bench_B_moran.R.
#
# pure_stability runs carry no __gp5pct suffix in their filenames.
# Only consica_object is stripped in the benchmark; S_matrix and M_matrix are always
# saved, so Moran is computable on every slide, NCBI776 included
# (store_consica_objects = FALSE). The stop() below guards that.

suppressPackageStartupMessages({
  library(spdep)
  library(dplyr)
})

source("scripts/benchmark/slide_coords.R")

runs_dir <- "output/ica_perturbations/"
out_dir <- "output/benchmark_moran/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

methods <- c("consICA", "CFS_As_Is", "icafast_no_Rmat", "icafast_with_Rmat")
K <- 20L
N_ITER <- 10L
ITER <- 1L # representative decomposition (deterministic methods: all 10 identical)
N_PERM <- 500L
PERM_THR <- 0.05 # threshold for "fraction of components with perm_p below it"

rows <- list()

for (sl in names(slides)) {
  coord <- load_coords(slides[[sl]]$file)
  k_nn <- if (slides[[sl]]$platform == "Xenium") 8L else 6L
  cat(sprintf("\n[%s] platform %s | k=%d\n", sl, slides[[sl]]$platform, k_nn))

  for (mth in methods) {
    f <- file.path(runs_dir, sprintf("runs__%s__%s__pure_stability__k%d__n%d.rds", sl, mth, K, N_ITER))
    if (!file.exists(f)) { message("  MISSING file, skipped: ", basename(f)); next }

    M <- readRDS(f)$runs[[ITER]]$M_matrix
    if (is.null(M)) stop("M_matrix ABSENT in ", basename(f), " (iter ", ITER, ") -- cannot compute Moran")

    spots <- intersect(colnames(M), rownames(coord))
    stopifnot(length(spots) > 0L)
    lw <- nb2listw(knn2nb(knearneigh(as.matrix(coord[spots, c("x", "y")]), k = k_nn)),
                   style = "W", zero.policy = TRUE)
    S0 <- Szero(lw); n0 <- length(spots)
    cat(sprintf("  %-18s %d components x %d units\n", mth, nrow(M), length(spots)))

    for (ic in seq_len(nrow(M))) {
      v <- M[ic, spots]
      mt <- moran.test(v, lw, zero.policy = TRUE)
      obsI <- mt$estimate[[1]]
      set.seed(2026)
      null <- replicate(N_PERM, moran(sample(v), lw, n0, S0, zero.policy = TRUE)$I)
      rows[[length(rows) + 1L]] <- data.frame(
        slide = sl,
        platform = slides[[sl]]$platform,
        method = mth,
        k_nn = k_nn,
        component = paste0("IC", ic),
        Moran_I = obsI,
        analytic_p = mt$p.value,
        perm_p = (sum(null >= obsI) + 1) / (N_PERM + 1),
        stringsAsFactors = FALSE)
    }
  }
}

per_comp <- do.call(rbind, rows)
write.csv(per_comp, file.path(out_dir, "moran_by_method__per_component__k20.csv"), row.names = FALSE)

summary_tbl <- per_comp %>%
  group_by(slide, platform, method) %>%
  summarise(n_comp = n(),
            median_Moran_I = median(Moran_I),
            min_Moran_I = min(Moran_I),
            max_Moran_I = max(Moran_I),
            frac_perm_sig = mean(perm_p < PERM_THR),
            .groups = "drop")
write.csv(summary_tbl, file.path(out_dir, "moran_by_method__summary__k20.csv"), row.names = FALSE)

cat("\n=== Moran by method (summary) ===\n")
print(as.data.frame(summary_tbl))
cat("\nWritten:\n ",
    file.path(out_dir, "moran_by_method__per_component__k20.csv"), "\n ",
    file.path(out_dir, "moran_by_method__summary__k20.csv"), "\n")
