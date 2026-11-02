# Moran's I for control B (consICA, ntry = 1) on the unperturbed condition,
# with EXACTLY the procedure used for the other methods in bench_moran_by_method.R:
#   kNN graph on unit coordinates, row-standardised (nb2listw style = "W")
#   k = 8 on Xenium, k = 6 on Visium; permutation null = 500 shuffles, set.seed(2026)
#   representative decomposition = pure_stability, iteration 1
# load_coords() and the slide table come from slide_coords.R, which the A script sources too,
# so the coordinate logic is identical by construction.
#
# ONE ADDITION. For A, iteration 1 is already an average of 20 runs; for B it is a single
# stochastic run, i.e. one draw. Running the 500-permutation null on all 10 repeats would cost
# 10x, so the null is computed on iteration 1 only (matching A), while the I statistic alone is
# computed for all 10 repeats to show how typical that draw is. That part is nearly free.
#
# READ-ONLY with respect to everything that exists: A's moran_by_method__*.csv are read, never
# written. B's outputs go to output/benchmark_update_2026-09-23/ under new names.

suppressPackageStartupMessages({
  library(spdep)
  library(dplyr)
})

COORDS <- "scripts/benchmark/slide_coords.R"
A_SUMMARY <- "output/benchmark_moran/moran_by_method__summary__k20.csv"
B_DIR <- "output/benchmark_update_2026-09-23"
METHOD <- "consICA_ntry1"
K <- 20L
N_ITER <- 10L
ITER <- 1L
N_PERM <- 500L
PERM_THR <- 0.05

# ---- the same slide table and coordinate logic the A script uses ----
source(COORDS)
stopifnot(is.function(load_coords), exists("slides"))

per_comp <- list(); spread <- list()

for (sl in names(slides)) {
  f <- file.path(B_DIR, sprintf("runs__%s__%s__pure_stability__k%d__n%d.rds", sl, METHOD, K, N_ITER))
  if (!file.exists(f)) { message("MISSING: ", basename(f)); next }

  coord <- load_coords(slides[[sl]]$file)
  k_nn <- if (slides[[sl]]$platform == "Xenium") 8L else 6L
  runs <- readRDS(f)$runs

  M <- runs[[ITER]]$M_matrix
  if (is.null(M)) stop("M_matrix absent in ", basename(f))
  spots <- intersect(colnames(M), rownames(coord))
  stopifnot(length(spots) > 0L)
  lw <- nb2listw(knn2nb(knearneigh(as.matrix(coord[spots, c("x", "y")]), k = k_nn)),
                 style = "W", zero.policy = TRUE)
  S0 <- Szero(lw); n0 <- length(spots)
  cat(sprintf("[%s] %s | k=%d | %d components x %d units\n",
              sl, slides[[sl]]$platform, k_nn, nrow(M), length(spots)))

  # representative decomposition: full procedure, identical to A
  for (ic in seq_len(nrow(M))) {
    v <- M[ic, spots]
    mt <- moran.test(v, lw, zero.policy = TRUE)
    obsI <- mt$estimate[[1]]
    set.seed(2026)
    null <- replicate(N_PERM, moran(sample(v), lw, n0, S0, zero.policy = TRUE)$I)
    per_comp[[length(per_comp) + 1L]] <- data.frame(
      slide = sl, platform = slides[[sl]]$platform, method = METHOD, k_nn = k_nn,
      component = paste0("IC", ic), Moran_I = obsI,
      analytic_p = mt$p.value, perm_p = (sum(null >= obsI) + 1) / (N_PERM + 1),
      stringsAsFactors = FALSE)
  }

  # how typical is that single draw: I only, all 10 repeats, no permutations
  for (it in seq_along(runs)) {
    Mi <- runs[[it]]$M_matrix
    sp <- intersect(colnames(Mi), rownames(coord))
    if (!identical(sp, spots)) next              # pure mode: the unit set is the same every repeat
    spread[[length(spread) + 1L]] <- data.frame(
      slide = sl, repeat_i = it,
      median_Moran_I = median(vapply(seq_len(nrow(Mi)),
                        function(ic) moran(Mi[ic, sp], lw, n0, S0, zero.policy = TRUE)$I,
                        numeric(1))),
      stringsAsFactors = FALSE)
  }
  rm(runs); gc()
}

pc <- do.call(rbind, per_comp)
write.csv(pc, file.path(B_DIR, "moran_B__per_component__k20.csv"), row.names = FALSE)

b_sum <- pc %>% group_by(slide, platform, method) %>%
  summarise(n_comp = n(), median_Moran_I = median(Moran_I), min_Moran_I = min(Moran_I),
            max_Moran_I = max(Moran_I), frac_perm_sig = mean(perm_p < PERM_THR), .groups = "drop")
write.csv(b_sum, file.path(B_DIR, "moran_B__summary__k20.csv"), row.names = FALSE)

cat("\n=== B (consICA, ntry = 1): Moran's I, representative decomposition ===\n")
print(as.data.frame(b_sum), row.names = FALSE, digits = 4)

spr <- do.call(rbind, spread)
if (!is.null(spr)) {
  sp_sum <- spr %>% group_by(slide) %>%
    summarise(n_repeats = n(), median_of_medians = median(median_Moran_I),
              min = min(median_Moran_I), max = max(median_Moran_I), .groups = "drop")
  write.csv(spr, file.path(B_DIR, "moran_B__across_repeats.csv"), row.names = FALSE)
  cat("\n=== how typical is iteration 1: per-repeat median Moran's I (B, 10 repeats) ===\n")
  print(as.data.frame(sp_sum), row.names = FALSE, digits = 4)
}

# ---- head-to-head against the existing A numbers (read only) ----
if (file.exists(A_SUMMARY)) {
  A <- read.csv(A_SUMMARY, stringsAsFactors = FALSE)
  A <- A[A$method == "consICA", c("slide", "median_Moran_I", "frac_perm_sig")]
  names(A)[2:3] <- c("A_median_I", "A_frac_sig")
  Bs <- as.data.frame(b_sum)[, c("slide", "median_Moran_I", "frac_perm_sig")]
  names(Bs)[2:3] <- c("B_median_I", "B_frac_sig")
  cmp <- merge(A, Bs, by = "slide")
  cmp$diff_I <- cmp$A_median_I - cmp$B_median_I
  cmp$favours <- ifelse(abs(cmp$diff_I) < 0.005, "~equal",
                        ifelse(cmp$diff_I > 0, "A (consensus)", "B (single run)"))
  write.csv(cmp, file.path(B_DIR, "AB_moran.csv"), row.names = FALSE)
  cat("\n=== A vs B: median Moran's I over the 20 components ===\n")
  print(cmp, row.names = FALSE, digits = 4)
  cat("\nHigher Moran's I means more spatial autocorrelation. It does not by itself",
      "\ndemonstrate a more correct biological signal.\n")
} else {
  cat("\nA summary not found at", A_SUMMARY, "- head-to-head skipped.\n")
}

cat("\nWrote moran_B__per_component__k20.csv, moran_B__summary__k20.csv,",
    "\n      moran_B__across_repeats.csv, AB_moran.csv  (all in", B_DIR, ")\n")
cat("A's moran_by_method__*.csv were read only, not modified.\n")
