# READ-ONLY. Measures intra_uniqueness on the 150 saved B decompositions and reports the
# A/B cost ratio per slide rather than pooled over the grid. Reads B's own rds (~31 MB
# each) and A's metrics CSVs; A's 43 GB runs__* files are not read, see the identity below.
#
# Two definitions of the same quantity:
#   legacy  = 1 - mean(abs(cor(S)) with zero diagonal, averaged over all k^2 cells)
#   offdiag = 1 - mean(abs(cor(S))[upper.tri])           <- mean over distinct pairs
# Since the zeroed diagonal contributes k zeros to a k^2 denominator,
#   legacy = 1 - ((k-1)/k) * mean_offdiag   =>   offdiag = 1 - (k/(k-1)) * (1 - legacy)
# an exact transform. A's offdiag is derived from its recorded legacy values instead of
# recomputing from saved decompositions. The identity is verified numerically on B below.

OUT_DIR <- "output/benchmark_update_2026-09-23"
RUN_DIR <- "output/ica_perturbations"
K <- 20L

# ---- metrics ----
iu_legacy <- function(S) { # compute_internal_uniqueness(), verbatim
  cm <- abs(cor(S)); cm[is.na(cm)] <- 0; diag(cm) <- 0
  1 - mean(cm)
}
iu_offdiag <- function(S) { # mean over distinct pairs
  cm <- abs(cor(S)); cm[is.na(cm)] <- 0
  1 - mean(cm[upper.tri(cm)])
}

# A. intra_uniqueness on B's 150 decompositions ----
files <- list.files(OUT_DIR, pattern = "^runs__.*consICA_ntry1__.*__n10.*\\.rds$", full.names = TRUE)
stopifnot(length(files) == 15L)

rows <- do.call(rbind, lapply(files, function(f) {
  obj <- readRDS(f); m <- obj$metadata
  out <- do.call(rbind, lapply(seq_along(obj$runs), function(i) {
    S <- obj$runs[[i]]$S_matrix
    data.frame(slide = m$slide_id, mode = m$mode, repeat_i = i, k = ncol(S),
               n_genes = nrow(S),
               iu_legacy  = iu_legacy(S),
               iu_offdiag = iu_offdiag(S),
               stringsAsFactors = FALSE)
  }))
  rm(obj); gc(); out
}))
write.csv(rows, file.path(OUT_DIR, "B_intra_uniqueness_per_run.csv"), row.names = FALSE)

# verify the algebraic identity on real data before using it on A
rows$offdiag_from_legacy <- 1 - (rows$k / (rows$k - 1)) * (1 - rows$iu_legacy)
id_err <- max(abs(rows$offdiag_from_legacy - rows$iu_offdiag))
cat(sprintf("Identity check  offdiag = 1 - (k/(k-1))(1 - legacy):  max abs error = %.3g\n", id_err))
cat(if (id_err < 1e-12) "  -> exact; A's offdiag can be derived from its recorded legacy values\n\n"
    else "  -> DOES NOT HOLD; do not derive A's offdiag, recompute it instead\n\n")

b_sum <- aggregate(cbind(iu_legacy, iu_offdiag) ~ slide + mode, rows,
                   function(v) c(mean = mean(v), sd = sd(v)))
b_sum <- do.call(data.frame, b_sum)
names(b_sum) <- sub("\\.", "_", names(b_sum))
cat("=== B (consICA, ntry = 1): intra_uniqueness over 10 repeats ===\n")
print(b_sum, row.names = FALSE, digits = 6)
cat(sprintf("\nB overall: legacy mean %.6f (min %.6f)  |  offdiag mean %.6f (min %.6f)\n",
            mean(rows$iu_legacy), min(rows$iu_legacy),
            mean(rows$iu_offdiag), min(rows$iu_offdiag)))
cat(sprintf("mean |cor| between distinct components, B: %.3g\n", 1 - mean(rows$iu_offdiag)))

# ---- A for reference: legacy from the metrics CSVs, offdiag by the identity ----
mf <- list.files(RUN_DIR, pattern = "^metrics_FINAL__.*__k20__n10__gp5pct\\.csv$", full.names = TRUE)
a <- do.call(rbind, lapply(mf, function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  d <- d[d$Method == "consICA", c("Method", "Mode", "Intra_Uniqueness")]
  d$slide <- sub("^metrics_FINAL__([^_]+)__.*$", "\\1", basename(f))
  d
}))
if (!is.null(a) && nrow(a)) {
  a$iu_legacy <- a$Intra_Uniqueness
  a$iu_offdiag <- 1 - (K / (K - 1)) * (1 - a$iu_legacy)
  a <- a[, c("slide", "Mode", "iu_legacy", "iu_offdiag")]
  cat("\n=== A (consICA, ntry = 20): legacy as recorded, offdiag derived ===\n")
  print(a, row.names = FALSE, digits = 6)
  write.csv(a, file.path(OUT_DIR, "A_intra_uniqueness_derived.csv"), row.names = FALSE)
}

# B. cost ratio A / B, PER SLIDE ----
man <- read.csv(file.path(OUT_DIR, "manifest_B.csv"), stringsAsFactors = FALSE)
b_cost <- aggregate(elapsed_sec ~ slide, man, sum)
names(b_cost)[2] <- "B_ica_sec"

a_cost <- do.call(rbind, lapply(mf, function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  d <- d[d$Method == "consICA", ]
  data.frame(slide = sub("^metrics_FINAL__([^_]+)__.*$", "\\1", basename(f)),
             A_ica_sec = sum(d$Total_Runtime_sec, na.rm = TRUE), stringsAsFactors = FALSE)
}))

cost <- merge(a_cost, b_cost, by = "slide")
cost$ratio_A_over_B <- cost$A_ica_sec / cost$B_ica_sec
cost <- cost[order(cost$ratio_A_over_B), ]
write.csv(cost, file.path(OUT_DIR, "AB_cost_by_slide.csv"), row.names = FALSE)

cat("\n=== Cost of 20-try consensus relative to a single run, per slide (ICA time only) ===\n")
print(cost, row.names = FALSE, digits = 4)
cat(sprintf("\npooled over slides: %.0f s / %.0f s = %.2fx\n",
            sum(cost$A_ica_sec), sum(cost$B_ica_sec),
            sum(cost$A_ica_sec) / sum(cost$B_ica_sec)))
cat("\nCAVEAT: A ran with ncores = 4 on 7 of its 15 cells (spot on TENX115/TENX117/MEND40/NCBI600,\n",
    "gene on TENX115/TENX117/NCBI600), so A's recorded time there is 4-core wall clock. The\n",
    "single-core ratio would be LARGER than shown for those slides. B was ncores = 1 throughout.\n", sep = "")
cat("\nWrote B_intra_uniqueness_per_run.csv, A_intra_uniqueness_derived.csv, AB_cost_by_slide.csv\n")
