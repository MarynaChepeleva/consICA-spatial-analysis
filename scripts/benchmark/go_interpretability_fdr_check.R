# READ-ONLY diagnostic: a Bonferroni bound on the §3.3 interpretability counts, from the
# saved go_benchmark__*.rds. Those tables hold the nominal weight01 Fisher p, uncorrected,
# and keep only the top 20 terms per component, so an exact BH cannot be reconstructed from
# them. The best p per component is exact, and Bonferroni is stricter than BH, so:
#   Frac at p < 0.05/m  (m = assumed number of BP terms tested)  <=  Frac under BH.
# Frac_bonf == 1 therefore settles the question; Frac_bonf < 1 does not.
# go_interpretability_bh.R computes the exact BH by re-running topGO.
#
# Reads the saved rds and writes go_interpretability_fdr_check.csv.

go_dir <- "output/go_enrichment/"
out_dir <- "output/benchmark_summary/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

N_BONF <- c(5000, 15000)          # conservative BP-term counts for the Bonferroni bound
files <- list.files(go_dir, pattern = "^go_benchmark__.*__k20\\.rds$", full.names = TRUE)


best_p_per_component <- function(go_results) {   # min nominal Fisher p per component
  vapply(go_results, function(res) {
    if (is.null(res) || !nrow(res)) return(NA_real_)
    min(res$Fisher, na.rm = TRUE)
  }, numeric(1))
}

rows <- list()
for (f in files) {
  slide <- sub("^go_benchmark__(.*)__k20\\.rds$", "\\1", basename(f))
  obj <- readRDS(f)
  pm <- obj$topN$per_method               # primary selection (top_n_genes)
  for (mth in names(pm)) {
    bp <- best_p_per_component(pm[[mth]]$go_results)
    n <- length(bp)
    row <- data.frame(
      Slide = slide, Method = mth, n_comp = n,
      Frac_nominal  = mean(bp < 0.05, na.rm = TRUE),          # should match go_summary
      Frac_bonf5000 = mean(bp < 0.05 / N_BONF[1], na.rm = TRUE),
      Frac_bonf15000= mean(bp < 0.05 / N_BONF[2], na.rm = TRUE),
      worst_best_p  = max(bp, na.rm = TRUE),                  # least-significant component
      stringsAsFactors = FALSE)
    rows[[length(rows) + 1L]] <- row
  }
  rm(obj); gc()
}

res <- do.call(rbind, rows)
res <- res[order(res$Slide, res$Method), ]
write.csv(res, file.path(out_dir, "go_interpretability_fdr_check.csv"), row.names = FALSE)

cat("=== interpretability under a conservative correction (READ-ONLY) ===\n")
cat("Frac_nominal = weight01 p < 0.05, uncorrected.\n")
cat("Frac_bonf*   = Bonferroni bound, stricter than BH: == 1 settles the question.\n\n")
print(res, digits = 4, row.names = FALSE)
cat("\nWritten:", file.path(out_dir, "go_interpretability_fdr_check.csv"), "\n")
