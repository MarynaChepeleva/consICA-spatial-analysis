# Assembles the head-to-head table, A against B.
#   A = consICA, ntry = 20 (consensus)
#   B = consICA, ntry = 1  (single run), computed by the same engine in bench_B_run.R
# Primary metric: Matched_Stability, a one-to-one Hungarian assignment on |cor| of gene
# signatures, pairwise between the 10 repeats within an arm. .pair_abs_cor intersects gene
# names before correlating, so gene dropout is handled by intersection.
#
# Reads both arms' metrics and joins them; new output goes under new names into
# output/benchmark_update_2026-09-23/.
#
# No Tissue column: tissue belongs to Table 1 of the manuscript, not to a generated file.
# No p-values: with 10 repeats the 45 pairwise comparisons inside an arm are dependent, so
# magnitude and spread are reported and significance is not claimed.

RUN_DIR <- "output/ica_perturbations"
OUT_DIR <- "output/benchmark_update_2026-09-23"
MODES <- c("pure_stability", "spot_perturbation", "gene_perturbation")

ARM <- c(consICA = "A (consensus, ntry=20)", consICA_ntry1 = "B (single run, ntry=1)")

# ---- B: $metrics from the 15 saved cells ----
bf <- list.files(OUT_DIR, pattern = "^runs__.*consICA_ntry1__.*__n10.*\\.rds$", full.names = TRUE)
stopifnot(length(bf) == 15L)
B <- do.call(rbind, lapply(bf, function(f) {
  o <- readRDS(f); m <- o$metrics; m$slide <- o$metadata$slide_id
  rm(o); gc(); m
}))

# ---- A: the k=20 / gp=5% metrics CSVs, consICA rows only ----
af <- list.files(RUN_DIR, pattern = "^metrics_FINAL__.*__k20__n10__gp5pct\\.csv$", full.names = TRUE)
A <- do.call(rbind, lapply(af, function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  d <- d[d$Method == "consICA", ]
  d$slide <- sub("^metrics_FINAL__([^_]+)__.*$", "\\1", basename(f))
  d
}))

keep <- c("slide", "Method", "Mode", "Mean_Stability", "SD_Stability",
          "Matched_Stability", "SD_Matched", "TopGene_Jaccard", "SD_TopGene",
          "N_Repro_50", "N_Repro_70", "Intra_Uniqueness")
long <- rbind(A[, intersect(keep, names(A))], B[, intersect(keep, names(B))])
long$arm <- unname(ARM[long$Method])
long$Mode <- factor(long$Mode, MODES)
long <- long[order(long$slide, long$Mode, long$Method), ]

# Intra_Uniqueness is carried for traceability only. Measured on B: exactly 1.000 on all 150
# decompositions (mean |cor| at machine precision), i.e. at the ceiling by construction
# for a single ICA run, so it cannot discriminate the arms and enters no A-vs-B claim.
long$intra_note <- "non-comparative (B = 1.000 by construction)"
write.csv(long, file.path(OUT_DIR, "AB_metrics_long.csv"), row.names = FALSE)

# ---- head-to-head on the primary metric ----
wide <- function(metric, sd_col) {
  a <- long[long$Method == "consICA", c("slide", "Mode", metric, sd_col)]
  b <- long[long$Method == "consICA_ntry1", c("slide", "Mode", metric, sd_col)]
  names(a)[3:4] <- c("A", "A_sd"); names(b)[3:4] <- c("B", "B_sd")
  w <- merge(a, b, by = c("slide", "Mode"))
  w$diff_A_minus_B <- w$A - w$B
  w$favours <- ifelse(abs(w$diff_A_minus_B) < 0.005, "~equal",
                      ifelse(w$diff_A_minus_B > 0, "A (consensus)", "B (single run)"))
  w[order(w$Mode, w$slide), ]
}

ms <- wide("Matched_Stability", "SD_Matched")
write.csv(ms, file.path(OUT_DIR, "AB_matched_stability.csv"), row.names = FALSE)

cat("=== PRIMARY: Matched stability, A (consensus) vs B (single run) ===\n")
cat("    mean over the 45 within-arm repeat pairs, +/- SD across those pairs\n\n")
for (md in MODES) {
  d <- ms[ms$Mode == md, ]
  cat("--", md, "--\n")
  print(data.frame(slide = d$slide,
                   A = sprintf("%.3f +/- %.3f", d$A, d$A_sd),
                   B = sprintf("%.3f +/- %.3f", d$B, d$B_sd),
                   `A-B` = sprintf("%+.3f", d$diff_A_minus_B),
                   favours = d$favours, check.names = FALSE), row.names = FALSE)
  cat("\n")
}

cat("=== direction summary (|difference| < 0.005 counted as equal) ===\n")
print(table(Mode = ms$Mode, favours = ms$favours))

# ---- secondary metrics, same layout ----
for (mm in list(c("Mean_Stability", "SD_Stability"), c("TopGene_Jaccard", "SD_TopGene"))) {
  w <- wide(mm[1], mm[2])
  write.csv(w, file.path(OUT_DIR, paste0("AB_", tolower(mm[1]), ".csv")), row.names = FALSE)
  cat(sprintf("\n=== SECONDARY: %s ===\n", mm[1]))
  print(data.frame(slide = w$slide, mode = as.character(w$Mode),
                   A = round(w$A, 3), B = round(w$B, 3),
                   `A-B` = sprintf("%+.3f", w$diff_A_minus_B),
                   favours = w$favours, check.names = FALSE), row.names = FALSE)
}

repro <- long[, c("slide", "Mode", "arm", "N_Repro_50", "N_Repro_70")]
write.csv(repro, file.path(OUT_DIR, "AB_n_reproducible.csv"), row.names = FALSE)

cat("\nWrote into", OUT_DIR, ":\n",
    " AB_metrics_long.csv, AB_matched_stability.csv,\n",
    " AB_mean_stability.csv, AB_topgene_jaccard.csv, AB_n_reproducible.csv\n")
cat("metrics_FINAL__*.csv were read only, not modified.\n")
