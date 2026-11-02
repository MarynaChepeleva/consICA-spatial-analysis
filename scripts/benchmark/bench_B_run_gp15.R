# Control B at 15% gene dropout, where arm A already exists.
# Scope: MEND40 and NCBI600, gene_perturbation only, 10 repeats = 20 decompositions.
#
# Sources bench_B_setup.R, the same setup bench_B_run.R uses: load_master with its
# dimension gate against A, the consICA version gate, and the benchmark engine. The gp05 run
# itself is untouched; it already ran and carries the provenance of those 150 decompositions.
#
# Seeds are seq_len(10) inside execute_benchmark_mode, identical to A's gp15 run, and the
# gene subsets are drawn from the same master, so B perturbs the same genes A did.
#
# Writes only into output/benchmark_update_2026-09-23/ under gp15 names. A's
# metrics_FINAL__*__gp15pct.csv are read for the head-to-head, never written.

source("scripts/benchmark/bench_B_setup.R")   # gives load_master, cfg_all, prov, METHOD, OUT_DIR
stopifnot(is.function(load_master), is.function(execute_benchmark_mode))

SLIDES_15 <- c("MEND40", "NCBI600")
MODE <- "gene_perturbation"
N_ITER <- 10L
GP_PCT <- 15

cat(sprintf("\n=== B at 15%% gene dropout: %d slides x 1 mode x %d repeats ===\n\n",
            length(SLIDES_15), N_ITER))

manifest <- list()
for (sl in SLIDES_15) {
  cfg <- cfg_all[cfg_all$run_id == paste0(sl, "__k20__gp15"), ]
  stopifnot(nrow(cfg) == 1L)
  stopifnot(abs(cfg$gene_frac - 0.15) < 1e-9) # guard: this must be the gp15 config
  cat("[", sl, "] loading\n", sep = "")
  M <- load_master(cfg)

  cat(sprintf("  -> %s (gene_frac=%.2f, n_iter=%d, ntry=1, ncores=1)\n",
              MODE, cfg$gene_frac, N_ITER))
  warns <- character(0); t0 <- Sys.time()
  result <- withCallingHandlers(
    tryCatch(
      execute_benchmark_mode(
        srt_raw = M$srt, mode = MODE, k = 20L, n_iter = N_ITER, method = METHOD,
        gene_perturb_frac = cfg$gene_frac, spot_perturb_frac = cfg$spot_frac,
        ncores = 1L, top_n = 50L, repro_thresholds = c(0.5, 0.7),
        keep_all_consica = FALSE),
      error = function(e) { warns <<- c(warns, paste("ERROR:", conditionMessage(e))); NULL }),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (is.null(result)) {
    manifest[[length(manifest) + 1L]] <- data.frame(
      slide = sl, mode = MODE, repeat_i = NA_integer_, seed = NA_integer_, status = "FAILED",
      n_genes_eff = NA_integer_, n_spots_eff = NA_integer_, n_subset = NA_integer_,
      elapsed_sec = elapsed, ntry = 1L, ncores = 1L, gene_perturb_pct = GP_PCT,
      messages = paste(unique(warns), collapse = " | "), stringsAsFactors = FALSE)
    cat("     FAILED:", paste(unique(warns), collapse = " | "), "\n")
    rm(M)
    gc()
    next
  }

  for (i in seq_along(result$runs)) result$runs[[i]]$consica_object <- NULL   # S and M only

  meta <- make_run_metadata(
    slide_id = sl, test_file = cfg$file, method = METHOD, mode = MODE, k = 20L,
    n_iter = N_ITER, gene_perturb_frac = cfg$gene_frac, spot_perturb_frac = cfg$spot_frac,
    ncores = 1L, top_n = 50L, repro_thresholds = c(0.5, 0.7),
    min_genes = M$min_genes, min_counts = M$min_counts, srt_master = M$srt)

  tag <- make_save_tag(sl, METHOD, MODE, 20L, N_ITER, GP_PCT)
  saveRDS(list(metadata = meta, metrics = result$metrics,
               runs = result$runs, diagnostics = result$diagnostics),
          file.path(OUT_DIR, paste0("runs__", tag, ".rds")))

  for (i in seq_along(result$runs)) {
    r <- result$runs[[i]]
    manifest[[length(manifest) + 1L]] <- data.frame(
      slide = sl, mode = MODE, repeat_i = i, seed = i, status = "OK",
      n_genes_eff = r$n_genes_eff, n_spots_eff = r$n_spots_eff,
      n_subset = length(r$subset), elapsed_sec = r$elapsed_sec,
      ntry = 1L, ncores = 1L, gene_perturb_pct = GP_PCT,
      messages = paste(unique(warns), collapse = " | "), stringsAsFactors = FALSE)
  }
  cat(sprintf("     saved %s (%.0f s)\n", tag, elapsed))
  rm(M, result); gc()
}

man <- do.call(rbind, manifest)
write.csv(man, file.path(OUT_DIR, "manifest_B_gp15.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo_B_gp15.txt"))

cat("\n=== DONE ===\n")
print(man[, c("slide", "repeat_i", "status", "n_genes_eff", "n_spots_eff", "elapsed_sec")],
      row.names = FALSE)
cat("\nDecompositions:", sum(man$status == "OK"), "| failed:", sum(man$status == "FAILED"), "\n")

# ---- head-to-head with A at 15% (A's files read only) ----
bm <- do.call(rbind, lapply(SLIDES_15, function(sl) {
  f <- file.path(OUT_DIR, sprintf("runs__%s__%s__%s__k20__n%d_gp%dpct.rds",
                                  sl, METHOD, MODE, N_ITER, GP_PCT))
  if (!file.exists(f)) return(NULL)
  m <- readRDS(f)$metrics; m$slide <- sl; m
}))
am <- do.call(rbind, lapply(SLIDES_15, function(sl) {
  f <- file.path("output/ica_perturbations",
                 sprintf("metrics_FINAL__%s__k20__n10__gp15pct.csv", sl))
  if (!file.exists(f)) { message("A gp15 metrics missing: ", basename(f)); return(NULL) }
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  d <- d[d$Method == "consICA" & d$Mode == MODE, ]; d$slide <- sl; d
}))

if (!is.null(bm) && !is.null(am)) {
  cmp <- merge(am[, c("slide", "Matched_Stability", "SD_Matched", "Mean_Stability", "TopGene_Jaccard")],
               bm[, c("slide", "Matched_Stability", "SD_Matched", "Mean_Stability", "TopGene_Jaccard")],
               by = "slide", suffixes = c("_A", "_B"))
  cmp$diff_matched <- cmp$Matched_Stability_A - cmp$Matched_Stability_B
  write.csv(cmp, file.path(OUT_DIR, "AB_matched_stability_gp15.csv"), row.names = FALSE)
  cat("\n=== A vs B at 15% gene dropout ===\n")
  print(cmp, row.names = FALSE, digits = 4)
  cat("\nFor reference, the same comparison at 5% gave A - B = +0.069 (MEND40),",
      "+0.026 (NCBI600), from AB_matched_stability.csv.\n")
  cat("CAVEAT: A's gp15 cells on both slides ran with ncores = 4; B is single-core throughout.\n")
}
cat("\nWrote manifest_B_gp15.csv, sessionInfo_B_gp15.txt, AB_matched_stability_gp15.csv in",
    OUT_DIR, "\n")
