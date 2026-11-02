# Computes control B: consICA with ntry = 1, a single-run ICA on the consICA core.
# Same engine, same orientation, same inputs, same seeds as A; the ONLY difference is one
# ICA run instead of 20 averaged, so A - B isolates the consensus procedure.
#
# Writes ONLY into output/benchmark_update_2026-09-23/. Nothing existing is overwritten:
# the driver's `methods` vector and output/ica_perturbations/ are untouched.
#
# Preflight (bench_B_preflight.R) established that B reproduces A's subsets exactly
# and that spot-mode SCTransform is bitwise deterministic. Run that first.
#
# The engine, arm A's provenance, the consICA version gate and load_master() are in
# bench_B_setup.R, sourced below.
#
# execute_benchmark_mode always computes the stability metrics as a by-product. They are
# saved but not analysed here.

source("scripts/benchmark/bench_B_setup.R")

# ---- what to run ----
TEST_RUN <- F    # TRUE = TENX117 / pure_stability / 2 repeats only. Flip to FALSE for all 150.
MODES <- c("pure_stability", "spot_perturbation", "gene_perturbation")
SLIDES <- c("TENX115", "TENX117", "MEND40", "NCBI776", "NCBI600")

# ---- run -----
slides_todo <- if (TEST_RUN) "TENX117" else SLIDES
modes_todo <- if (TEST_RUN) "pure_stability" else MODES
cat(sprintf("\n=== %s: %d slide(s) x %d mode(s) ===\n\n",
            if (TEST_RUN) "TEST RUN" else "FULL RUN B", length(slides_todo), length(modes_todo)))

manifest <- list()
for (sl in slides_todo) {
  cfg <- cfg_all[cfg_all$run_id == paste0(sl, "__k20__gp05"), ]
  stopifnot(nrow(cfg) == 1L)
  n_iter <- if (TEST_RUN) 2L else as.integer(cfg$n_iter)
  cat("[", sl, "] loading\n", sep = "")
  M <- load_master(cfg)

  for (md in modes_todo) {
    cat(sprintf("  -> %s (n_iter=%d, ntry=1, ncores=1)\n", md, n_iter))
    warns <- character(0)
    t0 <- Sys.time()
    result <- withCallingHandlers(
      tryCatch(
        execute_benchmark_mode(
          srt_raw = M$srt, mode = md, k = 20L, n_iter = n_iter, method = METHOD,
          gene_perturb_frac = cfg$gene_frac, spot_perturb_frac = cfg$spot_frac,
          ncores = 1L,                       # explicit: nothing to parallelise at ntry = 1
          top_n = 50L, repro_thresholds = c(0.5, 0.7),
          keep_all_consica = FALSE),
        error = function(e) { warns <<- c(warns, paste("ERROR:", conditionMessage(e))); NULL }),
      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (is.null(result)) {                    # record the failure, do not drop it
      manifest[[length(manifest) + 1L]] <- data.frame(
        slide = sl, mode = md, repeat_i = NA_integer_, seed = NA_integer_,
        status = "FAILED", n_genes_eff = NA_integer_, n_spots_eff = NA_integer_,
        n_subset = NA_integer_, elapsed_sec = elapsed, ntry = 1L, ncores = 1L,
        consICA = now_ver, messages = paste(unique(warns), collapse = " | "),
        stringsAsFactors = FALSE)
      cat("     FAILED:", paste(unique(warns), collapse = " | "), "\n")
      next
    }

    # B needs S and M only; the consICA object holds a dense copy of X (~1.4 GB on Visium)
    for (i in seq_along(result$runs)) result$runs[[i]]$consica_object <- NULL

    meta <- make_run_metadata(
      slide_id = sl, test_file = cfg$file, method = METHOD, mode = md, k = 20L,
      n_iter = n_iter, gene_perturb_frac = cfg$gene_frac, spot_perturb_frac = cfg$spot_frac,
      ncores = 1L, top_n = 50L, repro_thresholds = c(0.5, 0.7),
      min_genes = M$min_genes, min_counts = M$min_counts, srt_master = M$srt)

    gp <- if (md == "gene_perturbation") cfg$gene_frac * 100 else NA
    tag <- make_save_tag(sl, METHOD, md, 20L, n_iter, gp)
    saveRDS(list(metadata = meta, metrics = result$metrics,
                 runs = result$runs, diagnostics = result$diagnostics),
            file.path(OUT_DIR, paste0("runs__", tag, ".rds")))

    for (i in seq_along(result$runs)) {
      r <- result$runs[[i]]
      manifest[[length(manifest) + 1L]] <- data.frame(
        slide = sl, mode = md, repeat_i = i, seed = i, status = "OK",
        n_genes_eff = r$n_genes_eff, n_spots_eff = r$n_spots_eff,
        n_subset = if (length(r$subset) == 1L && is.na(r$subset)) NA_integer_ else length(r$subset),
        elapsed_sec = r$elapsed_sec, ntry = 1L, ncores = 1L, consICA = now_ver,
        messages = paste(unique(warns), collapse = " | "), stringsAsFactors = FALSE)
    }
    cat(sprintf("     saved %s (%.0f s)\n", tag, elapsed))
  }
  rm(M); gc()
}

man <- do.call(rbind, manifest)
man_file <- file.path(OUT_DIR, if (TEST_RUN) "manifest_B_test.csv" else "manifest_B.csv")
write.csv(man, man_file, row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, if (TEST_RUN) "sessionInfo_B_test.txt" else "sessionInfo_B.txt"))

cat("\n=== DONE ===\n")
print(man[, c("slide", "mode", "repeat_i", "status", "n_genes_eff", "n_spots_eff", "elapsed_sec")],
      row.names = FALSE)
cat("\nDecompositions:", sum(man$status == "OK"), "| failed:", sum(man$status == "FAILED"), "\n")
cat("Manifest:", normalizePath(man_file, mustWork = FALSE), "\n")
if (TEST_RUN) cat("\nTEST_RUN passed -> set TEST_RUN <- FALSE and re-run for the full 150.\n")
