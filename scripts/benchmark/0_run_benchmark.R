# Runs the benchmark end to end, one script per R session. Several of these hold a whole
# slide plus its decompositions in memory, so a fresh session per step keeps the peak
# bounded and makes a failure cost one step rather than the run.
#
# The order is not a preference. Each step reads a file an earlier one writes:
#   2 writes the decompositions everything below reads;
#   6 is gated against the CSV that 5 writes;
#   11, 12 read the Moran, semantic and BH summaries that 3, 4 and 6 write;
#   15 reads the control arm's compact runs, so it cannot precede 14.
#
# Run from the project root. Set ONLY to re-run part of it, for example
# ONLY <- STEPS[8:14] for the control arm alone.

STEPS <- c(
  # the benchmark
  "scripts/benchmark/config.R",                     # the run grid -> config/runs.csv
  "scripts/benchmark/bench_ica_perturbations.R",    # every decomposition, all configured runs
  "scripts/benchmark/bench_moran_by_method.R",      # Moran's I per method
  "scripts/benchmark/GO_pipeline.R",                # per-component GO and semantic uniqueness
  "scripts/benchmark/go_interpretability_fdr_check.R",  # Bonferroni check; its CSV gates the next step
  "scripts/benchmark/go_interpretability_bh.R",     # the same counts under BH
  "scripts/benchmark/go_semantic_NCBI600_spotaxis.R",   # the three NCBI600 cells GO_pipeline omits

  # the control arm: consICA at ntry = 1, to isolate the consensus step
  "scripts/benchmark/bench_B_preflight.R",          # provenance, subset and determinism gates
  "scripts/benchmark/bench_B_run.R",                # the 150 decompositions at ntry = 1
  "scripts/benchmark/bench_B_run_gp15.R",           # the same at 15% gene dropout
  "scripts/benchmark/bench_B_moran.R",              # Moran for the control arm
  "scripts/benchmark/bench_B_go_semantic.R",        # GO and semantic uniqueness for the control arm
  "scripts/benchmark/bench_B_intra_and_cost.R",     # intra-uniqueness and the cost ratio
  "scripts/benchmark/bench_AB_matched_stability.R", # the head-to-head table

  # analyses that need both
  "scripts/benchmark/component_gene_overlap.R",     # gene overlap between components, by configuration
  "scripts/benchmark/aggregate_results.R")          # ALL_*.csv over everything above

ONLY <- NULL     # character vector of paths to restrict to, or NULL for all
DRY_RUN <- TRUE  # TRUE prints the calls and runs nothing

todo <- if (is.null(ONLY)) STEPS else STEPS[STEPS %in% ONLY]
stopifnot(length(todo) > 0L, all(file.exists(todo)))
rscript <- file.path(R.home("bin"), "Rscript")

for (i in seq_along(todo)) {
  cat(sprintf("\n==== [%d/%d] %s\n", i, length(todo), todo[i]))
  if (DRY_RUN) next
  t0 <- Sys.time()
  if (system2(rscript, todo[i]) != 0) stop("failed at step ", i, ": ", todo[i])
  cat(sprintf("     done in %.0f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

if (DRY_RUN) cat("\nDRY_RUN is TRUE, nothing was run. Set it to FALSE to execute.\n")
