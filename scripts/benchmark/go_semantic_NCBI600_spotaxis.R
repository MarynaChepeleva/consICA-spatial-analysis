# Semantic uniqueness for the three spot-axis configurations on NCBI600, the slide
# GO_pipeline.R's table does not cover. It runs the same functions, from GO_functions.R,
# with the parameters that table gives the other Visium slides. Adding the slide there
# would replace this step.
#
# Writes:
#   output/go_enrichment/go_semantic_uniqueness__NCBI600__k20.csv   (the three icafast rows)
#   output/go_enrichment/go_semantic_spotaxis__NCBI600__k20.rds     (the term lists, for reuse)

suppressPackageStartupMessages({
  library(topGO)
  library(org.Hs.eg.db)
  library(GOSemSim)
})

GO_SRC <- "scripts/benchmark/GO_functions.R"
RUN_DIR <- "output/ica_perturbations"
GO_DIR <- "output/go_enrichment"
OUT_TXT <- file.path(GO_DIR, "go_semantic_NCBI600_spotaxis.log.txt")

SLIDE <- "NCBI600"
K <- 20L
N_ITER <- 10L
ITER <- 1L # the representative pure_stability iteration, as in GO_pipeline.R
TOP_N <- 100L # transcriptome-wide Visium, as the SLIDES table of GO_pipeline.R
Z_THRESH <- 3
P_THR <- 0.05
TOP_TERMS <- 10L # GOSemSim term set size, the default A used

# written out
METHODS <- c("CFS_As_Is", "icafast_no_Rmat", "icafast_with_Rmat")

# computed as a control and NOT written: bench_B_go_semantic.R already produced consICA on
# this slide, and tableS_comparison_full.R lets a file in GO_DIR win over that one. Recomputing it
# here only to confirm the path reproduces the published value keeps the table's existing
# number untouched whatever happens.
CONTROL <- "consICA"
CONTROL_EXPECT <- c(topN = 0.772978947368421, zscr = 0.777942105263158)
TOL <- 1e-9

# no on.exit here: it is a function-level construct, and at the top level of a sourced script
# it fires as soon as the expression that registered it finishes, which closes the connection
# immediately and leaves an empty log. Teardown is explicit at the foot of the file instead.
con <- file(OUT_TXT, open = "wt"); sink(con, split = TRUE)

cat(sprintf("semantic uniqueness, %s, k = %d, started %s\n\n", SLIDE, K, format(Sys.time())))

# ---- the GO protocol, definitions only ----
source(GO_SRC)
stopifnot(is.function(run_GO_benchmark), is.function(add_semantic_uniqueness),
          is.function(load_pure_stability_runs))
cat(sprintf("GO protocol loaded from %s (definitions only)\n", GO_SRC))

S_of <- function(method) {
  f <- file.path(RUN_DIR, sprintf("runs__%s__%s__pure_stability__k%d__n%d.rds",
                                  SLIDE, method, K, N_ITER))
  if (!file.exists(f)) stop("run file not found: ", f)
  S <- readRDS(f)$runs[[ITER]]$S_matrix
  gc(verbose = FALSE)
  S
}

cat("Building GO semantic data (slow, once)...\n")
t0 <- Sys.time()
semData <- godata("org.Hs.eg.db", ont = "BP")
cat(sprintf("  godata in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

rows <- list(); keep <- list()
for (m in c(METHODS, CONTROL)) {
  S <- S_of(m)
  t1 <- Sys.time()
  cat(sprintf("\n=== %s | %s | %d genes x %d components | top_n = %d ===\n",
              SLIDE, m, nrow(S), ncol(S), TOP_N))

  runs_list <- setNames(list(list(S_matrix = S, consica_object = NULL)), m)
  g_top <- run_GO_benchmark(runs_list, selection = "topN",   top_n_genes = TOP_N,
                            ontology = "BP", orgdb = "org.Hs.eg.db", p_threshold = P_THR)
  g_zsc <- run_GO_benchmark(runs_list, selection = "zscore", z_thresh = Z_THRESH,
                            ontology = "BP", orgdb = "org.Hs.eg.db", p_threshold = P_THR)

  rows[[m]] <- data.frame(
    Method       = m,
    SemUniq_topN = unname(add_semantic_uniqueness(g_top, semData, TOP_TERMS)),
    SemUniq_zscr = unname(add_semantic_uniqueness(g_zsc, semData, TOP_TERMS)),
    Slide        = SLIDE, stringsAsFactors = FALSE)
  keep[[m]] <- list(topN = g_top, zscore = g_zsc)

  cat(sprintf("  SemUniq topN %.6f | zscore %.6f | %.1f min\n",
              rows[[m]]$SemUniq_topN, rows[[m]]$SemUniq_zscr,
              as.numeric(difftime(Sys.time(), t1, units = "mins"))))
  rm(S, runs_list, g_top, g_zsc); gc(verbose = FALSE)
}

# ---- control: does this path reproduce the published consICA value on this slide? ---------
ctl <- rows[[CONTROL]]
dtop <- abs(ctl$SemUniq_topN - CONTROL_EXPECT[["topN"]])
dzsc <- abs(ctl$SemUniq_zscr - CONTROL_EXPECT[["zscr"]])
cat(sprintf("\ncontrol, %s on %s against bench_B_go_semantic.R:\n", CONTROL, SLIDE))
cat(sprintf("  topN   recomputed %.15f  published %.15f  difference %.2e\n",
            ctl$SemUniq_topN, CONTROL_EXPECT[["topN"]], dtop))
cat(sprintf("  zscore recomputed %.15f  published %.15f  difference %.2e\n",
            ctl$SemUniq_zscr, CONTROL_EXPECT[["zscr"]], dzsc))
ok <- dtop < TOL && dzsc < TOL
cat(sprintf("  %s\n", if (ok) "PASS: same path, same numbers; the three new values are comparable"
                      else "FAIL: the path differs from the one that produced the other slides"))

# ---- write -----
sem <- do.call(rbind, rows[METHODS])# the control is deliberately not written
rownames(sem) <- NULL
cat("\nwritten rows:\n"); print(sem, row.names = FALSE)

if (!ok) {
  cat("\nNOT WRITING: the control did not reproduce. Nothing on disk has changed.\n")
} else {
  out_csv <- file.path(GO_DIR, sprintf("go_semantic_uniqueness__%s__k%d.csv", SLIDE, K))
  if (file.exists(out_csv)) stop("refusing to overwrite an existing file: ", out_csv)
  write.csv(sem, out_csv, row.names = FALSE)
  cat(sprintf("\nwrote %s\n", out_csv))

  # the term lists, so this never has to run again. Deliberately NOT named go_benchmark__*,
  # which component_gene_overlap.R and go_interpretability_fdr_check.R pick up by glob and
  # would then meet an object without the consICA-native slots those two expect.
  out_rds <- file.path(GO_DIR, sprintf("go_semantic_spotaxis__%s__k%d.rds", SLIDE, K))
  saveRDS(list(per_method = keep, semantic_uniqueness = sem,
               top_n = TOP_N, z_thresh = Z_THRESH, top_terms = TOP_TERMS,
               iter = ITER, built = Sys.time()), out_rds)
  cat(sprintf("wrote %s\n", out_rds))
  cat("\nnext: rerun scripts/tableS_comparison_full.R, then scripts/assemble_supplement.R\n")
}

cat(sprintf("\nfinished %s\n", format(Sys.time())))
sink()
close(con)
message("log written: ", OUT_TXT)
