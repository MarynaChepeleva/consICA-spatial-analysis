# GO enrichment per component, the interpretability side of the ICA benchmark. Runs over the
# pure_stability decompositions saved by bench_ica_perturbations.R, one representative
# iteration per method and slide, and writes go_benchmark__*.rds and three flat CSVs per
# slide into output/go_enrichment/.
#
# NCBI600 is not in the slide table below. It entered the benchmark after this driver was
# written, and its GO was filled in by bench_B_go_semantic.R and go_semantic_NCBI600_spotaxis.R;
# adding it here would not reproduce those files.
#
# Gene selection, topGO and the metrics are in GO_functions.R.

source("scripts/benchmark/GO_functions.R")

# ---- configuration ----

output_dir <- "output/ica_perturbations"
go_dir <- "output/go_enrichment"
dir.create(go_dir, showWarnings = FALSE, recursive = TRUE)

methods <- c("consICA", "CFS_As_Is",
             "icafast_no_Rmat", "icafast_with_Rmat", "fastICA_Package")

# top_n_genes follows the platform: a Xenium targeted panel carries ~300 genes, so 40 is the
# usable top, while transcriptome-wide Visium uses 100.
SLIDES <- data.frame(
  slide_id = c("TENX115", "TENX117", "MEND40", "NCBI776"),
  k_components = c(20L, 20L, 20L, 20L),
  n_iterations = c(10L, 10L, 10L, 10L),
  top_n_genes = c(40L, 40L, 100L, 100L),
  stringsAsFactors = FALSE)

# The organism is chosen here and nowhere else: topGO and GOSemSim each load the package
# from this string, so a mouse xenograft only sets it to "org.Mm.eg.db" and the package has
# to be installed, not attached. consICA::getGO calls requireNamespace("org.Hs.eg.db")
# whatever the genome, so the human package stays a dependency either way.
orgdb <- "org.Hs.eg.db"
ontology <- "BP"
z_thresh <- 3
p_thr <- 0.05

# ---- one slide per row ----

for (s in seq_len(nrow(SLIDES))) {

  slide_id <- SLIDES$slide_id[s]
  k_components <- SLIDES$k_components[s]
  n_iterations <- SLIDES$n_iterations[s]
  top_n_genes <- SLIDES$top_n_genes[s]

  cat(sprintf("\n======== %s (%d/%d) ========\n", slide_id, s, nrow(SLIDES)))

  cat("Loading pure_stability decompositions...\n")
  method_runs <- load_pure_stability_runs(
    output_dir, slide_id, methods, k_components, n_iterations, rep_iter = 1L)

  # ---- (1) uniform topGO: primary (topN) + sensitivity (zscore) ----
  go_topN <- run_GO_benchmark(method_runs, selection = "topN",
                              top_n_genes = top_n_genes, ontology = ontology,
                              orgdb = orgdb, p_threshold = p_thr)

  go_zscore <- run_GO_benchmark(method_runs, selection = "zscore",
                                z_thresh = z_thresh, ontology = ontology,
                                orgdb = orgdb, p_threshold = p_thr)

  summary_all <- rbind(go_topN$summary, go_zscore$summary)
  cat("\n=== INTERPRETABILITY SUMMARY (provisional metrics) ===\n")
  print(summary_all)

  # ---- (2) native consICA getGO ----
  cat("\nRunning native consICA::getGO...\n")
  consica_go <- run_consica_native_GO(
    method_runs$consICA$consica_object, genome = orgdb, db = ontology, ncores = 4L)

  # ---- (3) semantic-similarity metrics (GOSemSim) ----
  cat("\nBuilding GO semantic data (godata)...\n")
  semData <- godata(orgdb, ont = ontology)          # reused across methods

  sem_uniq <- data.frame(
    Method       = names(go_topN$per_method),
    SemUniq_topN = add_semantic_uniqueness(go_topN,   semData),
    SemUniq_zscr = add_semantic_uniqueness(go_zscore, semData),
    row.names = NULL)
  cat("\n=== SEMANTIC UNIQUENESS (biology-level diversity) ===\n")
  print(sem_uniq)

  # semantic concordance (consICA): uniform vs native getGO$pos, per component
  concordance <- compute_go_concordance_semantic(
    go_topN$per_method$consICA$go_results, consica_go$GO$GOBP, semData)
  cat(sprintf("\nconsICA semantic concordance (uniform vs native getGO): mean = %.3f\n",
              mean(concordance$GO_Concordance_semantic, na.rm = TRUE)))

  # ---- save ----
  saveRDS(list(topN = go_topN,
               zscore = go_zscore,
               summary = summary_all,
               semantic_uniqueness = sem_uniq,
               consica_native = consica_go,
               concordance = concordance),
          file.path(go_dir, sprintf("go_benchmark__%s__k%d.rds",
                                    slide_id, k_components)))
  write.csv(summary_all,
            file.path(go_dir, sprintf("go_summary__%s__k%d.csv",
                                      slide_id, k_components)),
            row.names = FALSE)

  # flat CSVs for cross-slide aggregation
  sem_uniq$Slide <- slide_id
  write.csv(sem_uniq,
            file.path(go_dir, sprintf("go_semantic_uniqueness__%s__k%d.csv",
                                      slide_id, k_components)),
            row.names = FALSE)
  concordance$Slide <- slide_id
  write.csv(concordance,
            file.path(go_dir, sprintf("go_concordance__%s__k%d.csv",
                                      slide_id, k_components)),
            row.names = FALSE)

  cat("\nSaved GO results to", go_dir, "\n")
}

# ---- what a saved run holds ----
#   res <- readRDS("output/go_enrichment/go_benchmark__TENX115__k20.rds")
#   res$summary                              # cross-method interpretability (exact-ID)
#   res$semantic_uniqueness                  # biology-level diversity (GOSemSim)
#   res$topN$per_method$consICA$go_results   # per-component GO tables (primary)
#   res$consica_native                       # consICA cross-validated GO
#   res$concordance                          # semantic concordance uniform vs getGO
#   plots <- plot_GO_benchmark(res$topN); print(plots$p2)

# consica specific

# saveReport(consica_go, file = file.path(go_dir, sprintf("Report_cica_k%d_%s.pdf",
#                                                         k_components, slide_id)))
