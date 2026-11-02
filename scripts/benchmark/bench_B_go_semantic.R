# Semantic uniqueness and BH interpretability for control B (consICA, ntry = 1), for the
# A against B comparison rather than against the spot-axis configurations.
# Scope: one representative per slide, arm B, pure_stability, iteration 1. Five runs.
#
# top_n follows GO_pipeline.R's table: 40 on the Xenium panels, 100 on Visium. Both
# quantities come from the uniform topGO path over S_matrix, so getGO is not called.
#
# The functions are GO_functions.R's own, so orientation, selection, topGO settings and
# the semantic measure match A: positive tail by skewness, weight01 + Fisher, nodeSize = 10,
# BP, background = the slide's own gene set, org.Hs.eg.db, GOSemSim Wang/BMA over the 10
# strongest terms per component.
#
# A's files are read only. New output goes to output/benchmark_update_2026-09-23/.

suppressPackageStartupMessages({
  library(topGO)
  library(org.Hs.eg.db)
  library(GOSemSim)
})

GO_SRC <- "scripts/benchmark/GO_functions.R"
RUN_DIR <- "output/ica_perturbations"
B_DIR <- "output/benchmark_update_2026-09-23"
K <- 20L
N_ITER <- 10L
ITER <- 1L
P_THR <- 0.05
TOP_TERMS <- 10L # GOSemSim term set size, the default A used

# A never ran GO on NCBI600 (no slide_id line in GO_pipeline.R's driver), so its semantic
# value is missing. Computing it for A as well keeps the paired comparison on all five
# slides and costs one extra slide. Set FALSE to restrict B to A's existing four.
FILL_NCBI600_FOR_A <- TRUE

# per-slide top_n exactly as GO_pipeline.R
SLIDES <- list(TENX115 = 40L, TENX117 = 40L, MEND40 = 100L, NCBI776 = 100L, NCBI600 = 100L)

# ---- the GO protocol, definitions only ----
source(GO_SRC)
stopifnot(is.function(run_GO_benchmark), is.function(add_semantic_uniqueness))
cat(sprintf("Loaded GO protocol from %s (definitions only)\n", GO_SRC))

# ---- full per-term p vectors for BH (GenTable truncates to 20 rows, so a separate call) ----
topGO_full <- function(top_genes, all_genes) {
  if (!length(top_genes)) return(NULL)
  gl <- factor(as.integer(all_genes %in% top_genes)); names(gl) <- all_genes
  gd <- suppressMessages(new("topGOdata", ontology = "BP", allGenes = gl,
                             geneSel = function(x) x == 1, annot = annFUN.org,
                             mapping = "org.Hs.eg.db", ID = "symbol", nodeSize = 10L))
  cls <- suppressMessages(runTest(gd, algorithm = "classic",  statistic = "fisher"))
  w01 <- suppressMessages(runTest(gd, algorithm = "weight01", statistic = "fisher"))
  list(m = length(usedGO(gd)), p_classic = score(cls), p_weight = score(w01))
}

S_of <- function(dir, slide, method) {
  f <- file.path(dir, sprintf("runs__%s__%s__pure_stability__k%d__n%d.rds",
                              slide, method, K, N_ITER))
  if (!file.exists(f)) { message("  missing: ", basename(f)); return(NULL) }
  S <- readRDS(f)$runs[[ITER]]$S_matrix; gc(); S
}

cat("Building GO semantic data (slow, once)...\n")
semData <- godata("org.Hs.eg.db", ont = "BP")

jobs <- c(lapply(names(SLIDES), function(s) list(slide = s, method = "consICA_ntry1",
                                                 dir = B_DIR,  arm = "B")),
          if (FILL_NCBI600_FOR_A) list(list(slide = "NCBI600", method = "consICA",
                                            dir = RUN_DIR, arm = "A")))

sem_rows <- list()
comp_rows <- list()
for (j in jobs) {
  S <- S_of(j$dir, j$slide, j$method)
  if (is.null(S)) next
  tn <- SLIDES[[j$slide]]
  cat(sprintf("\n=== %s | arm %s | %d genes x %d comps | top_n = %d ===\n",
              j$slide, j$arm, nrow(S), ncol(S), tn))

  # --- semantic uniqueness, identical path to A ----
  runs_list <- setNames(list(list(S_matrix = S, consica_object = NULL)), j$method)
  g_top <- run_GO_benchmark(runs_list, selection = "topN",   top_n_genes = tn, p_threshold = P_THR)
  g_zsc <- run_GO_benchmark(runs_list, selection = "zscore", z_thresh = 3,     p_threshold = P_THR)
  sem_rows[[length(sem_rows) + 1L]] <- data.frame(
    Method       = j$method,
    SemUniq_topN = unname(add_semantic_uniqueness(g_top, semData, TOP_TERMS)),
    SemUniq_zscr = unname(add_semantic_uniqueness(g_zsc, semData, TOP_TERMS)),
    Slide        = j$slide, arm = j$arm, top_n = tn, stringsAsFactors = FALSE)

  # --- BH interpretability, primary selection only ----
  genes_all <- rownames(S)
  sets <- select_genes_per_component(S, selection = "topN", top_n = tn)
  for (cmp in names(sets)) {
    go <- tryCatch(topGO_full(sets[[cmp]], genes_all),
                   error = function(e) { message("  ", cmp, ": ", e$message); NULL })
    comp_rows[[length(comp_rows) + 1L]] <- if (is.null(go)) data.frame(
        Slide = j$slide, arm = j$arm, Method = j$method, Component = cmp,
        m_tested = NA_integer_, minBH_classic = NA_real_, minBH_weight01 = NA_real_,
        stringsAsFactors = FALSE) else data.frame(
        Slide = j$slide, arm = j$arm, Method = j$method, Component = cmp,
        m_tested = go$m,
        minBH_classic  = min(p.adjust(go$p_classic, "BH"), na.rm = TRUE),
        minBH_weight01 = min(p.adjust(go$p_weight,  "BH"), na.rm = TRUE),
        stringsAsFactors = FALSE)
  }
  rm(S, g_top, g_zsc); gc()
}

# ---- outputs, A's files untouched ----
sem <- do.call(rbind, sem_rows)
write.csv(sem, file.path(B_DIR, "go_semantic_uniqueness_B__k20.csv"), row.names = FALSE)

comp <- do.call(rbind, comp_rows)
comp$interp_BH_classic <- comp$minBH_classic < P_THR
write.csv(comp, file.path(B_DIR, "go_interpretability_bh_B__per_component__k20.csv"), row.names = FALSE)

summ <- do.call(rbind, lapply(split(comp, list(comp$Slide, comp$arm), drop = TRUE), function(d)
  data.frame(Slide = d$Slide[1], arm = d$arm[1], Method = d$Method[1], n_comp = nrow(d),
             Interp_Frac_BH_classic = mean(d$minBH_classic < P_THR, na.rm = TRUE),
             Interp_Frac_BH_weight01 = mean(d$minBH_weight01 < P_THR, na.rm = TRUE),
             worst_minBH_classic = max(d$minBH_classic, na.rm = TRUE),
             m_tested = d$m_tested[1], stringsAsFactors = FALSE)))
write.csv(summ, file.path(B_DIR, "go_interpretability_bh_B__summary__k20.csv"), row.names = FALSE)

cat("\n=== B: semantic uniqueness ===\n");        print(sem,  row.names = FALSE, digits = 4)
cat("\n=== B: BH interpretability + m ===\n");    print(summ, row.names = FALSE, digits = 4)

# ---- head-to-head against A where A exists ----
af <- list.files("output/go_enrichment", pattern = "^go_semantic_uniqueness__.*__k20\\.csv$",
                 full.names = TRUE)
if (length(af)) {
  A_old <- do.call(rbind, lapply(af, read.csv, stringsAsFactors = FALSE))
  A_old <- A_old[A_old$Method == "consICA", c("Slide", "SemUniq_topN", "SemUniq_zscr")]
  # arm-A slides computed here because GO_pipeline.R's driver never ran them (NCBI600).
  # Without this the comparison shows NA for them although the value exists.
  A_new <- sem[sem$arm == "A", c("Slide", "SemUniq_topN", "SemUniq_zscr")]
  A <- rbind(A_old, A_new[!A_new$Slide %in% A_old$Slide, , drop = FALSE])
  names(A)[2:3] <- c("A_topN", "A_zscr")
  Bs <- sem[sem$arm == "B", c("Slide", "SemUniq_topN", "SemUniq_zscr")]
  names(Bs)[2:3] <- c("B_topN", "B_zscr")
  cmp <- merge(A, Bs, by = "Slide", all = TRUE)
  cmp$diff_topN <- cmp$A_topN - cmp$B_topN
  cat("\n=== A vs B: semantic uniqueness (higher = more distinct biology) ===\n")
  print(cmp, row.names = FALSE, digits = 4)
  write.csv(cmp, file.path(B_DIR, "AB_semantic_uniqueness.csv"), row.names = FALSE)
}
ai <- "output/benchmark_summary/go_interpretability_bh__summary__k20.csv"
if (file.exists(ai)) {
  A2 <- read.csv(ai, stringsAsFactors = FALSE)
  A2 <- A2[A2$Method == "consICA", c("Slide", "Interp_Frac_BH_classic", "m_med")]
  names(A2)[2:3] <- c("A_Interp_BH", "A_m")
  B2 <- summ[summ$arm == "B", c("Slide", "Interp_Frac_BH_classic", "m_tested")]
  names(B2)[2:3] <- c("B_Interp_BH", "B_m")
  cmp2 <- merge(A2, B2, by = "Slide")
  cmp2$diff <- cmp2$A_Interp_BH - cmp2$B_Interp_BH
  cmp2$m_agrees <- cmp2$A_m == cmp2$B_m
  cat("\n=== A vs B: BH interpretability, and whether m matches ===\n")
  print(cmp2, row.names = FALSE, digits = 4)
  write.csv(cmp2, file.path(B_DIR, "AB_interpretability_bh.csv"), row.names = FALSE)
}
cat("\nA's go_enrichment / benchmark_summary files were read only, not modified.\n")
