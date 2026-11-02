# Writes the three cfs_mechanism__*.csv diagnostics from the saved runs.
#   1. overlap of the top-gene lists between components of one decomposition
#   2. how strongly a component's loadings track the gene's overall profile
#   3. how many distinct GO terms the annotated components of a configuration share
# Gene lists match the enrichment selection: 40 on Xenium, 100 on Visium.
# The gene axis is reported for both arms; arm A files reach 13.2 GB, and one that cannot
# be read is left out, leaving the arm B row alone.

suppressPackageStartupMessages({ library(matrixStats) })

DO_POINT2 <- FALSE      # FALSE = points 1 and 3 only (seconds); TRUE adds ~40 min of SCTransform
A_SIZE_CAP <- Inf

RUN_DIR <- "output/ica_perturbations"
B_DIR <- "output/benchmark_update_2026-09-23"
GO_DIR <- "output/go_enrichment"
OUT_DIR <- B_DIR
K <- 20L
ITER <- 1L

# ordered by the size of arm A's file, so the heaviest read comes last and a failure there
# costs only that slide: TENX115 70 MB, TENX117 61 MB, NCBI600 1.6 GB, MEND40 6.3 GB, NCBI776 13.2 GB
SLIDES <- list(TENX115 = 40L, TENX117 = 40L, NCBI600 = 100L, MEND40 = 100L, NCBI776 = 100L)
PROJ <- c("CFS_As_Is", "icafast_no_Rmat", "icafast_with_Rmat")   # projected loadings

# gene selection identical to the GO pipeline (orientation by skewness, positive tail)
source("scripts/benchmark/GO_functions.R")
stopifnot(is.function(select_genes_per_component))

S_from <- function(path) { if (!file.exists(path)) return(NULL)
  S <- readRDS(path)$runs[[ITER]]$S_matrix; gc(); S }

# every configuration of one slide, as a named list of S matrices
configs_of <- function(slide) {
  out <- lapply(PROJ, function(m)
    S_from(file.path(RUN_DIR, sprintf("runs__%s__%s__pure_stability__k20__n10.rds", slide, m))))
  names(out) <- PROJ
  out[["consICA_ntry1 (gene axis)"]] <-
    S_from(file.path(B_DIR, sprintf("runs__%s__consICA_ntry1__pure_stability__k20__n10.rds", slide)))
  # Arm A's pure cells carry the consensus object, which holds a dense copy of X: the file
  # size is linear in X (70 MB on TENX115 -> 13.2 GB on NCBI776, a 204-fold matrix). Reading
  # is attempted anyway, but a failure degrades to arm B for that slide instead of aborting.
  fa <- file.path(RUN_DIR, sprintf("runs__%s__consICA__pure_stability__k20__n10.rds", slide))
  if (file.exists(fa) && file.size(fa) < A_SIZE_CAP) {
    cat(sprintf("    reading arm A (%.1f GB on disk)...\n", file.size(fa) / 1e9))
    Sa <- tryCatch(S_from(fa), error = function(e) {
      message("    arm A not read for ", slide, ": ", conditionMessage(e),
              " -- arm B stands in for this slide"); NULL })
    if (!is.null(Sa)) out[["consICA (gene axis, arm A)"]] <- Sa
    gc()
  }
  out[!vapply(out, is.null, logical(1))]
}

# 1. overlap of the top-gene lists between components ----
p1 <- list()
for (sl in names(SLIDES)) {
  tn <- SLIDES[[sl]]
  for (nm in names(cfg <- configs_of(sl))) {
    S <- cfg[[nm]]
    sets <- select_genes_per_component(S, selection = "topN", top_n = tn)
    pr <- utils::combn(length(sets), 2)
    jac <- apply(pr, 2, function(p) {
      a <- sets[[p[1]]]; b <- sets[[p[2]]]
      length(intersect(a, b)) / length(union(a, b))
    })
    uni <- length(unique(unlist(sets)))
    p1[[length(p1) + 1L]] <- data.frame(
      slide = sl, config = nm, loadings = if (nm %in% PROJ) "projected" else "gene-axis source",
      k = ncol(S), top_n = tn, n_pairs = ncol(pr),
      mean_pairwise_jaccard = mean(jac), max_pairwise_jaccard = max(jac),
      union_size = uni, k_times_N = ncol(S) * tn,
      union_ratio = uni / (ncol(S) * tn), stringsAsFactors = FALSE)
  }
  rm(cfg); gc()
}
P1 <- do.call(rbind, p1)
write.csv(P1, file.path(OUT_DIR, "cfs_mechanism__1_topgene_overlap.csv"), row.names = FALSE)
cat("=== 1. top-gene list overlap between components (lower union_ratio = more shared) ===\n")
print(P1[, c("slide","config","loadings","mean_pairwise_jaccard","union_ratio")],
      row.names = FALSE, digits = 3)

# 3. distinct GO terms among the annotated components ----
# From the saved go_benchmark objects, which exist for three slides. Significance here is the
# NOMINAL weight01 Fisher p < 0.05 held in those tables, NOT the BH-on-classic basis behind the
# 13-vs-7 counts: term-level BH values were never saved, only the per-component minimum.
gb <- list.files(GO_DIR, pattern = "^go_benchmark__.*__k20\\.rds$", full.names = TRUE)
p3 <- do.call(rbind, lapply(gb, function(f) {
  sl <- sub("^go_benchmark__(.+)__k20\\.rds$", "\\1", basename(f))
  pm <- readRDS(f)$topN$per_method
  out <- do.call(rbind, lapply(names(pm), function(m) {
    gr <- pm[[m]]$go_results
    sig <- lapply(gr, function(r) if (is.null(r) || !nrow(r)) character(0)
                  else r$GO.ID[which(r$Fisher < 0.05)])
    ann <- sum(lengths(sig) > 0)
    data.frame(slide = sl, config = m,
               loadings = if (m %in% PROJ) "projected" else "gene-axis source",
               n_components = length(gr), n_annotated = ann,
               n_distinct_terms = length(unique(unlist(sig))),
               terms_per_annotated = if (ann) length(unique(unlist(sig))) / ann else NA_real_,
               stringsAsFactors = FALSE)
  }))
  rm(pm); gc(); out
}))
write.csv(p3, file.path(OUT_DIR, "cfs_mechanism__3_go_term_sharing.csv"), row.names = FALSE)
cat("\n=== 3. distinct GO terms among annotated components (nominal weight01 p < 0.05) ===\n")
print(p3, row.names = FALSE, digits = 3)

# 2. do the loadings track the gene's overall profile? ----
if (DO_POINT2) {
  slide_names <- names(SLIDES)
  stopifnot(length(slide_names) == 5L)

  source("scripts/benchmark/bench_B_setup.R")       # load_master + its gate against A
  cfg_all <- read.csv("config/runs.csv", stringsAsFactors = FALSE)

  p2 <- list()
  for (sl in slide_names) {
    cat("\n[", sl, "] loading normalised matrix\n", sep = "")
    M <- load_master(cfg_all[cfg_all$run_id == paste0(sl, "__k20__gp05"), ])
    X <- as.matrix(M$srt[["SCT"]]@scale.data)
    g_mean <- rowMeans(X); g_sd <- matrixStats::rowSds(X); names(g_sd) <- rownames(X)
    rm(M); gc()
    for (nm in names(cfg <- configs_of(sl))) {
      S <- cfg[[nm]]; g <- intersect(rownames(S), rownames(X))
      cm <- vapply(seq_len(ncol(S)), function(i) abs(cor(S[g, i], g_mean[g])), numeric(1))
      cs <- vapply(seq_len(ncol(S)), function(i) abs(cor(abs(S[g, i]), g_sd[g])), numeric(1))
      p2[[length(p2) + 1L]] <- data.frame(
        slide = sl, config = nm,
        loadings = if (nm %in% PROJ) "projected" else "gene-axis source",
        n_genes = length(g),
        mean_abs_cor_with_gene_mean = mean(cm),
        mean_abs_cor_with_gene_sd   = mean(cs), stringsAsFactors = FALSE)
    }
    rm(cfg, X, g_mean, g_sd); gc()
  }
  P2 <- do.call(rbind, p2)
  if (is.null(P2) || !nrow(P2))
    stop("point 2 produced no rows -- refusing to write an empty file")
  write.csv(P2, file.path(OUT_DIR, "cfs_mechanism__2_loading_vs_expression.csv"), row.names = FALSE)
  cat("\n=== 2. |cor| of loadings with the gene's overall profile ===\n")
  cat("    SCT residuals are near zero-mean per gene, so the mechanistically relevant\n",
      "   covariate is the row SD (the magnitude a projection picks up), reported alongside\n",
      "   the row mean.\n", sep = "")
  print(P2, row.names = FALSE, digits = 3)
} else cat("\n[2] skipped (DO_POINT2 = FALSE)\n")

# ---- verdict -----
agg <- aggregate(cbind(mean_pairwise_jaccard, union_ratio) ~ loadings, P1, mean)
cat("\n=== averaged over the five slides ===\n"); print(agg, row.names = FALSE, digits = 3)
cat("\nSummary of top-gene overlap, union coverage and GO-term diversity across configurations.\n")
cat("\nCFS_As_Is and icafast_no_Rmat are the same deterministic decomposition, so identical\n",
    "rows for those two are a correctness check, not a duplicate.\n", sep = "")
