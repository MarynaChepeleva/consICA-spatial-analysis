# Wrappers, metrics and the benchmark engine. Definitions only: sourcing this file runs
# nothing. bench_ica_perturbations.R sources it and drives it over config/runs.csv; the
# control-arm scripts source it to reuse the same wrappers without running the benchmark.
#
# Conventions (do not break):
#   S_matrix = genes x components   (feature loadings)
#   M_matrix = components x spots    (mixing / spot activity)
#   consICA:  cons_res$M is ALREADY components x spots -> DO NOT transpose on save.
# Every wrapper returns the same structure: list(S_matrix, M_matrix).
# consICA additionally returns consica_object (full object for getGO()).
#
# consICA cores: ncores = 1 on pure_stability (those decompositions feed GO),
#                ncores = consica_ncores (default 4) on perturbation modes.
#
# Metrics (all derived from S_matrix, identically for every method):
#   Mean_Stability     - mean over run pairs of mean(rowMaxs(abs(cor(S_i,S_j))))
#   Matched_Stability  - same but with optimal 1:1 component assignment (Hungarian)
#   Intra_Uniqueness   - 1 - mean(abs(cor(S))) with zero diagonal, averaged over runs
#   N_Repro_<thr>      - mean count of components with matched corr > threshold
#   TopGene_Jaccard    - mean Jaccard of top-N gene sets of matched components
#   *_Runtime_sec      - wall-clock cost of the decomposition (SCTransform excluded)


invisible(gc())

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(zellkonverter)
  library(Seurat)
  library(matrixStats)
  library(ggplot2)
  library(ica)
  library(fastICA)
  library(purrr)
  library(consICA)
  library(e1071)
  library(clue)
})

# Allow large Visium matrices through SCTransform's future workers
options(future.globals.maxSize = 4 * 1024^3)

# 1. WRAPPERS  (unified interface: (X, srt_obj, k, seed, ncores) -> list(S_matrix, M_matrix[, consica_object])) ----
#    X        : genes x spots  (SCT scale.data)
#    S_matrix : genes x components
#    M_matrix : components x spots
#    ncores   : used only by consICA; ignored by the other wrappers

.name_components <- function(S_matrix, M_matrix) {
  comp <- paste0("IC_", seq_len(ncol(S_matrix)))
  colnames(S_matrix) <- comp
  rownames(M_matrix) <- comp
  list(S_matrix = S_matrix, M_matrix = M_matrix)
}

wrapper_icafast_det <- function(X, srt_obj, k, seed, ncores = 1L) {
  res <- ica::icafast(X = t(X), nc = k, maxit = 600)   # res$S: spots x comp
  S_matrix <- X %*% res$S                                   # genes x comp
  M_matrix <- t(res$S)                                      # comp x spots
  rownames(S_matrix) <- rownames(X)
  colnames(M_matrix) <- colnames(X)
  .name_components(S_matrix, M_matrix)
}

wrapper_icafast_stoch <- function(X, srt_obj, k, seed, ncores = 1L) {
  set.seed(seed)
  rand_mat <- matrix(rnorm(k * k), nrow = k, ncol = k)
  rsvd <- svd(rand_mat)
  Rmat_random <- tcrossprod(rsvd$u, rsvd$v)
  res <- ica::icafast(X = t(X), nc = k, Rmat = Rmat_random, maxit = 600)
  S_matrix <- X %*% res$S
  M_matrix <- t(res$S)
  rownames(S_matrix) <- rownames(X)
  colnames(M_matrix) <- colnames(X)
  .name_components(S_matrix, M_matrix)
}

wrapper_cfs_native <- function(X, srt_obj, k, seed, ncores = 1L) {
  srt_res <- custom_ICASpatial(srt_obj, nics = k, seed.use = seed)
  S_matrix <- srt_res@reductions$ica@feature.loadings       # genes x comp
  M_matrix <- t(srt_res@reductions$ica@cell.embeddings)     # comp x spots
  list(S_matrix = S_matrix, M_matrix = M_matrix)
}

wrapper_fastica_pkg <- function(X, srt_obj, k, seed, ncores = 1L) {
  set.seed(seed)
  res <- fastICA::fastICA(X = t(X), n.comp = k, maxit = 600)  # res$S: spots x comp, res$A: comp x genes
  S_matrix <- t(res$A)                                             # genes x comp
  M_matrix <- t(res$S)                                             # comp x spots
  rownames(S_matrix) <- rownames(X)
  colnames(M_matrix) <- colnames(X)
  .name_components(S_matrix, M_matrix)
}

wrapper_consica_pipeline <- function(X, srt_obj, k, seed, ncores = 4L) {
  set.seed(seed)
  cons_res <- consICA::consICA(X = X, ncomp = k, ntry = 20,
                               ncores = ncores, verbose = FALSE)
  S_matrix <- cons_res$S            # genes x components
  M_matrix <- cons_res$M            # components x spots  -- DO NOT transpose
  rownames(S_matrix) <- rownames(X)
  colnames(M_matrix) <- colnames(X)
  list(S_matrix = S_matrix,
       M_matrix = M_matrix,
       consica_object = cons_res) # full object for getGO()
}

# CONTROL B: identical to wrapper_consica_pipeline except ntry = 1 —
# one ICA run, no consensus assembly. Same core, same orientation, same seed handling, so
# A - B isolates the consensus procedure. With a single try there is nothing to parallelise,
# hence ncores = 1. Used only by bench_B_run.R; not in the driver's `methods` vector.
wrapper_consica_single <- function(X, srt_obj, k, seed, ncores = 1L) {
  set.seed(seed)
  cons_res <- consICA::consICA(X = X, ncomp = k, ntry = 1,
                               ncores = ncores, verbose = FALSE)
  S_matrix <- cons_res$S            # genes x components
  M_matrix <- cons_res$M            # components x spots  -- DO NOT transpose
  rownames(S_matrix) <- rownames(X)
  colnames(M_matrix) <- colnames(X)
  list(S_matrix = S_matrix,
       M_matrix = M_matrix,
       consica_object = cons_res)
}

# 2. CFS BLUEPRINT ----
# custom_ICASpatial is a copy of CellsFromSpace's ICASpatial and carries that
# package's GPL (>= 2). It sits in its own file so the licence boundary is a
# file boundary.
source("scripts/benchmark/cfs_ICASpatial.R")

# 3. METRICS ----

# abs correlation between two S matrices on their shared genes -> comp_i x comp_j
.pair_abs_cor <- function(Si, Sj) {
  common <- intersect(rownames(Si), rownames(Sj))
  cm <- abs(cor(Si[common, , drop = FALSE], Sj[common, , drop = FALSE]))
  cm[is.na(cm)] <- 0 # constant component -> treat as no similarity
  cm
}

# top-N genes per component (by |loading|), returns list of length ncol
.top_genes <- function(S, top_n) {
  n <- min(top_n, nrow(S))
  lapply(seq_len(ncol(S)), function(ci)
    rownames(S)[order(abs(S[, ci]), decreasing = TRUE)[seq_len(n)]])
}

# Mean_Stability: rowMaxs matching (many-to-one allowed). Returns per-pair vector.
compute_cross_stability <- function(run_list) {
  n <- length(run_list)
  if (n < 2) return(list(mean = 1, sd = 0, scores = 1))
  pairs <- utils::combn(n, 2)
  scores <- apply(pairs, 2, function(p) {
    cm <- .pair_abs_cor(run_list[[p[1]]]$S_matrix, run_list[[p[2]]]$S_matrix)
    mean(matrixStats::rowMaxs(cm, na.rm = TRUE))
  })
  list(mean = mean(scores), sd = sd(scores), scores = scores)
}

# Matched_Stability: optimal 1:1 assignment (Hungarian). Also counts reproducible
# components (matched corr > threshold). Returns per-pair vectors for plotting.
compute_matched_stability <- function(run_list, thresholds = c(0.5, 0.7)) {
  n <- length(run_list)
  empty <- setNames(rep(NA_real_, length(thresholds)), thresholds)
  if (n < 2) return(list(mean = 1, sd = 0, scores = 1,
                         matched_by_pair = list(), n_repro = empty))
  pairs <- utils::combn(n, 2)
  per_pair <- rep(NA_real_, ncol(pairs))
  matched_list <- vector("list", ncol(pairs))
  repro_counts <- matrix(NA_real_, nrow = ncol(pairs), ncol = length(thresholds))
  for (idx in seq_len(ncol(pairs))) {
    Si <- run_list[[pairs[1, idx]]]$S_matrix
    Sj <- run_list[[pairs[2, idx]]]$S_matrix
    if (ncol(Si) != ncol(Sj)) next            # unequal component count -> skip pair
    cm <- .pair_abs_cor(Si, Sj)
    assign <- as.integer(clue::solve_LSAP(cm, maximum = TRUE))
    matched <- cm[cbind(seq_len(nrow(cm)), assign)]
    matched_list[[idx]] <- matched
    per_pair[idx] <- mean(matched)
    repro_counts[idx, ] <- vapply(thresholds,
                                  function(th) sum(matched > th), numeric(1))
  }
  list(mean = mean(per_pair, na.rm = TRUE),
       sd = sd(per_pair, na.rm = TRUE),
       scores = per_pair,
       matched_by_pair = matched_list,
       n_repro = setNames(colMeans(repro_counts, na.rm = TRUE), thresholds))
}

# TopGene_Jaccard: Jaccard of top-N gene sets of matched (1:1) components.
compute_topgene_stability <- function(run_list, top_n = 50) {
  n <- length(run_list)
  if (n < 2) return(list(mean = 1, sd = 0, scores = 1))
  pairs <- utils::combn(n, 2)
  scores <- apply(pairs, 2, function(p) {
    Si <- run_list[[p[1]]]$S_matrix
    Sj <- run_list[[p[2]]]$S_matrix
    if (ncol(Si) != ncol(Sj)) return(NA_real_)
    cm <- .pair_abs_cor(Si, Sj)
    assign <- as.integer(clue::solve_LSAP(cm, maximum = TRUE))
    topi <- .top_genes(Si, top_n)
    topj <- .top_genes(Sj, top_n)
    jac <- vapply(seq_along(topi), function(ci) {
      gi <- topi[[ci]]; gj <- topj[[assign[ci]]]
      length(intersect(gi, gj)) / length(union(gi, gj))
    }, numeric(1))
    mean(jac)
  })
  list(mean = mean(scores, na.rm = TRUE),
       sd = sd(scores, na.rm = TRUE),
       scores = scores)
}

# Intra_Uniqueness: 1 - mean(abs(cor(S))) with zero diagonal, averaged over runs.
compute_internal_uniqueness <- function(run_list) {
  scores <- purrr::map_dbl(run_list, function(res) {
    corr_matrix <- abs(cor(res$S_matrix))
    corr_matrix[is.na(corr_matrix)] <- 0
    diag(corr_matrix) <- 0
    1 - mean(corr_matrix)
  })
  mean(scores)
}

# 4. METADATA BUILDER  (identical fields for all methods) ----

make_run_metadata <- function(slide_id, test_file, method, mode,
                              k, n_iter, gene_perturb_frac, spot_perturb_frac,
                              ncores, top_n, repro_thresholds,
                              min_genes, min_counts, srt_master) {
  list(
    run_id = paste0(slide_id, "__", method, "__", mode, "__",
                               format(Sys.time(), "%Y%m%d_%H%M%S")),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    slide_id = slide_id,
    input_file = test_file,
    n_spots = ncol(srt_master), # master dims; per-run effective dims stored in runs
    n_genes = nrow(srt_master),
    assay_used = "SCT",
    normalization = "SCTransform (return.only.var.genes = FALSE)",
    spot_qc = list(min_genes = min_genes, min_counts = min_counts),
    method = method,
    mode = mode,
    k_components = k,
    n_iterations = n_iter,
    spot_perturb_frac = if (mode == "spot_perturbation") spot_perturb_frac else NA,
    gene_perturb_frac = if (mode == "gene_perturbation") gene_perturb_frac else NA,
    gene_perturb_pct = if (mode == "gene_perturbation") gene_perturb_frac * 100 else NA,
    benchmark_params = list(top_n_genes      = top_n,
                             repro_thresholds = repro_thresholds),
    method_params = switch(method,
                               "consICA" = list(
                                 ntry      = 20,
                                 ncores    = ncores,
                                 algorithm = "consensus ICA (consICA::consICA)"),
                               "consICA_ntry1" = list(
                                 ntry      = 1,
                                 ncores    = ncores,
                                 algorithm = "single-run ICA (consICA core, ntry = 1)"),
                               "CFS_As_Is" = list(
                                 seed      = 42,
                                 Rmat      = "diagonal (default icafast behaviour)",
                                 algorithm = "icafast deterministic - replicates CFS published implementation"),
                               "icafast_no_Rmat" = list(
                                 Rmat      = "diagonal (default)",
                                 algorithm = "ica::icafast deterministic"),
                               "icafast_with_Rmat" = list(
                                 Rmat      = "random orthogonal matrix",
                                 algorithm = "ica::icafast stochastic"),
                               "fastICA_Package" = list(
                                 algorithm = "fastICA::fastICA CRAN package, stochastic")),
    r_version = R.version$version.string,
    package_versions = list(
      consICA = as.character(packageVersion("consICA")),
      ica = as.character(packageVersion("ica")),
      Seurat = as.character(packageVersion("Seurat")),
      fastICA = as.character(packageVersion("fastICA")),
      clue = as.character(packageVersion("clue"))),
    seeds_used = seq_len(n_iter)
  )
}

# 5. Dynamic filename ----
#    runs__SLIDEID__METHOD__MODE__kK__nN[__gpXpct].rds

make_save_tag <- function(slide_id, method, mode, k, n_iter,
                          gene_perturb_pct = NA) {
  gp_str <- if (!is.na(gene_perturb_pct)) {
    paste0("_gp", gene_perturb_pct, "pct")
  } else ""
  sprintf("%s__%s__%s__k%d__n%d%s",
          slide_id, method, mode, k, n_iter, gp_str)
}

# 6. INPUT PREP PER ITERATION  (perturbation happens here; ICA timed separately) ----

prep_iteration_input <- function(srt_raw, mode, spot_perturb_frac, gene_perturb_frac) {
  if (mode == "pure_stability") {
    X <- as.matrix(srt_raw[["SCT"]]@scale.data)
    return(list(X = X, srt = srt_raw, subset = NA))

  } else if (mode == "spot_perturbation") {
    all_cells <- colnames(srt_raw)
    keep <- sample(all_cells, size = floor((1 - spot_perturb_frac) * length(all_cells)))
    srt_sub <- srt_raw[, keep]
    srt_sub <- suppressWarnings(SCTransform(
      srt_sub, assay = "originalexp",
      return.only.var.genes = FALSE, verbose = FALSE)) # re-normalise on subset
    X <- as.matrix(srt_sub[["SCT"]]@scale.data)
    return(list(X = X, srt = srt_sub, subset = keep))

  } else if (mode == "gene_perturbation") {
    all_genes <- rownames(srt_raw[["SCT"]])
    keep <- sample(all_genes, size = floor((1 - gene_perturb_frac) * length(all_genes)))
    srt_tmp <- srt_raw[keep, ] # subset SCT residuals (no re-fit)
    X <- as.matrix(srt_tmp[["SCT"]]@scale.data)
    return(list(X = X, srt = srt_tmp, subset = keep))

  } else {
    stop("Unknown mode: ", mode)
  }
}

# 7. MASTER BENCHMARK ENGINE  ->  list(metrics, runs, diagnostics) ----

execute_benchmark_mode <- function(srt_raw, mode, k, n_iter, method,
                                   gene_perturb_frac = 0.05,
                                   spot_perturb_frac = 0.05,
                                   ncores = 1L,
                                   top_n = 50L,
                                   repro_thresholds  = c(0.5, 0.7),
                                   keep_all_consica  = TRUE) {
  seeds <- seq_len(n_iter)
  run_algo <- switch(method,
                       "icafast_no_Rmat" = wrapper_icafast_det,
                       "icafast_with_Rmat" = wrapper_icafast_stoch,
                       "CFS_As_Is" = wrapper_cfs_native,
                       "fastICA_Package" = wrapper_fastica_pkg,
                       "consICA" = wrapper_consica_pipeline,
                       "consICA_ntry1" = wrapper_consica_single,
                       stop("Unknown method: ", method))

  results <- lapply(seeds, function(s) {
    set.seed(s)
    prep <- prep_iteration_input(srt_raw, mode, spot_perturb_frac, gene_perturb_frac)
    tm <- system.time(
      run <- run_algo(prep$X, prep$srt, k, s, ncores)) # ICA cost only (SCTransform excluded)
    run$elapsed_sec <- unname(tm[["elapsed"]])
    run$subset <- prep$subset # kept spots (spot mode) / genes (gene mode) / NA
    run$n_genes_eff <- nrow(prep$X)
    run$n_spots_eff <- ncol(prep$X)
    # consica_object stores dense copies of X (~1.4 GB on Visium). On large slides
    # keep only the GO representative (pure_stability, iter 1) to avoid OOM; on small
    # slides keep all (getGO flexibility). Controlled by keep_all_consica.
    if (!keep_all_consica && !(mode == "pure_stability" && s == seeds[1]))
      run$consica_object <- NULL
    run
  })

  cs <- compute_cross_stability(results)
  ms <- compute_matched_stability(results, thresholds = repro_thresholds)
  tg <- compute_topgene_stability(results, top_n = top_n)
  uq <- compute_internal_uniqueness(results)
  rt <- vapply(results, function(r) r$elapsed_sec, numeric(1))

  base_df <- data.frame(
    Method = method,
    Mode = mode,
    Gene_Perturb_Pct  = if (mode == "gene_perturbation") gene_perturb_frac * 100 else NA,
    Mean_Stability = cs$mean,
    SD_Stability = cs$sd,
    Matched_Stability = ms$mean,
    SD_Matched = ms$sd,
    Intra_Uniqueness = uq,
    TopGene_Jaccard = tg$mean,
    SD_TopGene = tg$sd,
    Mean_Runtime_sec = mean(rt),
    SD_Runtime_sec = sd(rt),
    Total_Runtime_sec = sum(rt),
    check.names = FALSE)
  repro_df <- as.data.frame(as.list(setNames(
    ms$n_repro, paste0("N_Repro_", repro_thresholds * 100))),
    check.names = FALSE)
  metrics_df <- cbind(base_df, repro_df)

  diagnostics <- list(
    cross_stability_scores = cs$scores,       # per-pair rowMaxs means
    matched_scores = ms$scores,       # per-pair 1:1 matched means
    matched_cor_by_pair = ms$matched_by_pair,
    n_repro_by_threshold = ms$n_repro,
    topgene_jaccard_scores = tg$scores,
    runtime_sec = rt,
    top_n_genes = top_n,
    repro_thresholds = repro_thresholds)

  list(metrics = metrics_df, runs = results, diagnostics = diagnostics)
}

# QC: is the slide already restricted to in-tissue spots? ----
check_in_tissue <- function(sce, slide_id) {
  cd <- SummarizedExperiment::colData(sce)
  it <- grep("in.?tissue", colnames(cd), ignore.case = TRUE, value = TRUE)
  n <- ncol(sce)
  if (!length(it)) {
    cat(sprintf("[%s] %d spots | no in_tissue column\n", slide_id, n)); return(invisible())
  }
  v <- cd[[it[1]]]
  ins <- sum(v %in% c(TRUE, 1, "1", "TRUE"))
  off <- n - ins
  verdict <- if (off == 0) "built-in column in_tissue"
  else sprintf("BACKGROUND PRESENT -> remove %d off-tissue spots", off)
  cat(sprintf("[%s] %d spots | in-tissue %d / off-tissue %d | %s\n",
              slide_id, n, ins, off, verdict))
}

# SLIDE OVERVIEW: Spatial map of spots coloured by total counts ----
plot_spatial_spots <- function(sce, slide_id) {
  cd <- as.data.frame(SummarizedExperiment::colData(sce))
  if ("spatial" %in% SingleCellExperiment::reducedDimNames(sce)) {
    xy <- as.data.frame(SingleCellExperiment::reducedDim(sce, "spatial"))[, 1:2]
  } else {
    xc <- intersect(c("pxl_col_in_fullres", "array_col", "imagecol"), colnames(cd))[1]
    yc <- intersect(c("pxl_row_in_fullres", "array_row", "imagerow"), colnames(cd))[1]
    xy <- cd[, c(xc, yc)]
  }
  colnames(xy) <- c("x", "y")
  xy$total_counts <- Matrix::colSums(SummarizedExperiment::assay(sce, "X"))
  ggplot(xy, aes(x, y, colour = total_counts)) +
    geom_point(size = 2) +
    scale_colour_viridis_c(trans = "log10") +
    coord_fixed() + scale_y_reverse() +
    theme_void() +
    ggtitle(sprintf("%s — spots by total counts (n=%d)", slide_id, ncol(sce)))
}

# QC: which spots the per-spot threshold keeps (tissue) vs removes (background)
plot_spot_qc <- function(sce, keep, slide_id) {
  cd <- as.data.frame(SummarizedExperiment::colData(sce))
  if ("spatial" %in% SingleCellExperiment::reducedDimNames(sce)) {
    xy <- as.data.frame(SingleCellExperiment::reducedDim(sce, "spatial"))[, 1:2]
  } else {
    xc <- intersect(c("pxl_col_in_fullres", "array_col", "imagecol"), colnames(cd))[1]
    yc <- intersect(c("pxl_row_in_fullres", "array_row", "imagerow"), colnames(cd))[1]
    xy <- cd[, c(xc, yc)]
  }
  colnames(xy) <- c("x", "y")
  xy$spot <- ifelse(keep, "tissue", "removed")
  ggplot(xy, aes(x, y, colour = spot)) +
    geom_point(size = 2) +
    scale_colour_manual(values = c(tissue = "#2c7fb8", removed = "grey80")) +
    coord_fixed() + scale_y_reverse() + theme_void() +
    ggtitle(sprintf("%s — QC: %d kept / %d removed", slide_id, sum(keep), sum(!keep)))
}

# Auto min_genes threshold: valley between background and tissue in log10(genes/spot).
# Guarded: if not clearly bimodal (no background), returns floor (no aggressive cut).
auto_spot_threshold <- function(counts, floor = 50L) {
  set.seed(1)                                  # reproducible threshold (kmeans is random)
  ng <- Matrix::colSums(counts > 0)
  lg <- log10(ng[ng > 0] + 1)
  km <- stats::kmeans(lg, centers = 2, nstart = 10)
  if (abs(diff(km$centers)) < 0.5) return(floor)                 # unimodal -> floor
  max(floor, round(10^max(lg[km$cluster == which.min(km$centers)]) - 1))
}
