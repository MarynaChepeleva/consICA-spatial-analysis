# Gene selection, topGO, the semantic-similarity metrics and the plots of the GO side of the
# ICA benchmark. Definitions only: sourcing this file runs nothing. GO_pipeline.R sources it
# and drives it over the slides; the control-arm scripts source it to reuse the same
# selection and enrichment without running the pipeline.
#
# Two layers:
#   (1) UNIFORM topGO on S_matrix for every method -> fair cross-method comparison
#   (2) NATIVE consICA::getGO on the consica_object -> consICA's own reporting
#
# Gene selection (S_matrix = genes x components), POSITIVE (signal) tail only:
#   components are oriented (skewness) so the heavy tail is positive, then
#   "topN"   primary      - top_n most-positive genes per component
#   "zscore" sensitivity  - genes with positive z-score of loading > z_thresh
# Positive tail only: pooling |loading| would mix opposing up/down programs.

suppressPackageStartupMessages({
  library(topGO)
  library(org.Hs.eg.db)
  library(consICA)
  library(GOSemSim)
  library(ggplot2)
  library(dplyr)
  library(purrr)
})

# 0. LOAD REPRESENTATIVE pure_stability RUNS FROM DISK ----

# Returns a named list: method -> list(S_matrix, consica_object).
# consica_object is NULL for every method except consICA.
load_pure_stability_runs <- function(output_dir, slide_id, methods,
                                     k, n_iter, rep_iter = 1L) {
  runs <- lapply(methods, function(m) {
    tag <- sprintf("%s__%s__pure_stability__k%d__n%d", slide_id, m, k, n_iter)
    path <- file.path(output_dir, paste0("runs__", tag, ".rds"))
    if (!file.exists(path)) {
      warning("Missing pure_stability file for method '", m, "': ", path)
      return(NULL)
    }
    obj <- readRDS(path)
    run <- obj$runs[[rep_iter]]
    list(S_matrix       = run$S_matrix,
         consica_object = run$consica_object)   # NULL unless consICA
  })
  setNames(runs, methods)
}

# 1. GENE SELECTION PER COMPONENT ----

# Orient a component so its heavy (signal) tail is positive: flip if skewness < 0.
# ICA sign is arbitrary; this makes the "positive tail" the informative program,
# uniformly across methods (approximates consICA's compsign).
.orient_positive <- function(w) {
  s <- mean(((w - mean(w)) / sd(w))^3)   # skewness
  if (is.finite(s) && s < 0) -w else w
}

# Positive tail only: pooling |loading| mixes opposing up/down programs and muddies
# GO. We orient, then take the positive side (top_n most positive, or z > z_thresh).
select_genes_per_component <- function(S_matrix, selection = "topN",
                                       top_n = 100L, z_thresh = 3, tail = "pos") {
  apply(S_matrix, 2, function(w) {
    wo <- if (tail == "pos") .orient_positive(w) else w
    if (selection == "topN") {
      names(sort(wo, decreasing = TRUE))[seq_len(min(top_n, length(wo)))]
    } else if (selection == "zscore") {
      z <- (wo - mean(wo)) / sd(wo)
      names(wo)[which(z > z_thresh)]      # positive tail only
    } else {
      stop("Unknown selection: ", selection)
    }
  }, simplify = FALSE)
}

# 2. topGO FOR ONE COMPONENT ----

# GenTable returns tiny p-values as strings like "< 1e-30" -> coerce safely.
.parse_fisher <- function(x) {
  suppressWarnings(as.numeric(gsub("[^0-9eE.+-]", "", as.character(x))))
}

run_topGO_one_component <- function(top_genes, all_genes, ontology = "BP",
                                    orgdb = "org.Hs.eg.db", id = "symbol",
                                    node_size = 10, top_n_terms = 20) {
  if (length(top_genes) == 0) return(NULL)

  gene_list <- factor(as.integer(all_genes %in% top_genes))
  names(gene_list) <- all_genes

  go_data <- suppressMessages(
    new("topGOdata",
        ontology = ontology,
        allGenes = gene_list,
        geneSel = function(x) x == 1,
        annot = annFUN.org,
        mapping = orgdb,
        ID = id,
        nodeSize  = node_size))

  # weight01 Fisher p is a NOMINAL per-term p. topGO weight01 p-values are not
  # independent and are not intended for multiple-testing correction, so this is
  # reported as nominal, NOT FDR; p.adjust is deliberately never applied.
  result <- suppressMessages(
    runTest(go_data, algorithm = "weight01", statistic = "fisher"))

  res_table <- GenTable(go_data, Fisher = result,
                        topNodes = top_n_terms, numChar = 100)
  res_table$Fisher <- .parse_fisher(res_table$Fisher)
  res_table$Enrichment_fold <- res_table$Significant / res_table$Expected
  res_table
}

# 3. topGO FOR ALL COMPONENTS OF ONE METHOD RUN ----

run_GO_for_method <- function(S_matrix, all_genes = NULL,
                              selection = "topN", top_n_genes = 100L,
                              z_thresh = 3, ontology = "BP",
                              orgdb = "org.Hs.eg.db", id = "symbol",
                              p_threshold = 0.05, verbose = TRUE) {
  if (is.null(all_genes)) all_genes <- rownames(S_matrix)

  comp_names <- colnames(S_matrix)
  gene_sets <- select_genes_per_component(S_matrix, selection = selection,
                                           top_n = top_n_genes, z_thresh = z_thresh)

  results <- list()
  for (comp in comp_names) {
    if (verbose) cat(sprintf("  GO [%s | %s]: %s (%d genes)\n",
                             selection, ontology, comp, length(gene_sets[[comp]])))
    results[[comp]] <- tryCatch({
      go_res <- run_topGO_one_component(
        top_genes = gene_sets[[comp]], all_genes = all_genes,
        ontology = ontology, orgdb = orgdb, id = id)
      if (!is.null(go_res)) {
        go_res$Component <- comp
        go_res$Significant_GO <- go_res$Fisher < p_threshold
      }
      go_res
    }, error = function(e) {
      if (verbose) cat(sprintf("    WARNING: GO failed for %s: %s\n", comp, e$message))
      NULL
    })
  }
  results
}

# 4. INTERPRETABILITY METRICS ----
#
# GO_Uniqueness is the biological analogue of the benchmark's Intra_Uniqueness.

compute_interpretability_metrics <- function(go_results_list, p_threshold = 0.05) {
  n_total <- length(go_results_list)

  n_significant <- sum(vapply(go_results_list, function(res) {
    if (is.null(res)) return(FALSE)
    any(res$Fisher < p_threshold, na.rm = TRUE)
  }, logical(1)))

  best_pvals <- vapply(go_results_list, function(res) {
    if (is.null(res)) return(NA_real_)
    min(res$Fisher, na.rm = TRUE)
  }, numeric(1))
  mean_neg_log_p <- mean(-log10(best_pvals), na.rm = TRUE)

  sig_terms_per_comp <- lapply(go_results_list, function(res) {
    if (is.null(res)) return(character(0))
    res$GO.ID[which(res$Fisher < p_threshold)]
  })

  go_uniqueness <- NA_real_
  if (n_total >= 2) {
    pairs <- combn(n_total, 2)
    jac <- apply(pairs, 2, function(p) {
      ti <- sig_terms_per_comp[[p[1]]]; tj <- sig_terms_per_comp[[p[2]]]
      if (length(ti) == 0 || length(tj) == 0) return(0)
      length(intersect(ti, tj)) / length(union(ti, tj))
    })
    go_uniqueness <- 1 - mean(jac, na.rm = TRUE)
  }

  data.frame(
    Frac_Interpretable = n_significant / n_total,
    N_Interpretable = n_significant,
    N_Total = n_total,
    Mean_neg_log10_p = mean_neg_log_p,
    GO_Uniqueness = go_uniqueness)   # 1 = every component enriches distinct terms
}

# 4b. SEMANTIC-SIMILARITY METRICS (GOSemSim) ----
# Exact-GO-ID Jaccard treats "epidermis development" and "skin development" as
# unrelated. Semantic similarity scores them as close, so these metrics reflect
# biological themes rather than exact labels. Build semData once (godata) and pass in.

# term set per component: top_n_terms strongest by p (Fisher for uniform tables)
.component_term_sets <- function(go_results_list, top_n_terms, p_col = "Fisher") {
  lapply(go_results_list, function(res) {
    if (is.null(res) || !nrow(res)) return(character(0))
    head(res$GO.ID[order(res[[p_col]])], top_n_terms)
  })
}

# Semantic uniqueness = 1 - mean pairwise mgoSim across components (higher = the
# components cover more distinct biological themes).
compute_semantic_uniqueness <- function(go_results_list, semData,
                                        top_n_terms = 10L,
                                        measure = "Wang", combine = "BMA") {
  ts <- .component_term_sets(go_results_list, top_n_terms)
  ts <- ts[lengths(ts) > 0]
  if (length(ts) < 2) return(NA_real_)
  pairs <- utils::combn(length(ts), 2)
  sims <- apply(pairs, 2, function(p) {
    s <- GOSemSim::mgoSim(ts[[p[1]]], ts[[p[2]]], semData = semData,
                          measure = measure, combine = combine)
    if (is.na(s)) 0 else s
  })
  1 - mean(sims)
}

# Convenience: semantic uniqueness per method for a run_GO_benchmark() result.
add_semantic_uniqueness <- function(go_benchmark, semData, top_n_terms = 10L) {
  vapply(names(go_benchmark$per_method), function(m)
    compute_semantic_uniqueness(go_benchmark$per_method[[m]]$go_results,
                                semData, top_n_terms = top_n_terms),
    numeric(1))
}

# Semantic concordance (consICA): mgoSim between uniform and native getGO term sets,
# per component (1:1 by index). ICA sign is arbitrary and the uniform (skewness) and
# native (compsign) orientations may disagree, so we match against whichever native
# tail (pos/neg) agrees better -> orientation-robust. Replaces exact-ID Jaccard.
compute_go_concordance_semantic <- function(uniform_go_results, native_go_bp,
                                            semData, top_n_terms = 10L,
                                            rank_col = "classicFisher",
                                            measure = "Wang", combine = "BMA") {
  .topset <- function(tbl, p_col) {
    if (is.null(tbl) || !nrow(tbl)) return(character(0))
    head(tbl$GO.ID[order(tbl[[p_col]])], top_n_terms)
  }
  .sim <- function(a, b) {
    if (length(a) == 0 || length(b) == 0) return(NA_real_)
    GOSemSim::mgoSim(a, b, semData = semData, measure = measure, combine = combine)
  }
  n <- min(length(uniform_go_results), length(native_go_bp))
  sim <- vapply(seq_len(n), function(i) {
    u_ids <- .topset(uniform_go_results[[i]], "Fisher")
    vals <- c(.sim(u_ids, .topset(native_go_bp[[i]]$pos, rank_col)),
               .sim(u_ids, .topset(native_go_bp[[i]]$neg, rank_col)))
    vals <- vals[!is.na(vals)]
    if (!length(vals)) NA_real_ else max(vals)   # best-matching tail
  }, numeric(1))
  data.frame(Component = seq_len(n), GO_Concordance_semantic = sim)
}

# 5. UNIFORM GO BENCHMARK ACROSS METHODS ----

run_GO_benchmark <- function(method_runs_list, selection = "topN",
                             top_n_genes = 100L, z_thresh = 3,
                             ontology = "BP", orgdb = "org.Hs.eg.db",
                             id = "symbol", p_threshold = 0.05,
                             metric_fn = compute_interpretability_metrics) {
  results_summary <- list()

  for (method_name in names(method_runs_list)) {
    run <- method_runs_list[[method_name]]
    if (is.null(run)) next
    cat(sprintf("\n=== GO Enrichment: %s (%s) ===\n", method_name, selection))

    go_results <- run_GO_for_method(
      S_matrix = run$S_matrix, selection = selection,
      top_n_genes = top_n_genes, z_thresh = z_thresh,
      ontology = ontology, orgdb = orgdb, id = id,
      p_threshold = p_threshold, verbose = TRUE)

    metrics <- metric_fn(go_results, p_threshold)
    metrics$Method <- method_name
    metrics$Selection <- selection

    results_summary[[method_name]] <- list(go_results = go_results, metrics = metrics)
  }

  metrics_df <- do.call(rbind, lapply(results_summary, `[[`, "metrics"))
  rownames(metrics_df) <- NULL
  list(per_method = results_summary, summary = metrics_df)
}

# 6. NATIVE consICA getGO  (consICA section only; cross-validated) ----
run_consica_native_GO <- function(consica_object, genome = "org.Hs.eg.db",
                                  db = "BP", ncores = 1L, ...) {
  # getGO signature: getGO(cica, alpha, genenames, genome, db, ncores, rotate)
  #   genome = annotation package ("org.Hs.eg.db" / "org.Mm.eg.db")
  #   db     = ontology ("BP" / "CC" / "MF")
  if (is.null(consica_object)) {
    warning("No consica_object supplied — skipping native getGO.")
    return(NULL)
  }
  consICA::getGO(consica_object, genome = genome, db = db, ncores = ncores, ...)
}

# 7. VISUALISATION ----

plot_GO_benchmark <- function(go_benchmark_result) {
  library(ggplot2)
  df <- go_benchmark_result$summary

  mk_bar <- function(yvar, title, ylab) {
    ggplot(df, aes(x = reorder(Method, .data[[yvar]]),
                   y = .data[[yvar]], fill = Method)) +
      geom_col() + coord_flip() +
      labs(title = title, x = NULL, y = ylab) +
      theme_minimal() + theme(legend.position = "none")
  }

  list(
    p1 = mk_bar("Frac_Interpretable", "Fraction of interpretable components",
                "Fraction (>=1 GO term, nominal weight01 p < 0.05)"),
    p2 = mk_bar("GO_Uniqueness", "GO term uniqueness across components",
                "1 - mean Jaccard of significant terms"),
    p3 = mk_bar("Mean_neg_log10_p", "Mean enrichment strength per component",
                "mean -log10(best p-value)"))
}

# Spatial intensity map of one consICA component (M row over spots)
plot_component_map <- function(cons_obj, comp, sce, slide_id = "") {
  act <- cons_obj$M[comp, ]                       # activity per spot (components x spots)

  cd <- as.data.frame(SummarizedExperiment::colData(sce))
  if ("spatial" %in% SingleCellExperiment::reducedDimNames(sce)) {
    xy <- as.data.frame(SingleCellExperiment::reducedDim(sce, "spatial"))[, 1:2]
  } else {
    xc <- intersect(c("pxl_col_in_fullres", "array_col", "imagecol"), colnames(cd))[1]
    yc <- intersect(c("pxl_row_in_fullres", "array_row", "imagerow"), colnames(cd))[1]
    xy <- cd[, c(xc, yc)]
  }
  colnames(xy) <- c("x", "y")
  rownames(xy) <- colnames(sce)

  common <- intersect(names(act), rownames(xy))    # match spots by barcode
  df <- data.frame(xy[common, ], intensity = act[common])

  ggplot(df, aes(x, y, colour = intensity)) +
    geom_point(size = 1.8) +
    scale_colour_gradient2(low = "#2166ac", mid = "grey92",
                           high = "#b2182b", midpoint = 0) +
    coord_fixed() + scale_y_reverse() + theme_void() +
    ggtitle(sprintf("%s — component %s activity", slide_id, comp))
}
