# READ-ONLY probe on what Intra_Uniqueness measures. The single-run comparators feed t(X) to
# icafast (observations = spots) and derive the gene signature by projection X %*% res$S,
# while consICA runs ICA with genes in rows and takes S directly; both wrappers are in
# bench_ica_functions.R. ICA whitens then applies an orthogonal rotation, so the returned S
# columns are uncorrelated across the observation axis by construction, and a single run with
# genes in rows sits at Intra_Uniqueness ~= 1.0. Table 2 of the paper carries this as a
# footnote.
#
# Nothing is re-preprocessed: X is the exact matrix the benchmark saw (consica_object$X_num).
# Does NOT touch bench_ica_perturbations.R, GO_pipeline.R, or any saved runs__*.rds.

suppressPackageStartupMessages({
  library(ica)
  library(fastICA)
})

run_dir <- "output/ica_perturbations"
out_dir <- "output/benchmark_summary"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_csv <- file.path(out_dir, "intra_uniqueness_axis_check.csv")

K <- 20L
SEEDS <- 1:5 # for the single-run objects: a value fixed by construction must not wander
MAXIT <- 600L # as the benchmark wrappers
SLIDES <- c("MEND40", "TENX117")  # axes differ by two orders of magnitude in n_obs

# EXACT benchmark metric, compute_internal_uniqueness() in bench_ica_functions.R: mean over
# diagonal zeroed and NA -> 0. mean_abs_cor is that mean; intra_uniqueness = 1 - it.
iu <- function(S) {
  cm <- abs(cor(S)); cm[is.na(cm)] <- 0; diag(cm) <- 0
  c(mean_abs_cor = mean(cm), intra_uniqueness = 1 - mean(cm))
}

one <- function(slide, object, axis, n_obs, S) {
  m <- iu(S)
  data.frame(Slide = slide, Object = object, axis = axis, n_obs = n_obs,
             mean_abs_cor = m[["mean_abs_cor"]], mean_abs_cor_sd = NA_real_,
             intra_uniqueness = m[["intra_uniqueness"]], intra_uniqueness_sd = NA_real_,
             stringsAsFactors = FALSE)
}

# single-run objects across SEEDS -> mean + SD (a by-construction value should have ~0 SD)
seeded <- function(slide, object, axis, n_obs, make_S) {
  v <- vapply(SEEDS, function(s) { set.seed(s); iu(make_S()) }, numeric(2))  # rows: mean_abs_cor, intra
  data.frame(Slide = slide, Object = object, axis = axis, n_obs = n_obs,
             mean_abs_cor = mean(v[1, ]), mean_abs_cor_sd = sd(v[1, ]),
             intra_uniqueness = mean(v[2, ]), intra_uniqueness_sd = sd(v[2, ]),
             stringsAsFactors = FALSE)
}

rows <- list()
for (slide in SLIDES) {
  path <- file.path(run_dir, sprintf("runs__%s__consICA__pure_stability__k20__n10.rds", slide))
  r1 <- readRDS(path)$runs[[1]]
  X <- r1$consica_object$X_num                      # genes x spots, exactly as the benchmark saw it
  if (is.null(X)) stop("no consica_object$X_num in ", path)
  ng <- nrow(X); ns <- ncol(X)
  cat(sprintf("%-8s X = %d genes x %d spots\n", slide, ng, ns))

  # ---- GENE axis (n_obs = genes): does a single run already sit at ~1.0? ----
  rows[[length(rows) + 1L]] <- seeded(slide, "1_single_icafast_geneaxis_Sdirect", "gene", ng,
                                      function() ica::icafast(X, nc = K, maxit = MAXIT)$S)
  rows[[length(rows) + 1L]] <- seeded(slide, "2_single_fastICA_geneaxis_Sdirect", "gene", ng,
                                      function() fastICA::fastICA(X, n.comp = K, maxit = MAXIT)$S)
  rows[[length(rows) + 1L]] <- one(slide, "3_consICA_S_matrix_geneaxis", "gene", ng, r1$S_matrix)          # expect ~0.984
  rows[[length(rows) + 1L]] <- one(slide, "4_comparator_Xprojected_geneaxis", "gene", ng,
                                   X %*% ica::icafast(t(X), nc = K, maxit = MAXIT)$S)                       # expect ~0.764

  # ---- SPOT axis (n_obs = spots): the mirror image -----
  rows[[length(rows) + 1L]] <- one(slide, "5_consICA_M_matrix_spotaxis", "spot", ns, t(r1$M_matrix))
  rows[[length(rows) + 1L]] <- one(slide, "6_comparator_icafast_spotaxis", "spot", ns,
                                   ica::icafast(t(X), nc = K, maxit = MAXIT)$S)
}

out <- do.call(rbind, rows)
out <- out[order(out$Slide, out$Object), ]
write.csv(out, out_csv, row.names = FALSE)

cat("\nR:", R.version.string,
    "| ica:", as.character(packageVersion("ica")),
    "| fastICA:", as.character(packageVersion("fastICA")), "\n\n")
print(out, row.names = FALSE, digits = 4)
cat("\nRead as:\n",
    " objects 1-2 ~0.98-1.0 with ~0 SD -> metric tracks the axis, not consensus (0.984 vs 0.764 is axis, not gain)\n",
    " objects 1-2 ~0.75-0.85          -> consensus really adds, and only the gene-axis single-run comparator was missing\n",
    " objects 5-6 show the mirror on the spot axis\n", sep = "")
cat("\nWrote", normalizePath(out_csv, mustWork = FALSE), "\n")
