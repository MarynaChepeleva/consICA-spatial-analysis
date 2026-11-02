# Driver of the benchmark. One row of config/runs.csv is one cell of the grid: a slide at
# one K with one gene-dropout fraction. A row runs four methods across three perturbation
# modes, n_iter decompositions each, and writes runs__TAG.rds per method and mode, then
# metrics_FINAL__*.csv and session_log__*.rds for the row itself. Everything lands in
# output/ica_perturbations/.

# Wrappers, metrics and the engine are in bench_ica_functions.R.

source("scripts/benchmark/bench_ica_functions.R")

# ---- configuration ----

df_params <- read.csv("config/runs.csv", stringsAsFactors = FALSE)

# All rows, as config.R writes them. To run a subset, list the ids here instead:
#   "TENX115__k20__gp05", "TENX117__k20__gp05",
#   "MEND40__k20__gp05",  "MEND40__k50__gp05",  "MEND40__k20__gp15",
#   "NCBI776__k20__gp05", "NCBI776__k50__gp05",
#   "NCBI600__k20__gp05", "NCBI600__k20__gp15"
RUN_IDS <- df_params$run_id

TOP_N_GENES <- 50L
REPRO_THRESHOLDS <- c(0.5, 0.7)
methods <- c("consICA", "CFS_As_Is", "icafast_no_Rmat", "icafast_with_Rmat")
modes <- c("pure_stability", "spot_perturbation", "gene_perturbation")

output_dir <- "output/ica_perturbations/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- one cell of the grid per row ----

for (RUN_ID in RUN_IDS) {

  cat(sprintf("\n======== %s (%d/%d) ========\n",
              RUN_ID, match(RUN_ID, RUN_IDS), length(RUN_IDS)))

  cfg_run <- df_params[df_params$run_id == RUN_ID, ]
  stopifnot(nrow(cfg_run) == 1L)

  slide_id <- cfg_run$slide
  test_file <- cfg_run$file

  # A row whose slide is not on disk is skipped rather than allowed to fail here.
  if (!file.exists(test_file)) {
    message(RUN_ID, ": ", test_file, " not found, skipped")
    next
  }

  k_components <- as.integer(cfg_run$k)
  min_genes <- if (is.na(cfg_run$min_genes)) NA_integer_ else as.integer(cfg_run$min_genes)
  min_counts <- 0L
  consica_ncores <- as.integer(cfg_run$ncores)
  keep_all_consica <- as.logical(cfg_run$store_consica_objects)

  run_config <- expand.grid(method = methods, mode = modes, stringsAsFactors = FALSE)
  run_config$slide_id <- slide_id
  run_config$input_file <- test_file
  run_config$k_components <- k_components
  run_config$n_iterations <- as.integer(cfg_run$n_iter)
  run_config$gene_perturb_frac <- cfg_run$gene_frac
  run_config$spot_perturb_frac <- cfg_run$spot_frac
  run_config$ncores_effective <- ifelse(
    run_config$method == "consICA" & run_config$mode != "pure_stability",
    consica_ncores, 1L)
  run_config <- run_config[, c("slide_id", "input_file", "method", "mode",
                               "k_components", "n_iterations",
                               "gene_perturb_frac", "spot_perturb_frac",
                               "ncores_effective")]

  # ---- data loading and preprocessing ----
  cat("\n[1/3] Loading slide:", slide_id, "\n")
  sce <- readH5AD(test_file)
  srt <- as.Seurat(sce, counts = "X", data = NULL)
  check_in_tissue(sce, slide_id)

  master_counts <- GetAssayData(srt, assay = "originalexp", layer = "counts")
  if (is.na(min_genes)) min_genes <- auto_spot_threshold(master_counts)
  keep_spots <- Matrix::colSums(master_counts) >= min_counts &
                Matrix::colSums(master_counts > 0) >= min_genes
  print(plot_spatial_spots(sce, slide_id))
  print(plot_spot_qc(sce, keep_spots, slide_id))
  rm(sce); gc()

  srt <- srt[, keep_spots]
  srt <- srt[rowSums(GetAssayData(srt, assay = "originalexp", layer = "counts")) > 0, ]
  cat(sprintf("    QC: kept %d / %d spots (min_genes=%d, min_counts=%d)\n",
              sum(keep_spots), length(keep_spots), min_genes, min_counts))

  # Drop technical control features (Xenium negative controls / blank codewords).
  # Harmless on Visium: the pattern matches nothing in real gene symbols.
  ctrl_pattern <- "^(BLANK|NegControlProbe|NegControlCodeword|antisense|DeprecatedCodeword|UnassignedCodeword|Intergenic)"
  n_ctrl <- sum(grepl(ctrl_pattern, rownames(srt)))
  srt <- srt[!grepl(ctrl_pattern, rownames(srt)), ]
  cat(sprintf("    Removed %d control features; %d real genes remain\n",
              n_ctrl, nrow(srt)))

  if (inherits(srt[["originalexp"]], "Assay5")) {
    srt[["originalexp"]]$counts <- as(srt[["originalexp"]]$counts, "dgCMatrix")
  } else {
    srt[["originalexp"]]@counts <- as(srt[["originalexp"]]@counts, "dgCMatrix")
  }

  cat("[2/3] SCTransform normalization...\n")
  srt_master <- suppressWarnings(
    SCTransform(srt, assay = "originalexp",
                return.only.var.genes = FALSE, verbose = FALSE))

  cat(sprintf("    Slide: %s | Spots: %d | Genes: %d\n",
              slide_id, ncol(srt_master), nrow(srt_master)))

  # Large slides: keep only the GO representative (pure_stability iter 1) to avoid OOM.
  # All other slides keep every consica_object (getGO flexibility). Add big slides here.
  # keep_all_consica is set from config/runs.csv (cfg_run$strip) in section 8.
  cat(sprintf("    consica_object: %s\n",
              if (keep_all_consica) "keep all iterations" else "keep GO representative only"))

  # ---- the twelve method and mode cells ----
  cat("[3/3] Running benchmark...\n\n")

  all_metrics <- list()
  all_metadata <- list()

  for (i in seq_len(nrow(run_config))) {
  #for (i in 5:nrow(run_config)) {

    cfg <- run_config[i, ]

    gp_pct <- if (cfg$mode == "gene_perturbation") cfg$gene_perturb_frac * 100 else NA
    tag <- make_save_tag(cfg$slide_id, cfg$method, cfg$mode,
                            cfg$k_components, cfg$n_iterations, gp_pct)

    cat(sprintf(" -> [%d/%d] Method: %-20s | Mode: %s\n",
                i, nrow(run_config), cfg$method, cfg$mode))

    result <- execute_benchmark_mode(
      srt_raw = srt_master,
      mode = cfg$mode,
      k = cfg$k_components,
      n_iter = cfg$n_iterations,
      method = cfg$method,
      gene_perturb_frac = cfg$gene_perturb_frac,
      spot_perturb_frac = cfg$spot_perturb_frac,
      ncores = cfg$ncores_effective,
      top_n = TOP_N_GENES,
      repro_thresholds = REPRO_THRESHOLDS,
      keep_all_consica = keep_all_consica)

    run_meta <- make_run_metadata(
      slide_id = cfg$slide_id,
      test_file = cfg$input_file,
      method = cfg$method,
      mode = cfg$mode,
      k = cfg$k_components,
      n_iter = cfg$n_iterations,
      gene_perturb_frac = cfg$gene_perturb_frac,
      spot_perturb_frac = cfg$spot_perturb_frac,
      ncores = cfg$ncores_effective,
      top_n = TOP_N_GENES,
      repro_thresholds  = REPRO_THRESHOLDS,
      min_genes = min_genes,
      min_counts = min_counts,
      srt_master = srt_master)

    # full object: metadata + metrics + runs (S_matrix, M_matrix, subset, timing,
    # and consica_object for consICA) + cross-run diagnostics (distributions)
    save_obj <- list(
      metadata = run_meta,
      metrics = result$metrics,
      runs = result$runs,
      diagnostics = result$diagnostics)

    runs_path <- file.path(output_dir, paste0("runs__", tag, ".rds"))
    saveRDS(save_obj, runs_path)
    cat(sprintf("    Saved: %s\n", runs_path))

    # accumulate metrics + metadata in memory (no re-reading of heavy .rds later)
    all_metrics[[tag]] <- result$metrics
    all_metadata[[tag]] <- run_meta

    # checkpoint CSV after every run
    metrics_so_far <- do.call(rbind, all_metrics)
    rownames(metrics_so_far) <- NULL
    write.csv(metrics_so_far,
              file.path(output_dir,
                        sprintf("metrics__checkpoint__%s__k%d.csv",
                                slide_id, k_components)),
              row.names = FALSE)
  }

  # ---- summary for this row ----
  final_dataframe <- do.call(rbind, all_metrics)
  rownames(final_dataframe) <- NULL

  final_metrics_path <- file.path(output_dir,
                                  sprintf("metrics_FINAL__%s__k%d__n%d__gp%dpct.csv",
                                          slide_id, k_components,
                                          run_config$n_iterations[1],
                                          round(run_config$gene_perturb_frac[1] * 100)))
  write.csv(final_dataframe, final_metrics_path, row.names = FALSE)

  # session log = metadata of all runs (accumulated in-memory)
  session_log_path <- file.path(output_dir,
                                sprintf("session_log__%s__k%d__%s.rds",
                                        slide_id, k_components,
                                        format(Sys.time(), "%Y%m%d_%H%M%S")))
  saveRDS(all_metadata, session_log_path)

  cat("\n=== BENCHMARK COMPLETE ===\n")
  print(final_dataframe)
  cat("\nFinal metrics :", final_metrics_path, "\n")
  cat("Session log   :", session_log_path,    "\n")

  rm(srt, srt_master); gc()
}

# ---- kept for reference ----

# The per-slide configuration as it was set by hand before config.R, and a snippet that
# rebuilds the summary table from the saved rds.

if (FALSE) {

TOP_N_GENES <- 50L          # top genes per component for Jaccard metric
REPRO_THRESHOLDS <- c(0.5, 0.7)  # matched-corr thresholds for reproducible-component count

methods <- c("consICA", "CFS_As_Is",
             "icafast_no_Rmat", "icafast_with_Rmat")#, "fastICA_Package")
modes <- c("pure_stability", "spot_perturbation", "gene_perturbation")

run_config <- expand.grid(method = methods, mode = modes,
                          stringsAsFactors = FALSE)

run_config$slide_id <- "TENX115"
run_config$input_file <- "data/HEST/skin/st/TENX115.h5ad"

run_config$slide_id <- "TENX117"
run_config$input_file <- "data/HEST/skin/st/TENX117.h5ad"

run_config$slide_id <- "MEND40"
run_config$input_file <- "data/HEST/skin/st/MEND40.h5ad"

run_config$slide_id <- "MEND38"
run_config$input_file <- "data/HEST/skin/st/MEND38.h5ad"

run_config$slide_id <- "NCBI600"
run_config$input_file <- "data/HEST/bladder/st/NCBI600.h5ad"

run_config$slide_id <- "NCBI776"
run_config$input_file <- "data/HEST/breast/st/NCBI776.h5ad"

run_config$k_components <- 50L
run_config$n_iterations <- 10L
run_config$gene_perturb_frac <- 0.05    # 0.05 = 5% (primary), 0.15 = 15% (supplementary)
run_config$spot_perturb_frac <- 0.05    # 0.05 = drop 5% of spots (keep 95%)

# One config file targets one slide. Derive scalars now (ncores lookup below needs slide_id).
stopifnot(length(unique(run_config$slide_id)) == 1L)
slide_id <- run_config$slide_id[1]
test_file <- run_config$input_file[1]
k_components <- run_config$k_components[1]

# ncores_effective = cores actually passed to the wrapper (source of truth).
# Only consICA uses it, and only in perturbation modes; pure_stability -> 1
# (its decompositions feed GO). Non-consICA methods run single-threaded.
# Per-slide consICA cores for perturbation modes (pure_stability is always 1).
# Large Visium needs 1: a 4-core PSOCK cluster runs out of memory.
CONSICA_NCORES <- c(TENX115 = 4L, TENX117 = 4L, MEND40 = 4L, NCBI600 = 4L, NCBI776 = 1L)
consica_ncores <- unname(CONSICA_NCORES[slide_id])
if (is.na(consica_ncores)) consica_ncores <- 1L
run_config$ncores_effective <- ifelse(
  run_config$method == "consICA" & run_config$mode != "pure_stability",
  consica_ncores, 1L)

run_config <- run_config[, c("slide_id", "input_file", "method", "mode",
                             "k_components", "n_iterations",
                             "gene_perturb_frac", "spot_perturb_frac",
                             "ncores_effective")]

output_dir <- "output/ica_perturbations/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# PER-SLIDE spot QC (verified once via plot_spot_qc; logged in metadata):
#   1L = keep-all (full-tissue slide)   NA = auto (slide has real background)
SPOT_MIN_GENES <- c(TENX115 = NA, TENX117 = NA,
                    MEND40 = 1L, MEND38 = 1L, NCBI600 = 1L, NCBI776 = 1L)
min_genes <- unname(SPOT_MIN_GENES[slide_id])   # NA (incl. unlisted slides) -> auto during loading
min_counts <- 0L

}  # end OLD config (disabled — see config/runs.csv above)

if(F){
  files <- list.files(output_dir, "^runs__NCBI776__.*\\.rds$", full.names = TRUE)
  length(files)
  final_dataframe <- do.call(rbind, lapply(files, function(f) readRDS(f)$metrics))
  rownames(final_dataframe) <- NULL
  print(final_dataframe)
}

# ---- what a saved run holds ----
#   obj <- readRDS("output/ica_perturbations/runs__TENX115__consICA__pure_stability__k20__n10.rds")
#   obj$metadata                         # run parameters, package versions, seeds
#   obj$metrics                          # one-row summary of all metrics
#   obj$runs[[1]]$S_matrix               # genes x components
#   obj$runs[[1]]$M_matrix               # components x spots
#   obj$runs[[1]]$consica_object         # full consICA object for getGO() (consICA only)
#   obj$diagnostics$cross_stability_scores   # per-pair distributions for figures
#   obj$diagnostics$matched_cor_by_pair      # matched correlations per run pair
