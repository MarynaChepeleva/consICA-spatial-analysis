# Setup shared by the control-arm scripts: the benchmark engine, arm A's provenance, the
# consICA version gate and load_master(). Sourced by bench_B_run.R, bench_B_run_gp15.R and
# component_gene_overlap.R.
#
# Sourcing this reads config/runs.csv and preflight_provenance_A.csv and runs the version
# gate, which is what the callers want. It starts no decomposition. What to run is left to
# the caller, so no slide or mode list is defined here.

suppressPackageStartupMessages({
  library(zellkonverter)
  library(Seurat)
  library(Matrix)
  library(consICA)
})

# A was computed on consICA 2.7.4; the library now holds 2.11.4, a later development build.
# Overriding the gate is a considered decision, checked against the package sources:
#   * coreICA.R is byte-identical between the two versions, so the ICA maths did not change;
#   * with ncores = 1 and `seed` not passed, 2.11.4 takes the single-core default branch,
#     the same loop as 2.7.4's one-core loop (same seq.int(1, ntry), same coreICA call);
#   * at ntry = 1 both versions return at the early `if (ntry == 1)` and never reach the
#     consensus assembly, the only part that was rewritten (which.max -> alignComponents,
#     in 2.11.1).
# So B's execution path is identical on both versions. A is NOT recomputed.
ALLOW_VERSION_MISMATCH <- TRUE

BENCH <- "scripts/benchmark/bench_ica_functions.R"
RUN_DIR <- "output/ica_perturbations"
OUT_DIR <- "output/benchmark_update_2026-09-23"
METHOD <- "consICA_ntry1"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- the benchmark's wrappers, metrics and engine ----
source(BENCH)
stopifnot(is.function(wrapper_consica_single), is.function(execute_benchmark_mode))

cfg_all <- read.csv("config/runs.csv", stringsAsFactors = FALSE)
prov <- read.csv(file.path(OUT_DIR, "preflight_provenance_A.csv"), stringsAsFactors = FALSE)

# ---- GATE 1: B must run on the same consICA as A -----------------------------------------
# A used 2.7.4. The greedy, non-bijective component matching (apply(r^2, 2, which.max)) was
# replaced by alignComponents() only in 2.11.1, so a newer package would change the consensus
# implementation itself and A - B would no longer isolate ntry.
a_ver <- unique(prov$consICA)
now_ver <- as.character(packageVersion("consICA"))
cat(sprintf("consICA: A ran on %s | installed now %s\n", paste(a_ver, collapse = "/"), now_ver))
if (length(a_ver) != 1L || now_ver != a_ver) {
  msg <- sprintf(paste0("VERSION MISMATCH: A ran on %s, installed is %s. B would use a different\n",
                        "consensus implementation, so A - B would not isolate ntry. Install %s, or\n",
                        "set ALLOW_VERSION_MISMATCH <- TRUE and record this in the limitations."),
                 paste(a_ver, collapse = "/"), now_ver, paste(a_ver, collapse = "/"))
  if (!ALLOW_VERSION_MISMATCH) stop(msg) else warning(msg)
}

# ---- load one slide exactly as the driver does, in its data loading section -------------
# Verbatim copy minus the two diagnostic plots; guarded by a dimension check against A.
load_master <- function(cfg) {
  sce <- readH5AD(cfg$file)
  srt <- as.Seurat(sce, counts = "X", data = NULL)
  mc <- GetAssayData(srt, assay = "originalexp", layer = "counts")
  min_genes <- if (is.na(cfg$min_genes)) auto_spot_threshold(mc) else as.integer(cfg$min_genes)
  min_counts <- 0L
  keep_spots <- Matrix::colSums(mc) >= min_counts & Matrix::colSums(mc > 0) >= min_genes
  rm(sce); gc()
  srt <- srt[, keep_spots]
  srt <- srt[rowSums(GetAssayData(srt, assay = "originalexp", layer = "counts")) > 0, ]
  ctrl <- "^(BLANK|NegControlProbe|NegControlCodeword|antisense|DeprecatedCodeword|UnassignedCodeword|Intergenic)"
  srt <- srt[!grepl(ctrl, rownames(srt)), ]
  if (inherits(srt[["originalexp"]], "Assay5")) {
    srt[["originalexp"]]$counts <- as(srt[["originalexp"]]$counts, "dgCMatrix")
  } else {
    srt[["originalexp"]]@counts <- as(srt[["originalexp"]]@counts, "dgCMatrix")
  }
  m <- suppressWarnings(SCTransform(srt, assay = "originalexp",
                                    return.only.var.genes = FALSE, verbose = FALSE))
  # GATE 2: reloaded master must match the dimensions A recorded for this slide
  a <- prov[prov$slide == cfg$slide, ][1, ]
  if (!(nrow(m) == a$n_genes_master && ncol(m) == a$n_spots_master))
    stop(sprintf("%s: reloaded master is %dx%d but A recorded %dx%d — loading drifted.",
                 cfg$slide, nrow(m), ncol(m), a$n_genes_master, a$n_spots_master))
  cat(sprintf("    master %d genes x %d spots (min_genes=%d) — matches A\n",
              nrow(m), ncol(m), min_genes))
  list(srt = m, min_genes = min_genes, min_counts = min_counts)
}
