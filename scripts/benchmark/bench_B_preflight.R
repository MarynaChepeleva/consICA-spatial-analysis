# Preflight for control B (consICA, ntry = 1): can B be compared with A at all?
# READ-ONLY: no ICA is run, no benchmark is re-run, no existing file is modified.
# Everything written goes to output/benchmark_update_2026-09-23/ only.
#
# A's provenance is taken from the tiny session_log__*.rds (metadata of every run of a
# slide, written by the driver) instead of the runs__*.rds, which total 43 GB and cannot be
# read partially (NCBI776 pure alone is 13 GB).
#
# The second half asks whether B can reproduce A's perturbation inputs. Two separate questions:
#   (a) do the SUBSETS match what A actually used?  -> compared against A's SAVED run$subset,
#       i.e. across sessions, which is the thing that decides A/B comparability. The subset is
#       the first RNG draw after set.seed(s), so it is reproduced here by replaying exactly that
#       draw from prep_iteration_input without paying for SCTransform.
#   (b) is the spot-mode matrix X itself deterministic? X is never saved, so this can only be
#       tested within one session: call prep_iteration_input twice under the same seed and
#       compare bitwise. gene_perturbation refits nothing, so only spot mode is at risk.

suppressPackageStartupMessages({
  library(zellkonverter)
  library(Seurat)
  library(Matrix)
})

BENCH <- "scripts/benchmark/bench_ica_functions.R"
RUN_DIR <- "output/ica_perturbations"
OUT_DIR <- "output/benchmark_update_2026-09-23"
SLIDE_X <- "TENX117"      # smallest slide, used for the determinism probe
SEEDS_X <- 1L             # seed(s) for the expensive bitwise-X check (each costs 2 SCTransform)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

nz <- function(x, default = NA) if (is.null(x) || length(x) != 1L) default else x

# ---- the benchmark's wrappers, metrics and engine ----
source(BENCH)
stopifnot(is.function(prep_iteration_input), is.function(auto_spot_threshold))
cat(sprintf("Loaded definitions from %s\n\n", BENCH))

cfg_all <- read.csv("config/runs.csv", stringsAsFactors = FALSE)

# Provenance of the existing A runs, from the session_log metadata ----
logs <- list.files(RUN_DIR, pattern = "^session_log__.*__k20__.*\\.rds$", full.names = TRUE)
prov <- do.call(rbind, lapply(logs, function(f) {
  ml <- readRDS(f)
  do.call(rbind, lapply(names(ml), function(tag) {
    m <- ml[[tag]]
    if (!identical(nz(m$method), "consICA")) return(NULL)    # A only
    pp <- m$method_params; pv <- m$package_versions
    data.frame(
      log_file   = basename(f), run_id = nz(m$run_id), computed = nz(m$timestamp),
      slide = nz(m$slide_id), mode = nz(m$mode),
      k = nz(m$k_components), n_iter = nz(m$n_iterations),
      ntry = nz(pp$ntry), ncores = nz(pp$ncores),
      gene_perturb_pct = nz(m$gene_perturb_pct), spot_perturb_frac = nz(m$spot_perturb_frac),
      assay = nz(m$assay_used), normalization = nz(m$normalization),
      min_genes = nz(m$spot_qc$min_genes), min_counts = nz(m$spot_qc$min_counts),
      n_genes_master = nz(m$n_genes), n_spots_master = nz(m$n_spots),
      seeds_used = paste(m$seeds_used, collapse = ","),
      r_version = nz(m$r_version),
      consICA = nz(pv$consICA), ica = nz(pv$ica), Seurat = nz(pv$Seurat),
      fastICA = nz(pv$fastICA), clue = nz(pv$clue),
      stringsAsFactors = FALSE)
  }))
}))

# expected per-iteration dims. DERIVED, not measured: the true n_genes_eff / n_spots_eff live
# in the 43 GB runs__*.rds. Spot mode refits SCTransform on the subset, which may itself drop
# genes, so its gene count is left NA rather than guessed.
prov$n_spots_eff_expected <- ifelse(prov$mode == "spot_perturbation",
                                    floor((1 - prov$spot_perturb_frac) * prov$n_spots_master),
                                    prov$n_spots_master)
prov$n_genes_eff_expected <- ifelse(
  prov$mode == "gene_perturbation", floor((1 - prov$gene_perturb_pct / 100) * prov$n_genes_master),
  ifelse(prov$mode == "spot_perturbation", NA_real_, prov$n_genes_master))

prov <- prov[order(prov$slide, prov$mode, prov$gene_perturb_pct), ]
write.csv(prov, file.path(OUT_DIR, "preflight_provenance_A.csv"), row.names = FALSE)

cat("=== STEP 1.1: provenance of A (consICA runs, k=20) ===\n")
print(prov[, c("slide", "mode", "gene_perturb_pct", "ntry", "ncores",
               "n_genes_master", "n_spots_master", "seeds_used", "consICA", "computed")],
      row.names = FALSE)
cat("\nDistinct consICA versions :", paste(unique(prov$consICA), collapse = ", "), "\n")
cat("Distinct R versions       :", paste(unique(prov$r_version), collapse = " | "), "\n")
cat("Cells run with ncores > 1 :", sum(prov$ncores > 1, na.rm = TRUE), "of", nrow(prov), "\n\n")

# Can B reproduce A's perturbation inputs? ----
cfg <- cfg_all[cfg_all$run_id == paste0(SLIDE_X, "__k20__gp05"), ]
stopifnot(nrow(cfg) == 1L)

# ---- load the slide exactly as the driver does, in its data loading section --------------
# Verbatim copy; the only omission is the two diagnostic plots, which do not touch srt_master.
# Guarded below by an assertion against A's recorded master dimensions.
cat("=== STEP 1.0: loading", SLIDE_X, "===\n")
sce <- readH5AD(cfg$file)
srt <- as.Seurat(sce, counts = "X", data = NULL)
master_counts <- GetAssayData(srt, assay = "originalexp", layer = "counts")
min_genes <- if (is.na(cfg$min_genes)) auto_spot_threshold(master_counts) else as.integer(cfg$min_genes)
min_counts <- 0L
keep_spots <- Matrix::colSums(master_counts) >= min_counts &
              Matrix::colSums(master_counts > 0) >= min_genes
rm(sce); gc()
srt <- srt[, keep_spots]
srt <- srt[rowSums(GetAssayData(srt, assay = "originalexp", layer = "counts")) > 0, ]
ctrl_pattern <- "^(BLANK|NegControlProbe|NegControlCodeword|antisense|DeprecatedCodeword|UnassignedCodeword|Intergenic)"
srt <- srt[!grepl(ctrl_pattern, rownames(srt)), ]
if (inherits(srt[["originalexp"]], "Assay5")) {
  srt[["originalexp"]]$counts <- as(srt[["originalexp"]]$counts, "dgCMatrix")
} else {
  srt[["originalexp"]]@counts <- as(srt[["originalexp"]]@counts, "dgCMatrix")
}
srt_master <- suppressWarnings(SCTransform(srt, assay = "originalexp",
                                           return.only.var.genes = FALSE, verbose = FALSE))
cat(sprintf("    srt_master: %d genes x %d spots (min_genes=%d)\n",
            nrow(srt_master), ncol(srt_master), min_genes))

# ---- GATE: the reloaded master must match the dimensions A recorded ----------------------
a_row <- prov[prov$slide == SLIDE_X & prov$mode == "pure_stability", ][1, ]
dim_ok <- nrow(srt_master) == a_row$n_genes_master && ncol(srt_master) == a_row$n_spots_master
cat(sprintf("    A recorded: %d genes x %d spots -> %s\n", a_row$n_genes_master,
            a_row$n_spots_master, if (dim_ok) "MATCH" else "MISMATCH"))
if (!dim_ok)
  stop("Loading drifted: reloaded master differs from the one A used. Resolve before running B.")

# ---- (a) subsets vs A's SAVED subsets, all 10 seeds, no SCTransform ----------------------
# replays only the first RNG draw of prep_iteration_input, the sample() call of each mode
seeds <- seq_len(as.integer(cfg$n_iter))
all_cells <- colnames(srt_master)
all_genes <- rownames(srt_master[["SCT"]])

sub_check <- do.call(rbind, lapply(c("spot_perturbation", "gene_perturbation"), function(md) {
  f <- file.path(RUN_DIR, sprintf("runs__%s__consICA__%s__k20__n10%s.rds", SLIDE_X, md,
                                  if (md == "gene_perturbation") "_gp5pct" else ""))
  if (!file.exists(f)) { message("  missing A run: ", f); return(NULL) }
  a_runs <- readRDS(f)$runs
  out <- do.call(rbind, lapply(seeds, function(s) {
    set.seed(s)
    keep <- if (md == "spot_perturbation")
      sample(all_cells, size = floor((1 - cfg$spot_frac) * length(all_cells)))
    else
      sample(all_genes, size = floor((1 - cfg$gene_frac) * length(all_genes)))
    a_keep <- a_runs[[s]]$subset
    data.frame(slide = SLIDE_X, mode = md, seed = s,
               n_regenerated = length(keep), n_saved_A = length(a_keep),
               identical_order = identical(as.character(keep), as.character(a_keep)),
               identical_as_set = setequal(as.character(keep), as.character(a_keep)),
               stringsAsFactors = FALSE)
  }))
  rm(a_runs); gc(); out
}))
write.csv(sub_check, file.path(OUT_DIR, "preflight_subset_match.csv"), row.names = FALSE)
cat("\n--- (a) regenerated subsets vs A's saved subsets ---\n")
print(sub_check, row.names = FALSE)

# ---- (b) is spot-mode X bitwise deterministic? (X is never saved -> in-session only) ------
cat("\n--- (b) spot_perturbation X determinism (2 full prep calls per seed) ---\n")
x_check <- do.call(rbind, lapply(SEEDS_X, function(s) {
  set.seed(s); p1 <- prep_iteration_input(srt_master, "spot_perturbation", cfg$spot_frac, cfg$gene_frac)
  set.seed(s); p2 <- prep_iteration_input(srt_master, "spot_perturbation", cfg$spot_frac, cfg$gene_frac)
  same_sub <- identical(p1$subset, p2$subset)
  same_dim <- identical(dim(p1$X), dim(p2$X))
  same_dnm <- identical(dimnames(p1$X), dimnames(p2$X))
  bitwise <- same_dim && same_dnm && identical(p1$X, p2$X)
  maxdiff <- if (same_dim) max(abs(p1$X - p2$X)) else NA_real_
  fracdiff <- if (same_dim) mean(p1$X != p2$X) else NA_real_
  d <- data.frame(slide = SLIDE_X, seed = s, subset_identical = same_sub,
                  dim_identical = same_dim, dimnames_identical = same_dnm,
                  X_bitwise_identical = bitwise, max_abs_diff = maxdiff,
                  frac_elements_differing = fracdiff,
                  n_genes_eff = nrow(p1$X), n_spots_eff = ncol(p1$X),
                  stringsAsFactors = FALSE)
  rm(p1, p2); gc(); d
}))
write.csv(x_check, file.path(OUT_DIR, "preflight_X_determinism.csv"), row.names = FALSE)
print(x_check, row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "preflight_sessionInfo.txt"))

# ---- verdict ------------------------------------------------------------------------------
sub_ok <- all(sub_check$identical_order)
x_ok <- all(x_check$X_bitwise_identical)
cat("\n=== VERDICT ===\n")
cat(" subsets reproduce A exactly :", sub_ok, "\n")
cat(" spot-mode X bitwise stable  :", x_ok, "\n")
if (sub_ok && x_ok) {
  cat(" -> unit and gene conditions are sound; proceed to step 1.2 (run B).\n")
} else {
  cat(" -> STOP. Do not run B yet.\n")
  if (!sub_ok) cat("    subsets differ: B would perturb different spots/genes than A.\n")
  if (!x_ok)   cat("    SCTransform is not bitwise reproducible: either drop the unit condition\n",
                   "   or save the normalised matrices explicitly for both A and B.\n", sep = "")
}
cat("\nWrote to", normalizePath(OUT_DIR, mustWork = FALSE), "\n")
if (!(sub_ok && x_ok)) stop("preflight failed; see the verdict above")
