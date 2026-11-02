# Recompute the §3.3 interpretability metric under Benjamini-Hochberg FDR, replacing
# the Bonferroni bound of go_interpretability_fdr_check.R. The saved go_benchmark rds
# keep only the top-20 GenTable rows per component (run_topGO_one_component in
# GO_functions.R), so the full p vector an exact BH needs is not recoverable from them.
# topGO is therefore re-run here from the saved pure_stability S_matrices; no consICA is
# re-run, since S_matrix is stored for iteration 1 of every method and slide.
# score(result) returns the untruncated per-term p vector.
#
# Two correction branches per component family, the terms tested for that component:
#   BH_classic  (HEADLINE)  = p.adjust(score(classic Fisher),  "BH")
#   BH_weight01 (sensitivity)= p.adjust(score(weight01 Fisher), "BH")
# classic is the headline because all four configurations get it identically. Neither
# branch reproduces the demonstration annotations, where consICA::getGO runs topGO's
# weight algorithm. topGO states weight01 p are not meant for multiple-testing, and
# classic p over the GO graph are dependent too: the correction is an order of
# magnitude, not an exact probability.
#
# READ-ONLY except the three CSVs + one versions file it writes into output/benchmark_summary/.
# Does NOT touch GO_pipeline.R, the saved go_benchmark__*.rds, or the demo annotations.

suppressPackageStartupMessages({
  library(topGO)
  library(org.Hs.eg.db)      # all five benchmark slides are human; xenograft is out of scope
})

run_dir <- "output/ica_perturbations"
out_dir <- "output/benchmark_summary"
old_csv <- file.path(out_dir, "go_interpretability_fdr_check.csv")   # Bonferroni check, for the sanity gate
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

K <- 20L
N_ITER <- 10L
REP_ITER <- 1L                        # the representative pure_stability iteration (as in GO_pipeline.R)
P_THR <- 0.05
N_BONF <- c(5000, 15000)            # reproduce the old Bonferroni columns for continuity/checking
NODE <- 10L                       # nodeSize, as GO_pipeline.R
ONT <- "BP"
ORGDB <- "org.Hs.eg.db"
IDTYPE <- "symbol"

# 4 methods present on disk (fastICA_Package excluded from the benchmark; no pure_stability run)
METHODS <- c("consICA", "CFS_As_Is", "icafast_no_Rmat", "icafast_with_Rmat")
# 5 slides; per-slide top_n genes == GO_pipeline.R (Xenium panel ~300 genes -> 40; Visium -> 100)
SLIDES <- list(TENX115 = 40L, TENX117 = 40L, MEND40 = 100L, NCBI776 = 100L, NCBI600 = 100L)

# ---- gene selection: byte-identical to select_genes_per_component in GO_functions.R ----
# duplicated on purpose: sourcing GO_pipeline.R whole would
# run its full driver and overwrite outputs. Keep these two in sync with GO_pipeline.R.
.orient_positive <- function(w) {                 # flip so the heavy (signal) tail is positive
  s <- mean(((w - mean(w)) / sd(w))^3)            # skewness
  if (is.finite(s) && s < 0) -w else w
}
select_genes_per_component <- function(S_matrix, top_n) {
  apply(S_matrix, 2, function(w) {
    wo <- .orient_positive(w)
    names(sort(wo, decreasing = TRUE))[seq_len(min(top_n, length(wo)))]
  }, simplify = FALSE)
}

# ---- topGO for one component: full p vectors for BOTH algorithms, plus m ----
topGO_one_component <- function(top_genes, all_genes) {
  if (length(top_genes) == 0) return(NULL)
  gene_list <- factor(as.integer(all_genes %in% top_genes)); names(gene_list) <- all_genes
  go_data <- suppressMessages(new(
    "topGOdata", ontology = ONT, allGenes = gene_list, geneSel = function(x) x == 1,
    annot = annFUN.org, mapping = ORGDB, ID = IDTYPE, nodeSize = NODE))
  w01 <- suppressMessages(runTest(go_data, algorithm = "weight01", statistic = "fisher"))
  cls <- suppressMessages(runTest(go_data, algorithm = "classic",  statistic = "fisher"))
  list(m         = length(usedGO(go_data)),   # BP terms tested = this component's family size
       p_weight  = score(w01),                # full named vector, untruncated
       p_classic = score(cls))
}

# ---- per-component recompute across every slide x method ----
rows_comp <- list()
for (slide in names(SLIDES)) {
  top_n <- SLIDES[[slide]]
  for (mth in METHODS) {
    path <- file.path(run_dir, sprintf("runs__%s__%s__pure_stability__k%d__n%d.rds",
                                       slide, mth, K, N_ITER))
    if (!file.exists(path)) { warning("missing: ", path); next }
    S <- readRDS(path)$runs[[REP_ITER]]$S_matrix
    genes_all <- rownames(S)
    gene_sets <- select_genes_per_component(S, top_n)
    cat(sprintf("%-8s %-18s  %d comp, top_n=%d, background=%d genes\n",
                slide, mth, ncol(S), top_n, length(genes_all)))
    for (comp in colnames(S)) {
      go <- tryCatch(topGO_one_component(gene_sets[[comp]], genes_all),
                     error = function(e) { message("  ", comp, ": ", e$message); NULL })
      if (is.null(go)) {
        rows_comp[[length(rows_comp) + 1L]] <- data.frame(
          Slide = slide, Method = mth, Component = comp, n_genes = length(gene_sets[[comp]]),
          m_tested = NA_integer_, best_p_weight01 = NA_real_, best_p_classic = NA_real_,
          minBH_weight01 = NA_real_, minBH_classic = NA_real_, stringsAsFactors = FALSE)
        next
      }
      rows_comp[[length(rows_comp) + 1L]] <- data.frame(
        Slide = slide, Method = mth, Component = comp, n_genes = length(gene_sets[[comp]]),
        m_tested        = go$m,
        best_p_weight01 = min(go$p_weight,  na.rm = TRUE),               # reproduces old Frac_nominal / bonf
        best_p_classic  = min(go$p_classic, na.rm = TRUE),
        minBH_weight01  = min(p.adjust(go$p_weight,  "BH"), na.rm = TRUE),
        minBH_classic   = min(p.adjust(go$p_classic, "BH"), na.rm = TRUE),
        stringsAsFactors = FALSE)
    }
  }
}
comp <- do.call(rbind, rows_comp)
# interpretable flags (>=1 term below threshold in that branch)
comp$interp_nominal <- comp$best_p_weight01 < P_THR                     # weight01 nominal
comp$interp_BH_classic <- comp$minBH_classic   < P_THR
comp$interp_BH_weight <- comp$minBH_weight01  < P_THR

write.csv(comp, file.path(out_dir, "go_interpretability_bh__per_component__k20.csv"), row.names = FALSE)

# ---- summary: one row per slide x method ----
agg <- function(v, f) tapply(v, list(comp$Slide, comp$Method), f)
frac_lt <- function(x, thr) mean(x < thr, na.rm = TRUE)

keys <- unique(comp[c("Slide", "Method")])
summ <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
  s <- keys$Slide[i]; m <- keys$Method[i]
  d <- comp[comp$Slide == s & comp$Method == m, ]
  data.frame(
    Slide = s, Method = m, n_comp = nrow(d),
    Frac_nominal          = frac_lt(d$best_p_weight01, P_THR),            # weight01, reproduces old CSV
    Frac_bonf5000         = frac_lt(d$best_p_weight01, P_THR / N_BONF[1]),
    Frac_bonf15000        = frac_lt(d$best_p_weight01, P_THR / N_BONF[2]),
    Interp_Frac_BH_classic  = frac_lt(d$minBH_classic,  P_THR),          # HEADLINE
    Interp_Frac_BH_weight01 = frac_lt(d$minBH_weight01, P_THR),          # sensitivity
    worst_best_p          = max(d$best_p_weight01, na.rm = TRUE),         # analog of old Interp_worst_best_p
    worst_minBH_classic   = max(d$minBH_classic,  na.rm = TRUE),         # least-interpretable comp, headline branch
    worst_minBH_weight01  = max(d$minBH_weight01, na.rm = TRUE),
    m_min = min(d$m_tested, na.rm = TRUE), m_med = median(d$m_tested, na.rm = TRUE),
    m_max = max(d$m_tested, na.rm = TRUE),
    stringsAsFactors = FALSE)
}))
summ <- summ[order(match(summ$Slide, names(SLIDES)), summ$Method), ]
write.csv(summ, file.path(out_dir, "go_interpretability_bh__summary__k20.csv"), row.names = FALSE)

# ---- package versions, reported in Methods ----
ver <- c(
  sprintf("date            : %s", format(Sys.Date())),
  sprintf("R               : %s", R.version.string),
  sprintf("topGO           : %s", as.character(packageVersion("topGO"))),
  sprintf("org.Hs.eg.db    : %s", as.character(packageVersion("org.Hs.eg.db"))),
  sprintf("GO.db           : %s", as.character(packageVersion("GO.db"))))
writeLines(ver, file.path(out_dir, "go_interpretability_bh__versions.txt"))
cat("\n", paste(ver, collapse = "\n"), "\n\n", sep = "")

# ---- SANITY GATE: weight01 recompute must reproduce the old Bonferroni CSV exactly ----
# (TENX115, TENX117, MEND40 x 4 methods). Mismatch => version drift; STOP before trusting BH.
if (file.exists(old_csv)) {
  old <- read.csv(old_csv, stringsAsFactors = FALSE)
  chk <- merge(old, summ, by = c("Slide", "Method"), suffixes = c(".old", ".new"))
  cols <- c("Frac_nominal", "Frac_bonf5000", "Frac_bonf15000")
  diffs <- do.call(rbind, lapply(cols, function(cc) {
    d <- abs(chk[[paste0(cc, ".old")]] - chk[[paste0(cc, ".new")]])
    data.frame(col = cc, Slide = chk$Slide, Method = chk$Method,
               old = chk[[paste0(cc, ".old")]], new = chk[[paste0(cc, ".new")]], d = d)
  }))
  bad <- diffs[is.na(diffs$d) | diffs$d > 1e-9, ]
  if (nrow(bad)) {
    print(bad)
    stop("SANITY GATE FAILED: weight01 recompute does not match the old Bonferroni CSV. ",
         "Do NOT trust the BH numbers. Check org.Hs.eg.db / topGO / GO.db versions (written above) ",
         "against those used for the saved go_benchmark rds.")
  }
  cat(sprintf("SANITY GATE PASSED: weight01 fracs reproduce the old CSV on %d slide x method rows.\n",
              nrow(chk)))
} else {
  cat("SANITY GATE SKIPPED: old", old_csv, "not found.\n")
}

cat("\nWrote:\n",
    " ", file.path(out_dir, "go_interpretability_bh__summary__k20.csv"), "\n",
    " ", file.path(out_dir, "go_interpretability_bh__per_component__k20.csv"), "\n",
    " ", file.path(out_dir, "go_interpretability_bh__versions.txt"), "\n", sep = "")
print(summ, row.names = FALSE, digits = 4)
