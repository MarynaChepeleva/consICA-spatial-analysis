# Supplementary Table S3, the full method comparison, and the source of Table 2.
#
# One row per slide, configuration and condition. Nothing is recomputed: every number is
# read from an output that already exists.
#
# Three kinds of absence, as the S3 caption promises to distinguish:
#   1. not applicable to the condition -> empty cell. Moran, BH annotation, semantic uniqueness
#      and gene overlap are computed on one unperturbed decomposition, so they are empty under
#      perturbation and at k = 50. This is the bulk of the blanks and is explained in the caption.
#   2. configuration not run -> no row at all. The single-run control does not exist at k = 50,
#      nor outside the gene condition at 15%, so the caption states where it is missing: with
#      no column recording it, an absent row is visible only as a short block of four.
#   3. computed for other configurations, absent here -> the literal token "not computed".
#      No cell carries it at present; the columns that can are written as text so that one
#      could be marked without turning its neighbours into text as well.
# Case 1 is written as R's NA in the CSV, which read.csv returns as a missing value and
# openxlsx then writes as a genuinely blank cell, so the workbook shows the empty cell.

suppressPackageStartupMessages({ library(dplyr) })
select <- dplyr::select

RUN_DIR <- "output/ica_perturbations"
B_DIR <- "output/benchmark_update_2026-09-23"
OUT_DIR <- "output/tables"
MOR_A <- "output/benchmark_moran/moran_by_method__summary__k20.csv"
BH_A <- "output/benchmark_summary/go_interpretability_bh__summary__k20.csv"
SEM_DIR <- "output/go_enrichment"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SLIDES <- c("TENX115", "TENX117", "MEND40", "NCBI776", "NCBI600")
NOTC <- "not computed"

CONFIG <- c(consICA = "consICA",
            consICA_ntry1 = "consICA, single run",
            icafast_no_Rmat = "icafast, fixed init.",
            icafast_with_Rmat = "icafast, random init.",
            CFS_As_Is = "ICASpatial (CellsFromSpace)")
COND <- c(pure_stability = "Full data",
            spot_perturbation = "95% units",
            gene_perturbation = "95% genes")

# what the configuration estimates, written out where Table 2 abbreviates it to GS and SP.
# Read off the method wrappers in bench_ica_functions.R: consICA takes X itself and
# returns the gene signatures directly, while the three others take t(X), estimate on the
# unit axis, and recover genes by projecting X onto that estimate. The CFS configuration does the
# same inside custom_ICASpatial, which is why it reproduces icafast exactly.
AXIS <- c(consICA = "Gene signatures",
            consICA_ntry1 = "Gene signatures",
            icafast_no_Rmat = "Spatial profiles",
            icafast_with_Rmat = "Spatial profiles",
            CFS_As_Is = "Spatial profiles")

# legacy intra-uniqueness averages |cor| over all k^2 cells with the diagonal zeroed; the
# reported measure is the mean over distinct pairs. Exact, and the factor depends on k.
to_offdiag <- function(legacy, k) 1 - (k / (k - 1)) * (1 - legacy)

# ---- rows: the four original configurations from metrics_FINAL, the control from its own runs ----
mf <- list.files(RUN_DIR, pattern = "^metrics_FINAL__.*__n10__gp(5|15)pct\\.csv$", full.names = TRUE)
mf <- mf[grepl(paste(SLIDES, collapse = "|"), basename(mf))]        # drops MEND38

A <- do.call(rbind, lapply(mf, function(f) {
  b <- basename(f)
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  d$slide <- sub("^metrics_FINAL__([^_]+)__.*$", "\\1", b)
  d$k <- as.integer(sub("^.*__k([0-9]+)__.*$", "\\1", b))
  d$gp_pct <- as.integer(sub("^.*__gp([0-9]+)pct\\.csv$", "\\1", b))
  d
}))
A <- A[A$Method %in% names(CONFIG), ]
A <- A[!(A$gp_pct == 15 & A$Mode != "gene_perturbation"), ] # 15% is reported for genes only

B <- do.call(rbind, lapply(
  list.files(B_DIR, pattern = "^runs__.*consICA_ntry1__.*__n10.*\\.rds$", full.names = TRUE),
  function(f) {
    o <- readRDS(f)
    d <- o$metrics
    d$slide <- o$metadata$slide_id
    d$k <- as.integer(o$metadata$k_components)
    d$gp_pct <- if (grepl("gp15", basename(f))) 15L else 5L
    rm(o)
    gc()
    return(d)
  }))

keep <- c("slide", "Method", "k", "gp_pct", "Mode", "Matched_Stability", "SD_Matched",
          "TopGene_Jaccard", "SD_TopGene", "N_Repro_50", "N_Repro_70", "Intra_Uniqueness")
tab <- rbind(A[, keep], B[, keep])

tab <- tab %>%
  transmute(
    Slide = slide,
    Configuration = unname(CONFIG[Method]),
    Method_id = Method,
    Source_axis = unname(AXIS[Method]), # a property of the configuration: filled on every row
    k = k,
    # the gpXpct file tag travels with every mode, but genes are excluded only in the gene
    # condition; carrying 5 into the other two would contradict the Condition column
    gene_excluded_pct = ifelse(Mode == "gene_perturbation", gp_pct, 0L),
    Condition = unname(COND[Mode]),
    Matched_stability = Matched_Stability, SD_matched = SD_Matched,
    TopGene_Jaccard = TopGene_Jaccard, SD_TopGene = SD_TopGene,
    N_matched_above_0.5 = N_Repro_50, N_matched_above_0.7 = N_Repro_70,
    Intra_uniqueness = to_offdiag(Intra_Uniqueness, k)
  )

# the label shows the share actually kept, so the 15% rows read "85% genes"
gp <- tab$Condition == COND[["gene_perturbation"]]
tab$Condition[gp] <- sprintf("%d%% genes", 100L - tab$gene_excluded_pct[gp])

# ordering key: ranked by the mode rather than by the label, which is no longer a fixed string
cond_rank <- function(x) ifelse(x == COND[["pure_stability"]], 1L,
                        ifelse(x == COND[["spot_perturbation"]], 2L, 3L))

# ---- single-decomposition measures: unperturbed k = 20 only, empty elsewhere by design ----
mor <- rbind(
  read.csv(MOR_A, stringsAsFactors = FALSE) %>%
    transmute(Slide = slide, Method_id = method, Moran_I_median = median_Moran_I,
              Moran_I_min = min_Moran_I, Moran_I_max = max_Moran_I,
              frac_components_autocorrelated = frac_perm_sig),
  read.csv(file.path(B_DIR, "moran_B__summary__k20.csv"), stringsAsFactors = FALSE) %>%
    transmute(Slide = slide, Method_id = method, Moran_I_median = median_Moran_I,
              Moran_I_min = min_Moran_I, Moran_I_max = max_Moran_I,
              frac_components_autocorrelated = frac_perm_sig))

# the control is a single stochastic run, so its representative draw is placed against the
# spread of the same statistic over all ten repeats
rep_spread <- read.csv(file.path(B_DIR, "moran_B__across_repeats.csv"), stringsAsFactors = FALSE) %>%
  group_by(Slide = slide) %>%
  summarise(Moran_I_median_over_repeats = median(median_Moran_I),
            Moran_I_median_min_over_repeats = min(median_Moran_I),
            Moran_I_median_max_over_repeats = max(median_Moran_I), .groups = "drop") %>%
  mutate(Method_id = "consICA_ntry1")

# n_comp is the denominator of the annotated fraction and equals k on every row, so only the
# numerator is carried out; the reader divides by the k already in the table
bh <- rbind(
  read.csv(BH_A, stringsAsFactors = FALSE) %>%
    transmute(Slide, Method_id = Method,
              n_components_annotated = round(Interp_Frac_BH_classic * n_comp)),
  read.csv(file.path(B_DIR, "go_interpretability_bh_B__summary__k20.csv"), stringsAsFactors = FALSE) %>%
    filter(arm == "B") %>%
    transmute(Slide, Method_id = Method,
              n_components_annotated = round(Interp_Frac_BH_classic * n_comp)))

# SEM_DIR holds four slides; the fifth lives in the control's own file under arm = "A".
# Filtering that file to arm == "B" would drop the NCBI600 consensus value, so both are read
# whole and the duplicate keys resolved afterwards.
sem <- rbind(
  do.call(rbind, lapply(list.files(SEM_DIR, pattern = "^go_semantic_uniqueness__.*__k20\\.csv$",
                                   full.names = TRUE), read.csv, stringsAsFactors = FALSE)) %>%
    transmute(Slide, Method_id = Method, SemUniq_rank = SemUniq_topN, SemUniq_zscore = SemUniq_zscr),
  read.csv(file.path(B_DIR, "go_semantic_uniqueness_B__k20.csv"), stringsAsFactors = FALSE) %>%
    transmute(Slide, Method_id = Method, SemUniq_rank = SemUniq_topN, SemUniq_zscore = SemUniq_zscr))
sem <- sem[!duplicated(sem[, c("Slide", "Method_id")]), ]   # SEM_DIR wins where both have a key

# Overlap of the leading-gene lists between the components of one decomposition. Two Jaccards
# therefore share the table: TopGene_Jaccard is between repeats, over components paired by the
# Hungarian assignment, while Overlap_mean_jaccard is within one decomposition, between
# components that were never meant to pair. Overlap_top_n travels as a column because the two
# platforms are read off different list lengths, 40 against 100, and so are not comparable.
ov <- read.csv(file.path(B_DIR, "cfs_mechanism__1_topgene_overlap.csv"), stringsAsFactors = FALSE)
ov$Method_id <- ifelse(grepl("^consICA \\(gene axis", ov$config), "consICA",
                ifelse(grepl("^consICA_ntry1", ov$config), "consICA_ntry1", ov$config))
ov <- ov %>%
  filter(Method_id %in% names(CONFIG)) %>%
  transmute(Slide = slide, Configuration = unname(CONFIG[Method_id]), Method_id, k = k,
            top_n_genes = top_n, n_component_pairs = n_pairs,
            mean_pairwise_jaccard = as.character(round(mean_pairwise_jaccard, 4)),
            max_pairwise_jaccard  = as.character(round(max_pairwise_jaccard, 4)),
            union_of_lists = union_size, k_times_top_n = k_times_N,
            union_ratio = as.character(round(union_ratio, 4)))

# keyed on what the overlap file holds rather than on a fixed pair of slides, so a later run
# that covers more of them closes those cells instead of being shadowed by a synthetic row
missing_ov <- expand.grid(Slide = SLIDES, Method_id = names(CONFIG), stringsAsFactors = FALSE)
missing_ov <- missing_ov[!paste(missing_ov$Slide, missing_ov$Method_id) %in%
                           paste(ov$Slide, ov$Method_id), , drop = FALSE]
if (nrow(missing_ov)) {
  add <- missing_ov %>%
    transmute(Slide, Configuration = unname(CONFIG[Method_id]), Method_id, k = 20L,
              top_n_genes = ifelse(Slide %in% c("TENX115", "TENX117"), 40L, 100L),
              n_component_pairs = 190L,
              mean_pairwise_jaccard = NOTC, max_pairwise_jaccard = NOTC,
              union_of_lists = NA_integer_, k_times_top_n = NA_integer_, union_ratio = NOTC)
  ov <- rbind(ov, add)
}
ov <- ov[order(match(ov$Slide, SLIDES), match(ov$Method_id, names(CONFIG))), ]

# n_component_pairs is 190 on every row, being k(k - 1)/2 at the only k the overlap was
# computed on; k_times_top_n is k times the column beside it, and union_ratio is the union
# divided by that product. All three are left to the caption and the raw export.
ov_join <- ov %>%
  transmute(Slide, Method_id,
            Overlap_top_n        = top_n_genes,
            Overlap_mean_jaccard = mean_pairwise_jaccard,
            Overlap_max_jaccard  = max_pairwise_jaccard,
            # text, not integer: in the two consensus rows that were never computed it has to
            # say so alongside its neighbours, where an empty cell would read as inapplicable
            Overlap_union        = ifelse(is.na(union_of_lists), NOTC,
                                          as.character(union_of_lists)))

unpert <- tab$k == 20 & tab$Condition == "Full data"
S3 <- tab %>%
  left_join(mor,        by = c("Slide", "Method_id")) %>%
  left_join(rep_spread, by = c("Slide", "Method_id")) %>%
  left_join(bh,         by = c("Slide", "Method_id")) %>%
  left_join(sem,        by = c("Slide", "Method_id")) %>%
  left_join(ov_join,    by = c("Slide", "Method_id"))

# a join cannot tell "not applicable" from "absent", so blank every single-decomposition
# column outside the unperturbed k = 20 rows explicitly, then mark the two genuine gaps
sd_cols <- c("Moran_I_median", "Moran_I_min", "Moran_I_max", "frac_components_autocorrelated",
             "Moran_I_median_over_repeats", "Moran_I_median_min_over_repeats",
             "Moran_I_median_max_over_repeats", "n_components_annotated",
             "SemUniq_rank", "SemUniq_zscore",
             "Overlap_top_n", "Overlap_mean_jaccard", "Overlap_max_jaccard", "Overlap_union")
S3[!unpert, sd_cols] <- NA

S3$SemUniq_rank <- as.character(round(S3$SemUniq_rank, 4))
S3$SemUniq_zscore <- as.character(round(S3$SemUniq_zscore, 4))
# Keyed on the value being absent rather than on the slide, so the token disappears by itself
# once the cells are filled. All 25 unperturbed rows carry a number; the rule stays.
gap_sem <- unpert & is.na(S3$SemUniq_rank)
S3$SemUniq_rank[gap_sem] <- NOTC
S3$SemUniq_zscore[gap_sem] <- NOTC

S3 <- S3[order(match(S3$Slide, SLIDES), S3$k, S3$gene_excluded_pct,
               cond_rank(S3$Condition), match(S3$Method_id, names(CONFIG))), ]

# ---- deliver ----
# The overlap export is the raw form of those five columns, for checking against the source;
# it is not read by assemble_supplement.R and reaches the workbook only through the columns.
write.csv(S3, file.path(OUT_DIR, "TableS_comparison_full.csv"), row.names = FALSE)
write.csv(ov, file.path(OUT_DIR, "Table_gene_overlap.csv"), row.names = FALSE)

cat(sprintf("\nS3: %d rows x %d cols (expected 109 x 28) | overlap export: %d rows (expected 25)\n",
            nrow(S3), ncol(S3), nrow(ov)))
cat("rows by condition:\n")
print(table(S3$Condition, S3$k))
cat("\nsingle-decomposition columns filled on", sum(unpert), "rows,",
    "empty on", sum(!unpert), "by design\n")
cat('cells marked "', NOTC, '": semantic ', sum(S3$SemUniq_rank == NOTC, na.rm = TRUE),
    ', overlap ', sum(S3$Overlap_mean_jaccard == NOTC, na.rm = TRUE), "\n", sep = "")

# the five overlap columns must land on exactly the rows the other single-decomposition
# measures occupy; a key that failed to match would show up here as a smaller filled set
f_ov <- which(!is.na(S3$Overlap_mean_jaccard))
f_mor <- which(!is.na(S3$Moran_I_median))
cat(sprintf("overlap filled on %d rows, identical to the Moran rows: %s\n",
            length(f_ov), identical(f_ov, f_mor)))

# A blank inside the 25 is a fourth kind of absence the caption does not describe, and a join
# that silently loses a key produces exactly that. Every one of these columns must therefore
# hold either a number or the token on every unperturbed row; the control's three Moran
# spread columns are the one exception, since they exist for the control alone.
ctrl_only <- grep("_over_repeats$", sd_cols, value = TRUE)
u <- S3$k == 20 & S3$Condition == "Full data"   # recomputed: S3 has been reordered since unpert
holes <- vapply(setdiff(sd_cols, ctrl_only),
                function(cc) sum(is.na(S3[[cc]][u])), integer(1))
cat("unexplained blanks among the", sum(u), "unperturbed rows:", sum(holes), "\n")
if (any(holes > 0))
  for (cc in names(holes)[holes > 0])
    cat(sprintf('  %s empty on: %s\n', cc,
                paste(sprintf("%s/%s", S3$Slide[u][is.na(S3[[cc]][u])],
                              S3$Method_id[u][is.na(S3[[cc]][u])]), collapse = ", ")))
cat(sprintf("overlap list length: %s\n",
            paste(sprintf("%s %d", S3$Slide[f_ov][c(1, 6, 11, 16, 21)],
                          S3$Overlap_top_n[f_ov][c(1, 6, 11, 16, 21)]), collapse = " | ")))

# ---- acceptance: every number the manuscript prints must come back out of this file -------
chk <- function(label, got, want, tol = 5e-4) {
  ok <- !is.na(got) && abs(got - want) < tol
  cat(sprintf("  %-58s %8.4f vs %8.4f  %s\n", label, got, want, if (ok) "ok" else "MISMATCH"))
  ok }
g <- function(sl, m, cond, col) {
  v <- S3[S3$Slide == sl & S3$Method_id == m & S3$Condition == cond &
            S3$k == 20 & S3$gene_excluded_pct %in% c(0, 5), col]   # 0 or 5 skips the 15% rows
  if (length(v)) as.numeric(v[1]) else NA_real_ }
d <- function(sl, cond) g(sl, "consICA", cond, "Matched_stability") -
                        g(sl, "consICA_ntry1", cond, "Matched_stability")

cat("\n=== acceptance ===\n")
res <- c(
  chk("matched stability gain, full data, min (0.048)",
      min(vapply(SLIDES, d, numeric(1), "Full data")), 0.048),
  chk("matched stability gain, full data, max (0.060)",
      max(vapply(SLIDES, d, numeric(1), "Full data")), 0.060),
  chk("matched stability gain, perturbed, min (0.021)",
      min(c(vapply(SLIDES, d, numeric(1), "95% units"),
            vapply(SLIDES, d, numeric(1), "95% genes"))), 0.021),
  chk("matched stability gain, perturbed, max (0.072)",
      max(c(vapply(SLIDES, d, numeric(1), "95% units"),
            vapply(SLIDES, d, numeric(1), "95% genes"))), 0.072),
  chk("intra-uniqueness consICA, full data, min (0.980)",
      min(vapply(SLIDES, g, numeric(1), "consICA", "Full data", "Intra_uniqueness")), 0.980, 1e-3),
  chk("intra-uniqueness consICA, full data, max (0.989)",
      max(vapply(SLIDES, g, numeric(1), "consICA", "Full data", "Intra_uniqueness")), 0.989, 1e-3),
  chk("intra-uniqueness icafast fixed, full data, min (0.679)",
      min(vapply(SLIDES, g, numeric(1), "icafast_no_Rmat", "Full data", "Intra_uniqueness")), 0.679, 1e-3),
  chk("intra-uniqueness icafast fixed, full data, max (0.927)",
      max(vapply(SLIDES, g, numeric(1), "icafast_no_Rmat", "Full data", "Intra_uniqueness")), 0.927, 1e-3),
  chk("autocorrelated fraction, NCBI600, consICA (0.65)",
      g("NCBI600", "consICA", "Full data", "frac_components_autocorrelated"), 0.65),
  chk("autocorrelated fraction, NCBI600, icafast fixed (0.90)",
      g("NCBI600", "icafast_no_Rmat", "Full data", "frac_components_autocorrelated"), 0.90))

g15 <- function(sl, m) {
  v <- S3[S3$Slide == sl & S3$Method_id == m & S3$gene_excluded_pct == 15, "Matched_stability"]
  if (length(v)) as.numeric(v[1]) else NA_real_ }
res <- c(res,
  chk("matched stability gain at 15%, MEND40 (0.044)",
      g15("MEND40", "consICA") - g15("MEND40", "consICA_ntry1"), 0.044, 1e-3),
  chk("matched stability gain at 15%, NCBI600 (0.008)",
      g15("NCBI600", "consICA") - g15("NCBI600", "consICA_ntry1"), 0.008, 1e-3))

ann <- function(sl) g(sl, "icafast_no_Rmat", "Full data", "n_components_annotated") -
                    g(sl, "consICA", "Full data", "n_components_annotated")
cat("\nannotated-component difference, icafast fixed minus consICA:\n")
for (sl in c("NCBI776", "MEND40", "NCBI600", "TENX115", "TENX117"))
  cat(sprintf("  %-9s %2.0f\n", sl, ann(sl)))
cat("  expected:  NCBI776 0, MEND40 1, NCBI600 4, TENX115 6, TENX117 9\n")

# read out of the assembled table, not out of ov, so that the join is what is being checked
ovj <- function(sl, m) g(sl, m, "Full data", "Overlap_mean_jaccard")
cat("\ngene-list overlap:\n")
res <- c(res,
  chk("  TENX115 consICA (0.0984)",       ovj("TENX115", "consICA"), 0.0984, 1e-3),
  chk("  TENX115 icafast fixed (0.1522)", ovj("TENX115", "icafast_no_Rmat"), 0.1522, 1e-3),
  chk("  NCBI600 consICA (0.0122)",       ovj("NCBI600", "consICA"), 0.0122, 1e-3),
  chk("  NCBI600 icafast fixed (0.0351)", ovj("NCBI600", "icafast_no_Rmat"), 0.0351, 1e-3),
  nrow(S3) == 109L, ncol(S3) == 28L, length(f_ov) == 25L, identical(f_ov, f_mor),
  sum(holes) == 0L, !anyNA(S3$Source_axis), length(unique(S3$Source_axis)) == 2L)

cat("\n", if (all(res)) "all headline checks reproduce" else
        "SOME CHECKS FAILED - either the assembly or the text is wrong, do not smooth it", "\n")
