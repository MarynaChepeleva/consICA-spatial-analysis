# Every number quoted in the xenograft section, consICA at k = 12
#
# Species comes from heap/gene_species.rds, which records which of the two
# Space Ranger matrices each row was merged from. Run xeno_data_build.R first.

library(dplyr)
library(mgcv)

# IRanges and AnnotationDbi, attached by consICA, mask the dplyr select() and
# slice() generics, so the dplyr versions are bound here.
select <- dplyr::select
slice <- dplyr::slice
rename <- dplyr::rename


# ---- settings ----

cica_dir <- "heap/cica_ntry20"
f_counts <- "heap/counts_mat.rds"
f_seurat <- "heap/merged_seurat.rds"
f_species <- "heap/gene_species.rds"
f_tumour <- "data/tumor_spots.csv"
out_dir <- "output/ST_paper"

ncomp <- 12
ref_slide <- 6
HAVE_SLIDE1_HM <- FALSE   # have human alignment for slide 1
rep_slides <- if (HAVE_SLIDE1_HM) c(6, 5, 7, 8, 1) else c(6, 5, 7, 8)   # 1 is the tumour-free control
qc_slides <- if (HAVE_SLIDE1_HM) c(1, 2, 3, 4, 5, 6, 7, 8) else c(5, 6, 7, 8)

z_thresh <- 3
min_side <- 10
n_perm <- 1000
seed <- 2026
ntry <- 20

# the null takes about a minute per side, so the full table is optional
null_every_side <- TRUE

# the three components are identified by behaviour further down, never by index:
# a rerun renumbers them

# melanocytic: the set used in the manuscript, which keeps this section
# comparable with the others. recognised: the melanocytic and matrix-remodelling
# genes in the human tail of ic.12 of the reference slide. Both are reported.
melanocytic <- c("PMEL", "TYR", "TYRP1", "MLANA", "SOX10", "TRPM1", "DCT")
recognised <- c("MLANA", "PRAME", "SOX10", "S100B", "ERBB3", "CDK4", "POSTN",
                 "THBS2", "LOXL4", "SERPINE2", "CTHRC1", "TIMP2", "SPP1",
                 "FABP7")

# Used both to count and to pick the component, and the two cannot be separated:
# changing the list changes which component the rule selects.
markers <- melanocytic


# ---- functions ----

get_sides <- function(S, ic) {
  s <- as.vector(scale(S[, ic]))
  list("+" = rownames(S)[s > z_thresh],
       "-" = rownames(S)[s < -z_thresh])
}

# naive orthologue map, GAPDH to Gapdh; renamed pairs such as HBB and Hbb-bs
# are missed, so the paired share is a slight underestimate
to_mouse <- function(g) paste0(substr(g, 1, 1), tolower(substr(g, 2, nchar(g))))

pct_human <- function(genes) {
  if (!length(genes)) return(NA)
  100 * mean(genes %in% human_genes)
}

# of the human genes of a side, the share whose mouse counterpart is in the same
# side; the counterpart has to be a gene that really came from the mouse matrix
pct_paired <- function(genes) {
  h <- intersect(genes, human_genes)
  if (!length(h)) return(NA)
  100 * mean(to_mouse(h) %in% intersect(genes, mouse_genes))
}

# Human genes are detected in fewer spots, so the reference for a side is a draw
# from the same detection strata, not the flat share of human genes. Detection is
# a proxy for the selection, so the null is a floor, not an exact reference. It
# is drawn from the stream seeded once at the top, so two blocks asking for the
# same tail differ. Everything delivered reads its null from export_metrics.csv.
null_human <- function(genes, det_bin) {
  obs <- mean(genes %in% human_genes)
  idx <- match(genes, names(det_bin))
  draws <- replicate(n_perm, {
    picked <- unlist(lapply(split(idx, det_bin[idx]), function(i)
      sample(which(det_bin == det_bin[i[1]]), length(i))))
    mean(names(det_bin)[picked] %in% human_genes)
  })
  c(mean = 100 * mean(draws), sd = 100 * sd(draws),
    p = (1 + sum(draws >= obs)) / (n_perm + 1))
}

side_table <- function(S) {
  lapply(colnames(S), function(ic) {
    sides <- get_sides(S, ic)
    data.frame(ic = ic,
               side = names(sides),
               n = lengths(sides),
               n_human = sapply(sides, function(g) sum(g %in% human_genes)),
               human = sapply(sides, pct_human),
               n_paired = sapply(sides, function(g) {
                 h <- intersect(g, human_genes)
                 sum(to_mouse(h) %in% intersect(g, mouse_genes))
               }),
               paired = sapply(sides, pct_paired),
               markers = sapply(sides, function(g) sum(markers %in% g)),
               row.names = NULL,
               stringsAsFactors = FALSE)
  }) %>%
    bind_rows() %>%
    filter(n >= min_side)
}

# p floors at 1/(n_perm + 1) and is not extrapolated through a normal tail,
# which above z = 10 gives unsupportable numbers. Report z for the size of the
# departure and p as "< 0.001" on the floor.
add_null <- function(tab, S, det_bin) {
  stats <- lapply(seq_len(nrow(tab)), function(i) {
    null_human(get_sides(S, tab$ic[i])[[tab$side[i]]], det_bin)
  })
  tab$null <- sapply(stats, "[", "mean")
  tab$z <- (tab$human - tab$null) / sapply(stats, "[", "sd")
  tab$p_perm <- sapply(stats, "[", "p")
  tab$fdr <- p.adjust(tab$p_perm, "BH")
  tab
}

# The tail carrying the melanocytic programme, by one rule applied everywhere.
# Not "most markers": a larger tail collects more of any list, which here picks
# the abundance axis. The list filters, the human share orders. Where no tail
# holds a gene of the list, the most human tail is returned.
pick_programme <- function(tab) {
  hits <- tab %>% filter(markers > 0)
  if (!nrow(hits)) hits <- tab
  hits %>% arrange(desc(human)) %>% dplyr::slice(1)
}

cica_path <- function(slide, k = ncomp) {
  hits <- list.files(cica_dir, full.names = TRUE,
                     pattern = sprintf("^cica%d_slide_%d_[0-9]+x[0-9]+\\.rds$",
                                       k, slide))
  if (!length(hits)) stop("no decomposition for slide ", slide)
  hits[1]
}

detection_bins <- function(counts) {
  det <- Matrix::rowSums(counts > 0)
  setNames(cut(rank(det, ties.method = "first"), breaks = 20, labels = FALSE),
           names(det))
}

# every table that reaches the manuscript is written out as well as printed, so
# a number can be traced to a file instead of to a console that has scrolled
scalars <- list()
keep <- function(name, value) {
  scalars[[name]] <<- value
  value
}
save_table <- function(x, name) {
  write.csv(as.data.frame(x), file.path(out_dir, paste0(name, ".csv")),
            row.names = FALSE)
  print(as.data.frame(x), digits = 4, row.names = FALSE)
  invisible(x)
}


# ---- species and data ----
set.seed(seed)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

counts_mat <- readRDS(f_counts)
gene_species <- readRDS(f_species)[rownames(counts_mat)]
stopifnot(!anyNA(gene_species))
human_genes <- names(gene_species)[gene_species == "human"]
mouse_genes <- names(gene_species)[gene_species == "mouse"]

cat("=== gene universe ===\n")
keep("genes_human", length(human_genes))
keep("genes_mouse", length(mouse_genes))
keep("genes_human_pct", 100 * length(human_genes) / length(gene_species))
cat("human", length(human_genes),
    sprintf("(%.1f%%)", scalars$genes_human_pct),
    " mouse", length(mouse_genes), "\n")

# the melanocytic set is only usable to the extent its genes survived the
# detection filter, so state which of them are present before relying on it
cat("\nmelanocytic set in the matrix:",
    paste(intersect(melanocytic, rownames(counts_mat)), collapse = " "), "\n")
cat("absent:",
    paste(setdiff(melanocytic, rownames(counts_mat)), collapse = " "), "\n")
cat("recognised set present:",
    length(intersect(recognised, rownames(counts_mat))), "of",
    length(recognised), "\n")

# Genes absent from the matrix are dropped from the counting set.
markers <- intersect(markers, rownames(counts_mat))
cat("counting against", length(markers), "genes:",
    paste(markers, collapse = " "), "\n")

merged_seurat <- readRDS(f_seurat)
meta <- merged_seurat@meta.data
rm(merged_seurat)
gc()

tumour_barcodes <- as.character(read.csv(f_tumour)[[1]])
spots_of <- function(s) rownames(meta)[meta$sample_id == s]


# ---- sequencing depth and human read fraction per slide ----
# Count medians are rounded at source, so the section table and the depth audit
# carry one value.

cat("\n=== depth and human read fraction per slide ===\n")
qc <- lapply(qc_slides, function(s) {
  sp_s <- spots_of(s)
  if (!length(sp_s)) return(NULL)
  cs <- counts_mat[, sp_s]
  total <- Matrix::colSums(cs)
  hum <- Matrix::colSums(cs[human_genes, , drop = FALSE])
  keep <- total > 0
  frac <- (hum / total)[keep]
  data.frame(slide = s,
             spots = ncol(cs),
             counts_median = round(median(total)),
             genes_median = round(median(Matrix::colSums(cs > 0))),
             human_min = min(frac),
             human_median = median(frac),
             human_max = max(frac),
             spots_over_40 = sum(frac > 0.40),
             stringsAsFactors = FALSE)
}) %>% bind_rows()
save_table(qc, "slide_qc")


# ---- reference slide ----

cica <- readRDS(cica_path(ref_slide))
S <- cica$S
M <- cica$M
spots <- colnames(M)
is_tumour <- spots %in% tumour_barcodes

counts_ref <- counts_mat[, spots]
det_bin <- detection_bins(counts_ref)
human_frac <- Matrix::colSums(counts_ref[human_genes, ]) /
  Matrix::colSums(counts_ref)

cat("\n=== reference slide", ref_slide, "===\n")
keep("ref_spots", length(spots))
keep("ref_tumour_spots", sum(is_tumour))
keep("ref_human_frac_median", median(human_frac))
keep("ref_human_frac_max", max(human_frac))
keep("ref_spots_over_40", sum(human_frac > 0.40))
cat("spots", length(spots), " tumour-annotated", sum(is_tumour), "\n")

comp <- side_table(S) %>% arrange(desc(human))
comp <- if (null_every_side) {
  add_null(comp, S, det_bin)
} else {
  add_null(comp[seq_len(4), ], S, det_bin) %>% bind_rows(comp[-seq_len(4), ])
}
cat("\ncomposition of every side against a detection-matched null\n")
save_table(comp, "composition_reference")

# What pooling both tails of a component would report. Taken over every gene at
# |z| > 3, so it includes the tails that did not clear min_side.
pooled <- lapply(colnames(S), function(ic) {
  sig <- rownames(S)[abs(scale(S[, ic])) > z_thresh]
  data.frame(ic = ic, pooled_human = pct_human(sig), stringsAsFactors = FALSE)
}) %>% bind_rows()
cat("\nthe same components with the two tails pooled\n")
save_table(comp %>% select(ic, side, human) %>% left_join(pooled, by = "ic"),
           "pooled_reference")


# ---- which component plays which role ----
# Identified by behaviour, because a rerun renumbers the components. Programme:
# most markers. Abundance: best explained by the human read fraction. Orthologue
# pairs: the most human of the remaining heavily paired sides.

dev_all <- sapply(rownames(M), function(ic)
  summary(gam(M[ic, spots] ~ s(human_frac, bs = "cs"),
              method = "REML"))$dev.expl)

prog_row <- pick_programme(comp)
abun_row <- comp %>% filter(ic == names(which.max(dev_all))) %>%
  arrange(desc(human)) %>% dplyr::slice(1)
shared_row <- comp %>%
  filter(!ic %in% c(prog_row$ic, abun_row$ic), paired > 60) %>%
  arrange(desc(human)) %>% dplyr::slice(1)

# the two must not be the same component; if they are, the rules disagree with
# the data and that has to be looked at rather than reported
if (prog_row$ic == abun_row$ic)
  warning("programme and abundance landed on the same component: ",
          prog_row$ic, call. = FALSE)

roles <- bind_rows(mutate(prog_row, role = "programme"),
                   mutate(abun_row, role = "abundance"),
                   mutate(shared_row, role = "orthologue pairs")) %>%
  select(role, everything())

ic_prog <- prog_row$ic
side_prog <- prog_row$side
ic_abun <- abun_row$ic

cat("\ncomponents identified by role\n")
save_table(roles, "components")


# Ordered gene lists and GO tables are written by xeno_export.R across every
# object. The reference section is its slide 6, k = 12 rows.


# ---- redundancy with the per-spot human read fraction ----
redundancy <- rbind(
  data.frame(ic = ic_abun, weight = M[ic_abun, spots], frac = human_frac),
  data.frame(ic = ic_prog, weight = M[ic_prog, spots], frac = human_frac))

gam_fit <- redundancy %>%
  group_by(ic) %>%
  summarise(dev = summary(gam(weight ~ s(frac, bs = "cs"),
                              method = "REML"))$dev.expl,
            r2 = cor(weight, frac)^2,
            rho = cor(weight, frac, method = "spearman"),
            .groups = "drop")

cat("\n=== how much of a component the human read fraction accounts for ===\n")
cat("mgcv::gam(weight ~ s(frac, bs = \"cs\"), method = \"REML\")\n")
save_table(gam_fit, "redundancy_reference")


# ---- does the result depend on k ----

cat("\n=== the same slide at k = 20 ===\n")
S20 <- readRDS(cica_path(ref_slide, 20))$S
save_table(side_table(S20) %>%
             arrange(desc(human)) %>%
             head(6) %>%
             add_null(S20, det_bin),
           "ksensitivity_reference")


# ---- control and replication ----
# The component is chosen on each section by the same rule and never named in
# advance: most markers first, then highest human share. On the control the rule
# finds no markers and returns its most human side.

cat("\n=== control and replication at k =", ncomp, "===\n")
replication <- lapply(rep_slides, function(s) {
  obj <- readRDS(cica_path(s))
  sp <- colnames(obj$M)
  tab <- pick_programme(side_table(obj$S)) %>%
    add_null(obj$S, detection_bins(counts_mat[, sp]))
  ref_side <- get_sides(S, ic_prog)[[side_prog]]
  own_side <- get_sides(obj$S, tab$ic)[[tab$side]]
  rm(obj)
  gc()
  cbind(slide = s, tab,
        shared = length(intersect(own_side, intersect(ref_side, human_genes))),
        jaccard = length(intersect(own_side, ref_side)) /
          length(union(own_side, ref_side)),
        melanocytic_genes = paste(intersect(melanocytic, own_side),
                                  collapse = " "),
        recognised_genes = paste(intersect(recognised, own_side),
                                 collapse = " "))
}) %>% bind_rows()
save_table(replication, "replication")


# The main replication table is assembled by xeno_tables.R, so that its expected
# percentages come from the same draw of the null as the supplementary workbook.


# ---- every tail of every object, so the selection can be audited ----
# The rule filters on markers before ordering by human share, so a tail could be
# purer than the selected one and carry none. Every object in cica_dir is here,
# whatever its k, so a row is identified by slide and k together.

all_files <- list.files(cica_dir, full.names = TRUE,
                        pattern = "^cica[0-9]+_slide_[0-9]+_[0-9]+x[0-9]+\\.rds$")

metrics <- lapply(all_files, function(f) {
  kk <- as.integer(sub("^cica([0-9]+)_.*$", "\\1", basename(f)))
  sl <- as.integer(sub("^cica[0-9]+_slide_([0-9]+)_.*$", "\\1", basename(f)))
  obj <- readRDS(f)
  out <- side_table(obj$S) %>%
    arrange(desc(human)) %>%
    add_null(obj$S, detection_bins(counts_mat[, colnames(obj$M)]))
  rm(obj)
  gc()
  cbind(slide = sl, k = kk, out)
}) %>% bind_rows()

cat("\nevery tail of every object, k =",
    paste(sort(unique(metrics$k)), collapse = " and "), "\n")
save_table(metrics, "export_metrics")

cat("\ntails purer than the selected one on the same section\n")
audit <- metrics %>%
  filter(k == ncomp) %>%
  inner_join(replication %>%
               select(slide, sel_ic = ic, sel_side = side, sel_human = human),
             by = "slide") %>%
  filter(human > sel_human, !(ic == sel_ic & side == sel_side)) %>%
  select(slide, ic, side, n, human, paired, markers,
         sel_ic, sel_side, sel_human)
if (nrow(audit)) {
  save_table(audit, "purer_unselected_tails")
} else {
  cat("none: on every section the selected tail is also the most human\n")
}


# ---- audit: is sequencing depth a confounder ----
# Ruled out: the deepest section is the second purest and the smallest tail
# belongs to the shallowest. Depth stays in slide_qc.csv and audit_depth.csv.

cat("\naudit: depth against the selected tail\n")
save_table(replication %>%
             select(slide, ic, side, n, human, paired) %>%
             left_join(qc %>% select(slide, counts_median, genes_median),
                       by = "slide") %>%
             arrange(counts_median),
           "audit_depth")


# ---- GO of both tails of the tumour component, anchored to the loadings ----
# getGO reorients a component on every call, so $pos and $neg cannot be assumed
# to match the sign in S. The orientation is recovered by correlating the
# component in the GO object with the same component in the reference S.

go_of <- function(path, ic, label) {
  obj <- readRDS(path)
  go_id <- sprintf("ic%02d", as.integer(sub("^ic\\.", "", ic)))
  flip <- sign(cor(obj$S[rownames(S), ic], S[, ic]))
  out <- lapply(c("pos", "neg"), function(tl) {
    df <- obj$GO$GOBP[[go_id]][[tl]]
    if (is.null(df) || !nrow(df)) return(NULL)
    tail_in_S <- if ((tl == "pos") == (flip > 0)) "+" else "-"
    head(data.frame(annotation = label, go_tail = tl, s_tail = tail_in_S,
                    term = df$Term,
                    p = df[[intersect(c("classicFisher", "pval"),
                                      names(df))[1]]],
                    fdr = df$FDR,
                    stringsAsFactors = FALSE), 6)
  })
  rm(obj)
  gc()
  bind_rows(out)
}

cat("\n=== GO of the programme component, both tails, both annotations ===\n")
f_go <- c(human = sub("\\.rds$", "_GO.rds", cica_path(ref_slide)),
          mouse = sub("\\.rds$", "_GO_mm.rds", cica_path(ref_slide)))
if (all(file.exists(f_go))) {
  save_table(bind_rows(go_of(f_go[["human"]], ic_prog, "human"),
                       go_of(f_go[["mouse"]], ic_prog, "mouse")),
             "go_programme")
} else {
  cat("GO objects not present yet, skipped:\n")
  cat(paste0("  ", f_go[!file.exists(f_go)], collapse = "\n"), "\n")
}


cat("\nboth gene sets, and whether each gene survived the detection filter\n")
save_table(data.frame(
  gene = c(melanocytic, setdiff(recognised, melanocytic)),
  set = c(rep("melanocytic", length(melanocytic)),
          rep("recognised", length(setdiff(recognised, melanocytic)))),
  in_matrix = c(melanocytic, setdiff(recognised, melanocytic)) %in%
    rownames(counts_mat)),
  "marker_sets")


# ---- scalars, so single numbers are traceable to a file too ----

keep("component_programme", paste0(ic_prog, side_prog))
keep("component_abundance", ic_abun)
keep("deviance_abundance", gam_fit$dev[gam_fit$ic == ic_abun])
keep("deviance_programme", gam_fit$dev[gam_fit$ic == ic_prog])
if (HAVE_SLIDE1_HM) {
  keep("control_human_frac_max", qc$human_max[qc$slide == 1])
  keep("control_human_frac_median", qc$human_median[qc$slide == 1])
}
keep("panel_b_sides", nrow(comp))
keep("n_perm", n_perm)
keep("z_thresh", z_thresh)
keep("seed", seed)

keep("ntry", ntry)
keep("consICA_version", as.character(packageVersion("consICA")))

save_table(data.frame(name = names(scalars),
                      value = as.character(unlist(scalars)),
                      row.names = NULL),
           "key_numbers")

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo.txt"))

cat("\ntables written to", out_dir, "\n")
