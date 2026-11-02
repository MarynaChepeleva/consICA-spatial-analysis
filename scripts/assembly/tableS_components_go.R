# Two supplementary spreadsheets, one sheet each; no LaTeX output, no page-width truncation.
#   TableS_components.xlsx   one row per component, all three demonstration slides
#   TableS_go_all.xlsx       every slide, component, tail and term
# Table numbers move around, so file names are content-based, never S{n}.
# Reads only, over output/.
#
# The xenograft half of the GO table is written by the xenograft scripts. Absent, it is
# skipped and the table is the human slides alone, which the closing message states.

suppressPackageStartupMessages({
  library(dplyr)
  library(GO.db)
  library(AnnotationDbi) })
have_xlsx <- requireNamespace("openxlsx", quietly = TRUE)

demo_dir <- "output/demo/"
out_dir <- "output/tables/"
# xenograft deliverables, written by the xenograft scripts; also checked in out_dir if copied
xeno_dir <- "output/ST_paper/"
find_xeno <- function(fname) {
  for (p in c(file.path(out_dir, fname), file.path(xeno_dir, fname))) if (file.exists(p)) return(p)
  NA_character_
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat("Working dir :", getwd(), "\nOutput dir  :", normalizePath(out_dir, mustWork = FALSE), "\n")
cat("openxlsx installed:", have_xlsx, "\n")

SLIDES <- c("TENX115", "TENX117", "MEND40", "NCBI776", "NCBI600")
# Platform sets which GO statistic the component table quotes, as Methods describes: Xenium
# takes the nominal classicFisher, written "p = ", and Visium the adjusted "FDR = ".
DEMOS <- list(list(slide = "TENX117", postfix = "SCT", K = 8L,  platform = "Xenium"),
              list(slide = "MEND40",  postfix = "SCT", K = 20L, platform = "Visium"),
              list(slide = "NCBI776", postfix = "SCT", K = 20L, platform = "Visium CytAssist"))

# XLSX one sheet; fall back to CSV (announced) if no openxlsx
write_full <- function(df, name) {
  if (have_xlsx) {
    openxlsx::write.xlsx(df, file.path(out_dir, name))
    cat(sprintf("  [xlsx] %s  (%d rows x %d cols)\n", name, nrow(df), ncol(df)))
  } else {
    alt <- sub("\\.xlsx$", ".csv", name)
    write.csv(df, file.path(out_dir, alt), row.names = FALSE)
    message("openxlsx not installed -> wrote ", alt, " (", nrow(df), " rows x ", ncol(df), " cols) instead of ", name)
  }
}

# ---- component table ----
# One file, all three demonstration slides (TENX117 k8, MEND40 k20,
# NCBI776 k20), one readable row per IC, 48 rows sorted by slide then IC.
# Columns: Slide, Platform, k, IC, Morans_I, perm_p, cor_depth, top_pos, GO_pos, top_neg, GO_neg.
# GO_pos / GO_neg each hold the five leading terms of that tail in file (rank) order, "; "-
# separated, the platform's statistic named in each parenthesis: Xenium -> "(p = ...)" from
# nominal classicFisher, Visium -> "(FDR = ...)". Signatures were oriented to the positive
# tail for enrichment, so GO_neg comes from that file's tail == "neg" rows. GO.ID is not
# carried here (readable view; the machine-readable one is TableS_go_all). Values are written
# as returned — classicFisher/FDR read as text so nothing is reformatted; no analytic moran_p.
# full GO term names from GO.db (topGO truncates long Terms via numChar); NA -> keep source
full_term <- function(ids) {
  ids <- as.character(ids)
  uid <- unique(ids)
  map <- setNames(rep(NA_character_, length(uid)), uid)
  ok <- uid[uid %in% AnnotationDbi::keys(GO.db::GOTERM)]
  if (length(ok)) {
    tt <- AnnotationDbi::Term(GO.db::GOTERM[ok])
    map[names(tt)] <- tt
  }
  unname(map[ids])
}
fmt2 <- function(v) { # 2 significant figures for a statistic inside a cell
  n <- suppressWarnings(as.numeric(v))
  ifelse(is.na(n), v, formatC(n, format = "e", digits = 1))
}
# five leading terms of one tail as "Term (stat = value); ...", in file (rank) order
go_cell <- function(ic, go, tl, stat_col, stat_lab, n = 5L) {
  sub <- go[go$IC == ic & go$tail == tl, ]
  sub <- head(sub[order(sub$rank), ], n)
  if (!nrow(sub)) return("")
  paste(sprintf("%s (%s = %s)", sub$Term, stat_lab, fmt2(sub[[stat_col]])), collapse = "; ")
}

comp <- list()
for (d in DEMOS) {
  gene_f <- sprintf("%sannotation_%s_%s_k%d.csv", demo_dir, d$slide, d$postfix, d$K)
  go_f <- sprintf("%sannotation_GO_%s_%s_k%d.csv", demo_dir, d$slide, d$postfix, d$K)
  ics_f <- sprintf("%sics_stat_%s_%s_k%d.csv", demo_dir, d$slide, d$postfix, d$K)
  miss <- c(gene_f, go_f, ics_f)[!file.exists(c(gene_f, go_f, ics_f))]
  if (length(miss)) {
    message(sprintf("  %s components skipped — missing: %s  [run regen_annotations.R]",
                    d$slide, paste(basename(miss), collapse = ", ")))
    next
  }
  is_xen <- d$platform == "Xenium"
  stat_col <- if (is_xen) "classicFisher" else "FDR"
  stat_lab <- if (is_xen) "p" else "FDR"

  genes <- read.csv(gene_f, stringsAsFactors = FALSE, check.names = FALSE)
  go <- read.csv(go_f, stringsAsFactors = FALSE, check.names = FALSE,
                    colClasses = c(classicFisher = "character", FDR = "character"))  # keep as returned
  ics <- read.csv(ics_f, stringsAsFactors = FALSE, check.names = FALSE)
  ft <- full_term(go$GO.ID)
  go$Term <- ifelse(is.na(ft), go$Term, ft) # restore truncated names

  tab <- merge(genes[, c("IC", "top_pos", "top_neg")],
               ics[, c("IC", "moran_I", "perm_p", "depth_cor")], by = "IC", all.x = TRUE)
  tab$GO_pos <- vapply(tab$IC, go_cell, character(1), go = go, tl = "pos", stat_col = stat_col, stat_lab = stat_lab)
  tab$GO_neg <- vapply(tab$IC, go_cell, character(1), go = go, tl = "neg", stat_col = stat_col, stat_lab = stat_lab)
  tab$Slide <- d$slide
  tab$Platform <- d$platform
  tab$k <- d$K
  names(tab)[match(c("moran_I", "depth_cor"), names(tab))] <- c("Morans_I", "cor_depth")

  ord <- c("Slide", "Platform", "k", "IC", "Morans_I", "perm_p", "cor_depth",
           "top_pos", "GO_pos", "top_neg", "GO_neg")
  comp[[d$slide]] <- tab[, ord]
}

# one file for all three demonstration slides, sorted by slide then IC number
if (length(comp)) {
  allc <- do.call(rbind, comp)
  allc$Slide <- factor(allc$Slide, c("TENX117", "MEND40", "NCBI776"))
  allc <- allc[order(allc$Slide, as.integer(sub("IC", "", allc$IC))), ]
  allc$Slide <- as.character(allc$Slide)
  write_full(allc, "TableS_components.xlsx")
}

# ---- GO enrichment, all slides ----
# Every component, tail and term, with both p statistics as returned.
# Human part from output/demo/annotation_GO_*.csv; the xenograft GO (Slide XG6, human/mouse
# annotation_db) is produced elsewhere and appended if output/tables/TableS_go_xenograft.csv
# is present, so the whole paper's enrichment lives in one table.
common13 <- c("Slide", "k", "IC", "tail", "annotation_db", "rank", "GO.ID", "Term",
              "Annotated", "Significant", "Expected", "classicFisher", "FDR")
go_files <- list.files(demo_dir, "^annotation_GO_.*\\.csv$", full.names = TRUE)
xeno_found <- NA
if (length(go_files)) {
  human <- do.call(rbind, lapply(go_files, function(f) {
    m <- regmatches(basename(f), regexec("annotation_GO_(.+)_([^_]+)_k(\\d+)\\.csv", basename(f)))[[1]]
    d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    d$Slide <- m[2]
    d$k <- as.integer(m[4])
    d$annotation_db <- "human"
    d[, common13]
  }))
  ft <- full_term(human$GO.ID)
  human$Term <- ifelse(is.na(ft), human$Term, ft)  # restore full names
  human <- human %>% mutate(Slide = factor(Slide, SLIDES)) %>%
    arrange(Slide, k, as.integer(sub("IC", "", IC)), tail, rank) %>%
    mutate(Slide = as.character(Slide))
  # append the xenograft GO half, delivered with full term names, as-is
  xg_path <- find_xeno("TableS_go_xenograft.csv")
  xeno_found <- !is.na(xg_path)
  s_go <- if (xeno_found)
    rbind(human, read.csv(xg_path, stringsAsFactors = FALSE, check.names = FALSE)[, common13]) else human
  write_full(s_go, "TableS_go_all.xlsx")
} else message("  GO table skipped — no annotation_GO_*.csv yet [run regen_annotations.R]")

cat("\nDone. Component table + GO-all in", out_dir,
    "\nxenograft GO appended:", isTRUE(xeno_found), "\n")
