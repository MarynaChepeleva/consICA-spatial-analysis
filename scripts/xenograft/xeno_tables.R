# Manuscript-facing tables for the xenograft section
#
# The computing scripts write their tables in full precision under internal
# names. This one converts them for delivery: sections become XG6 and components
# IC12, matching Table 1 and the prose, percentages are rounded to the precision
# the text quotes, and the genes and enriched terms of a tail are folded into the
# row that describes it.
#
# Tail composition is recomputed here rather than read, because min_side keeps
# tails under ten genes out of every export. The threshold is still respected:
# z and the permutation p stay empty on such a row. The recomputation is checked
# against export_metrics on every row the two share.

library(dplyr)

select <- dplyr::select
slice <- dplyr::slice
rename <- dplyr::rename

cica_dir <- "heap/cica_ntry20"
f_species <- "heap/gene_species.rds"
out_dir <- "output/ST_paper"

ref_slide <- 6
ref_k <- 12
HAVE_SLIDE1_HM <- FALSE   # have human alignment for slide 1
slide_order <- if (HAVE_SLIDE1_HM) c(6, 5, 7, 8, 1) else c(6, 5, 7, 8)
top_genes_n <- 15
top_terms_n <- 5
z_thresh <- 3
pct_digits <- 1
z_digits <- 2

cica_path <- function(slide, k = ref_k) {
  hits <- list.files(cica_dir, full.names = TRUE,
                     pattern = sprintf("^cica%d_slide_%d_[0-9]+x[0-9]+\\.rds$",
                                       k, slide))
  if (!length(hits)) stop("no decomposition for slide ", slide)
  hits[1]
}

# Identifiers are converted only on the way out, so the internal names used
# throughout the computing scripts survive a change of manuscript identifier.
xg_label <- function(slide) paste0("XG", slide)
ic_label <- function(component) paste0("IC", sub("^ic\\.", "", component))
ic_number <- function(component) as.integer(sub("^ic\\.", "", component))
tail_label <- function(side) ifelse(side == "+", "pos", "neg")

read_out <- function(name) {
  read.csv(file.path(out_dir, paste0(name, ".csv")), stringsAsFactors = FALSE)
}

deliver <- function(x, file) {
  if (grepl("\\.xlsx$", file)) {
    openxlsx::write.xlsx(list(xenograft = x), file.path(out_dir, file))
  } else {
    write.csv(x, file.path(out_dir, file), row.names = FALSE)
  }
  cat(sprintf("%-32s %6d rows %3d columns\n", file, nrow(x), ncol(x)))
}


# ---- sources ----
metrics <- read_out("export_metrics")
gene_export <- read_out("export_genes")
go <- read_out("export_go")

wanted <- c("go_id", "term", "annotated", "significant", "expected", "p", "fdr")
absent <- setdiff(wanted, names(go))
if (length(absent)) {
  stop("export_go.csv lacks ", paste(absent, collapse = ", "),
       ": run xeno_export.R first")
}

base <- metrics %>% filter(k == ref_k)

gene_species <- readRDS(f_species)
human_genes <- names(gene_species)[gene_species == "human"]
mouse_genes <- names(gene_species)[gene_species == "mouse"]

# the set the marker count is taken over, read rather than restated so that the
# two scripts cannot drift apart
markers <- read_out("marker_sets") %>%
  filter(set == "melanocytic") %>%
  pull(gene)

# naive orthologue map, GAPDH to Gapdh; renamed pairs such as HBB and Hbb-bs are
# missed, so the paired share is a slight underestimate
to_mouse <- function(g) paste0(substr(g, 1, 1), tolower(substr(g, 2, nchar(g))))


# ---- every tail of every component, with no threshold on its length ----
tail_row <- function(g, z) {
  h <- intersect(g, human_genes)
  ordered <- head(g[order(abs(z[g]), decreasing = TRUE)], top_genes_n)
  data.frame(
    n = length(g),
    n_human = length(h),
    human = if (length(g)) 100 * mean(g %in% human_genes) else NA_real_,
    n_paired = sum(to_mouse(h) %in% intersect(g, mouse_genes)),
    paired = if (length(h)) {
      100 * mean(to_mouse(h) %in% intersect(g, mouse_genes))
    } else NA_real_,
    markers = sum(markers %in% g),
    top_genes = if (length(g)) {
      paste0(ordered, " (", substr(unname(gene_species[ordered]), 1, 1), ")",
             collapse = ", ")
    } else NA_character_,
    row.names = NULL, stringsAsFactors = FALSE)
}

census <- lapply(slide_order, function(s) {
  S <- readRDS(cica_path(s))$S
  out <- lapply(colnames(S), function(ic) {
    z <- setNames(as.vector(scale(S[, ic])), rownames(S))
    tails <- list("+" = names(z)[z > z_thresh], "-" = names(z)[z < -z_thresh])
    cbind(slide = s, ic = ic, side = names(tails),
          bind_rows(lapply(tails, tail_row, z = z)),
          # the union of the two tails, which is what pooling them would report
          pooled_human = 100 * mean(unlist(tails) %in% human_genes))
  }) %>% bind_rows()
  rm(S)
  gc()
  out
}) %>% bind_rows()


# ---- the recomputation must agree with what was published ----

shared <- census %>%
  inner_join(base %>% transmute(slide, ic, side, e_n = n, e_nh = n_human,
                                e_h = human, e_np = n_paired, e_p = paired,
                                e_m = markers),
             by = c("slide", "ic", "side"))
same <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) < 1e-9)
bad <- shared %>%
  filter(n != e_n | n_human != e_nh | n_paired != e_np | markers != e_m |
           !same(human, e_h) | !same(paired, e_p))
if (nrow(bad)) {
  print(as.data.frame(bad))
  stop("recomputed tails disagree with export_metrics on ", nrow(bad), " rows")
}

export_cell <- gene_export %>%
  filter(k == ref_k, rank <= top_genes_n) %>%
  group_by(slide, ic = component, side = tail) %>%
  summarise(cell = paste0(gene, " (", substr(species, 1, 1), ")",
                          collapse = ", "), .groups = "drop")
cmp <- census %>% inner_join(export_cell, by = c("slide", "ic", "side"))
stopifnot(nrow(cmp) == nrow(shared), identical(cmp$top_genes, cmp$cell))

if (file.exists(file.path(out_dir, "pooled_reference.csv"))) {
  chk <- read_out("pooled_reference") %>%
    select(ic, pooled_human) %>%
    distinct() %>%
    left_join(census %>% filter(slide == ref_slide) %>%
                transmute(ic, recomputed = pooled_human) %>% distinct(),
              by = "ic")
  stopifnot(all(abs(chk$pooled_human - chk$recomputed) < 1e-9))
}
cat("recomputed tails reproduce export_metrics, export_genes and",
    "pooled_reference on all", nrow(shared), "published rows\n")


# ---- Table S4, one row per tail of every section at k = 12 ----
# The null and z of every tail are the single draw held in export_metrics, which
# main Table 3 rests on too, so the column is one quantity throughout.
#
# Every tail gets a row, including one no gene passed, so the sheet is
# rectangular at sections x components x 2. On such a row the counts are zero and
# the shares, the draw and the genes empty. pooled_human_pct is the exception,
# taken over the union of the two tails. fdr is read from export_metrics, so the
# correction covers the same tails as before.

go_cell <- function(db) {
  go %>%
    filter(k == ref_k, annotation == db) %>%
    group_by(slide, ic = component, side = tail_in_S) %>%
    slice(seq_len(top_terms_n)) %>%
    summarise(cell = paste0(term, " (FDR = ", formatC(fdr, format = "g",
                                                      digits = 2), ")",
                            collapse = "; "),
              .groups = "drop")
}

s4 <- census %>%
  left_join(base %>% select(slide, ic, side, null, z, p_perm, fdr),
            by = c("slide", "ic", "side")) %>%
  left_join(go_cell("human") %>% rename(GO_human = cell),
            by = c("slide", "ic", "side")) %>%
  left_join(go_cell("mouse") %>% rename(GO_mouse = cell),
            by = c("slide", "ic", "side")) %>%
  arrange(match(slide, slide_order), ic_number(ic), side) %>%
  transmute(Slide = xg_label(slide), k = ref_k, IC = ic_label(ic),
            tail = tail_label(side),
            n, n_human, human_pct = round(human, pct_digits),
            n_paired, paired_pct = round(paired, pct_digits),
            pooled_human_pct = round(pooled_human, pct_digits), markers,
            null_mean = round(null, pct_digits), z = round(z, z_digits),
            p_species = p_perm, fdr_species = fdr,
            top_genes, GO_human, GO_mouse)

stopifnot(nrow(s4) >= nrow(base), !anyNA(s4$n), !anyNA(s4$pooled_human_pct))

# the sheet is rectangular: every section carries both tails of every component
per_section <- table(s4$Slide)
keys <- paste(s4$Slide, s4$IC, s4$tail)
stopifnot(nrow(s4) == length(slide_order) * ref_k * 2,
          all(per_section == ref_k * 2),
          !anyDuplicated(keys))
cat(sprintf("Table S4: %d rows, %d per section, %d of them an empty tail\n",
            nrow(s4), unique(per_section), sum(s4$n == 0)))

# The rows added above must leave the species test alone: as many rows carry an
# fdr as export_metrics has tails, and each equals the value it holds.
src <- base %>% transmute(key = paste(xg_label(slide), ic_label(ic), tail_label(side)),
                          src_fdr = fdr)
cmp <- merge(src, data.frame(key = keys, fdr_species = s4$fdr_species), by = "key")
stopifnot(nrow(cmp) == nrow(base),
          sum(!is.na(s4$fdr_species)) == nrow(base),
          all(same(cmp$src_fdr, cmp$fdr_species)))
cat(sprintf("fdr_species: %d rows carry one, exactly the %d tails export_metrics holds\n",
            sum(!is.na(s4$fdr_species)), nrow(base)),
    sprintf("            %d carry none: %d empty tails and %d under min_side\n",
            sum(is.na(s4$fdr_species)), sum(s4$n == 0), sum(s4$n > 0 & is.na(s4$z))))

deliver(s4, "TableS_xenograft.xlsx")


# ---- the main replication table, emitted rather than typed ----
# Assembled here rather than in xeno_numbers.R so that its expected percentages
# come from the same draw of the null as the workbook. Composition comes from
# export_metrics, the overlap with the reference from replication.csv. The tex
# file is the table body only, in the row format the benchmark tables use. The
# Jaccard overlap was dropped from the table and is still in replication.csv.

section_role <- c("6" = "discovery", "5" = "replicate", "7" = "replicate",
                  "8" = "replicate", "1" = "control")

tab3 <- read_out("replication") %>%
  select(slide, ic, side) %>%
  left_join(base %>% select(slide, ic, side, n, human, null, paired, markers),
            by = c("slide", "ic", "side")) %>%
  arrange(match(slide, slide_order)) %>%
  transmute(section = xg_label(slide),
            role = unname(section_role[as.character(slide)]),
            component = paste(ic_label(ic), tail_label(side)),
            genes = n,
            human_observed = formatC(human, format = "f", digits = 1),
            human_expected = formatC(null, format = "f", digits = 1),
            paired = formatC(paired, format = "f", digits = 1),
            markers)

# every selected tail must be found in export_metrics, or a cell would be empty
stopifnot(nrow(tab3) == length(slide_order), !anyNA(tab3))
deliver(tab3, "Table_replication_xenograft.csv")

mtx <- as.matrix(tab3)
# the control slide is absent unless HAVE_SLIDE1_HM; 0L then matches no row and the
# extra rule above it is simply not drawn
control_row <- match("control", tab3$role)
before_control <- if (is.na(control_row)) 0L else control_row - 1L
tex <- vapply(seq_len(nrow(mtx)), function(i)
  paste0(paste(mtx[i, ], collapse = " & "),
         if (i == before_control) " \\\\[2pt]" else " \\\\"),
  character(1))
writeLines(tex, file.path(out_dir, "Table_replication_xenograft.tex"))
cat("Table_replication_xenograft.tex        body only, ", nrow(mtx), " rows\n",
    sep = "")


# ---- the xenograft part of the paper's GO table ----
# rank is the position within one section, component, tail and database, in the
# order the source holds. tail is the tail of S, not the pos/neg label of the GO
# object, since getGO reorients a component on every call. k = 20 appears nowhere
# in the paper, so its enrichment is left out.

go_out <- go %>%
  filter(k == ref_k) %>%
  group_by(slide, k, component, tail_in_S, annotation) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  arrange(match(slide, slide_order), k, ic_number(component), tail_in_S,
          annotation, rank) %>%
  transmute(Slide = xg_label(slide), k, IC = ic_label(component),
            tail = tail_label(tail_in_S), annotation_db = annotation, rank,
            GO.ID = go_id, Term = term, Annotated = annotated,
            Significant = significant, Expected = expected,
            classicFisher = p, FDR = fdr)
deliver(go_out, "TableS_go_xenograft.csv")


# ---- report ----

cat("\nrows per section, S4 then GO:\n")
print(s4 %>% count(Slide) %>%
        left_join(go_out %>% count(Slide, name = "go"), by = "Slide"),
      row.names = FALSE)
cat("\nrows no export carries:", nrow(s4) - nrow(base),
    "- of them, tails under min_side:", sum(s4$n > 0 & is.na(s4$z)),
    "and tails with no gene at all:", sum(s4$n == 0),
    "\nrows with no null and no z, which is those same two kinds:", sum(is.na(s4$z)),
    "\nterms still ending in an ellipsis:", sum(grepl("\\.\\.\\.$", go_out$Term)),
    "\ntails with no GO against either database:",
    sum(is.na(s4$GO_human) & is.na(s4$GO_mouse)), "\n")
# The tails that carry no draw, against the tail Table 3 shows for the same
# section. A tail under min_side has genes and so a human share, which can be
# compared. A tail with no gene has no share, cannot be more human than
# anything, and is listed but not compared.
restored <- s4 %>%
  filter(is.na(z)) %>%
  select(Slide, IC, tail, n, n_human, human_pct, paired_pct, top_genes) %>%
  left_join(read_out("replication") %>%
              transmute(Slide = xg_label(slide),
                        shown_IC = ic_label(ic),
                        shown_human_pct = round(human, pct_digits)),
            by = "Slide")

cat("\ntails carrying no draw, beside the tail Table 3 shows",
    "(n = 0 is an empty tail, n > 0 one under min_side):\n")
print(as.data.frame(restored %>% select(-top_genes)), row.names = FALSE)
cat("\ntheir genes:\n")
print(as.data.frame(restored %>% select(Slide, IC, tail, top_genes)),
      row.names = FALSE)
comparable <- restored %>% filter(n > 0)
cat("\nany such tail more human than the tail shown for its section, over the",
    nrow(comparable), "with genes:",
    if (any(comparable$human_pct > comparable$shown_human_pct)) {
      "YES - the threshold has to be stated in the text"
    } else "no", "\n")
cat("\nheader and first two rows of S4:\n")
print(as.data.frame(head(s4, 2)))
cat("\nthe replication table as delivered:\n")
print(as.data.frame(tab3), row.names = FALSE)
cat("\nwritten to", out_dir, "\n")
