# Gene and GO tables for every consICA object of the xenograft section
#
# Two long tables, one row per gene and one row per GO term, covering whatever
# objects exist. Nothing is filtered on significance and no component is
# singled out: the selection belongs downstream.

library(dplyr)

cica_dir <- "heap/cica_ntry20"
f_species <- "heap/gene_species.rds"
out_dir <- "output/ST_paper"

z_thresh <- 3    # |z| on the loadings that defines a tail
min_side <- 10   # smallest tail kept
top_go <- 20     # GO terms per tail, already ordered by significance

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---- the objects present ----

files <- list.files(cica_dir, pattern = "^cica[0-9]+_slide_[0-9]+_[0-9]+x[0-9]+\\.rds$",
                    full.names = TRUE)
objects <- data.frame(
  file = files,
  k = as.integer(sub("^cica([0-9]+)_.*$", "\\1", basename(files))),
  slide = as.integer(sub("^cica[0-9]+_slide_([0-9]+)_.*$", "\\1", basename(files))),
  stringsAsFactors = FALSE) %>%
  arrange(slide, k)

cat("objects found:\n")
print(objects %>% mutate(file = basename(file)), row.names = FALSE)

gene_species <- readRDS(f_species)


# ---- genes, every tail of every component ----
# Ordered by |z|, the same quantity the tail cut-off uses. A subset may be
# quoted downstream but not reordered, so the rank is carried explicitly.

tails_of <- function(S, ic) {
  z <- setNames(as.vector(scale(S[, ic])), rownames(S))
  lapply(c("+", "-"), function(sg) {
    g <- if (sg == "+") names(z)[z > z_thresh] else names(z)[z < -z_thresh]
    if (length(g) < min_side) return(NULL)
    g <- g[order(abs(z[g]), decreasing = TRUE)]
    data.frame(component = ic, tail = sg, rank = seq_along(g), gene = g,
               species = unname(gene_species[g]),
               loading = unname(S[g, ic]), z = unname(z[g]),
               row.names = NULL, stringsAsFactors = FALSE)
  }) %>% bind_rows()
}

# These files hold every object in cica_dir, so a component index belongs to one
# decomposition: ic.12 at k = 12 and ic.12 at k = 20 are different components
# with the same label, each at |z| > 3 in its own. Reading without filtering on
# k returns a union of unrelated tails that passes every check made on z.

genes_out <- lapply(seq_len(nrow(objects)), function(i) {
  S <- readRDS(objects$file[i])$S
  out <- lapply(colnames(S), function(ic) tails_of(S, ic)) %>% bind_rows()
  rm(S)
  gc()
  cbind(slide = objects$slide[i], k = objects$k[i], out)
}) %>% bind_rows()

write.csv(genes_out, file.path(out_dir, "export_genes.csv"), row.names = FALSE)
cat("\nexport_genes.csv:", nrow(genes_out), "rows,",
    length(unique(paste(genes_out$slide, genes_out$k))), "objects, k =",
    paste(sort(unique(genes_out$k)), collapse = " and "),
    "\nfilter on k before reading a gene list out of it\n")


# ---- GO, every component, both tails, both annotations ----
# getGO reorients a component on every call, so $pos and $neg need not match the
# sign in S. The orientation is recovered by correlation, and the tail as it
# appears in S is reported alongside.

# topGO truncates Term to forty characters, so the full name is recovered from
# GO.db by identifier. An identifier GO.db no longer carries keeps its string.
full_terms <- function(ids) {
  u <- unique(ids)
  tab <- suppressMessages(
    AnnotationDbi::select(GO.db::GO.db, keys = u, columns = "TERM",
                          keytype = "GOID"))
  unname(setNames(tab$TERM, tab$GOID)[ids])
}

go_of <- function(dec_file, go_file, annotation) {
  if (!file.exists(go_file)) return(NULL)
  S <- readRDS(dec_file)$S
  obj <- readRDS(go_file)
  out <- lapply(colnames(S), function(ic) {
    go_id <- sprintf("ic%02d", as.integer(sub("^ic\\.", "", ic)))
    if (!go_id %in% names(obj$GO$GOBP)) return(NULL)
    flip <- sign(cor(obj$S[rownames(S), ic], S[, ic]))
    lapply(c("pos", "neg"), function(tl) {
      df <- obj$GO$GOBP[[go_id]][[tl]]
      if (is.null(df) || !nrow(df)) return(NULL)
      p_col <- intersect(c("classicFisher", "pval", "p"), names(df))[1]
      head(data.frame(component = ic,
                      annotation = annotation,
                      go_tail = tl,
                      tail_in_S = if ((tl == "pos") == (flip > 0)) "+" else "-",
                      go_id = df$GO.ID,
                      term = df$Term,
                      annotated = df$Annotated,
                      significant = df$Significant,
                      expected = if ("Expected" %in% names(df)) df$Expected else NA_real_,
                      p = df[[p_col]],
                      fdr = df$FDR,
                      row.names = NULL, stringsAsFactors = FALSE), top_go)
    }) %>% bind_rows()
  }) %>% bind_rows()
  rm(S, obj)
  gc()
  out
}

go_out <- lapply(seq_len(nrow(objects)), function(i) {
  dec <- objects$file[i]
  hs <- go_of(dec, sub("\\.rds$", "_GO.rds", dec), "human")
  mm <- go_of(dec, sub("\\.rds$", "_GO_mm.rds", dec), "mouse")
  if (is.null(hs) && is.null(mm)) {
    cat("  no GO for", basename(dec), "\n")
    return(NULL)
  }
  cbind(slide = objects$slide[i], k = objects$k[i], bind_rows(hs, mm))
}) %>% bind_rows()

truncated <- grepl("\\.\\.\\.$", go_out$term)
full <- full_terms(go_out$go_id)
go_out$term <- ifelse(is.na(full), go_out$term, full)

write.csv(go_out, file.path(out_dir, "export_go.csv"), row.names = FALSE)
cat("\nexport_go.csv:", nrow(go_out), "rows,",
    length(unique(paste(go_out$slide, go_out$k))), "objects with GO, k =",
    paste(sort(unique(go_out$k)), collapse = " and "),
    "\nfilter on k before reading a term list out of it\n")
cat("terms stored truncated by topGO:", sum(truncated), "of", nrow(go_out),
    "\nrecovered from GO.db:", sum(truncated & !is.na(full)),
    "\nstill truncated, identifier absent from GO.db:",
    sum(truncated & is.na(full)), "\n")

cat("\nconsICA floors any p that topGO prints as < 1e-30 to 1e-31 before",
    "adjusting,\nso the smallest FDR values encode rank rather than strength.",
    "A repeated\nmantissa across unrelated terms is the tell.\n")
cat("\nwritten to", out_dir, "\n")
