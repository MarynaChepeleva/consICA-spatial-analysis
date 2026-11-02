# Annotation export shared by the two demo scripts. From a getGO-annotated consICA
# object it writes two files into output_folder:
# - annotation_GO_<slide>_<postfix>_k<K>.csv   per component, tail and term: the BP
#                                              terms, classicFisher and FDR
# - annotation_<slide>_<postfix>_k<K>.csv      per component: the top genes of both
#                                              tails, and an optional label column
#                                              filled by hand
#
# Tails match: GO uses getGO's $pos and $neg, genes the sign of cica$S. Tables come
# ordered by topGO, which the rank column keeps. consICA floors p at 1e-31 before BH,
# so for the strongest terms FDR encodes rank, not strength.

# Number of exported genes and GO terms per component tail, by slide.
ANNOT_EXPORT <- list(TENX117 = c(genes = 20, terms = 20L),
                     MEND40  = c(genes = 20L, terms = 20L),
                     NCBI776 = c(genes = 20L, terms = 20L))
ANNOT_DEFAULT <- c(genes = 20L, terms = 20L)

annot_param <- function(slide_id, what) {
  p <- ANNOT_EXPORT[[slide_id]]
  if (is.null(p)) {
    message("annotation export: ", slide_id, " is not in ANNOT_EXPORT, using ",
            ANNOT_DEFAULT[["genes"]], " genes and ", ANNOT_DEFAULT[["terms"]], " terms")
    p <- ANNOT_DEFAULT
  }
  unname(p[[what]])
}

write_annotation <- function(cica, slide_id, postfix, ncomp, output_folder,
                             top_n_genes = NULL, top_n_terms = NULL) {
  if (is.null(top_n_genes)) top_n_genes <- annot_param(slide_id, "genes")
  if (is.null(top_n_terms)) top_n_terms <- annot_param(slide_id, "terms")
  S <- cica$S

  # ---- genes: one row per component, both tails ----
  genes <- do.call(rbind, lapply(seq_len(ncomp), function(ic) {
    data.frame(
      IC      = paste0("IC", ic),
      label   = "",                                             # optional, filled by hand
      top_pos = paste(head(names(sort(S[, ic], decreasing = TRUE)), top_n_genes), collapse = ", "),
      top_neg = paste(head(names(sort(S[, ic])),                    top_n_genes), collapse = ", "),
      stringsAsFactors = FALSE)
  }))
  write.csv(genes,
            paste0(output_folder, "annotation_", slide_id, "_", postfix, "_k", ncomp, ".csv"),
            row.names = FALSE)

  # ---- GO: one row per component, tail and term, BP only ----
  keep_cols <- c("GO.ID", "Term", "Annotated", "Significant", "Expected",
                 "classicFisher", "FDR")
  go_long <- do.call(rbind, lapply(seq_len(ncomp), function(ic) {
    node <- cica$GO$GOBP[[sprintf("ic%02d", ic)]]
    do.call(rbind, lapply(c("pos", "neg"), function(tl) {
      tab <- node[[tl]]
      if (is.null(tab) || !nrow(tab)) return(NULL)
      tab <- head(tab, top_n_terms)
      cols <- intersect(keep_cols, names(tab))
      data.frame(IC = paste0("IC", ic), tail = tl, rank = seq_len(nrow(tab)),
                 tab[, cols, drop = FALSE],
                 stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
    }))
  }))
  write.csv(go_long,
            paste0(output_folder, "annotation_GO_", slide_id, "_", postfix, "_k", ncomp, ".csv"),
            row.names = FALSE)

  invisible(list(genes = genes, go = go_long))
}
