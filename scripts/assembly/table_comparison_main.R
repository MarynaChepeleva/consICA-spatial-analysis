# Table 2 of the main text, cut out of the supplementary comparison table so that the two
# cannot disagree. Run after tableS_comparison_full.R.
#
# Twenty rows: five slides by four configurations at k = 20 on the unperturbed data and the
# two 5% perturbations. CFS_As_Is is left out, being a deterministic instance of icafast
# identical to the fixed initialisation; Table S3 keeps both so that the identity is visible.
# The three conditions become three columns of matched stability; everything else is read off
# the unperturbed row.
#
# The two footnote markers follow from the numbers rather than from the method names:
#   *  a value of exactly 1 with an SD of exactly 0, so the repeats were identical
#   †  an intra-uniqueness of 1.000, which a single fit reaches on its own source axis

OUT_DIR <- "output/tables"
S3_CSV <- file.path(OUT_DIR, "TableS_comparison_full.csv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

METHODS <- c(icafast_no_Rmat   = "icafast, fixed init.",
             icafast_with_Rmat = "icafast, random init.",
             consICA_ntry1     = "consICA, single run",
             consICA           = "consICA")
SLIDES <- c("TENX115", "TENX117", "MEND40", "NCBI776", "NCBI600")
AXIS <- c("Gene signatures" = "GS", "Spatial profiles" = "SP")
COND <- c("Full data", "95% units", "95% genes")

S3 <- read.csv(S3_CSV, stringsAsFactors = FALSE)
S3 <- S3[S3$k == 20 & S3$gene_excluded_pct %in% c(0, 5) &
           S3$Method_id %in% names(METHODS) & S3$Condition %in% COND, ]

# mean +/- SD, with the three exceptions the caption describes: a zero SD is omitted, an SD
# that is nonzero but rounds below 0.001 says so, and an identical pair of repeats is starred
f3 <- function(v) formatC(round(v, 3), format = "f", digits = 3)
cell <- function(m, s, star = FALSE) {
  out <- paste0(f3(m), " $\\pm$ ", formatC(round(s, 3), format = "f", digits = 3))
  tiny <- s != 0 & round(s, 3) == 0
  out[tiny] <- paste0(f3(m[tiny]), " $\\pm <0.001$")
  bare <- s == 0
  out[bare] <- paste0(f3(m[bare]), if (star) "$^{*}$" else "")
  out
}

rows <- do.call(rbind, lapply(SLIDES, function(sl) {
  do.call(rbind, lapply(names(METHODS), function(mt) {
    r <- S3[S3$Slide == sl & S3$Method_id == mt, ]
    g <- function(cond, col){
      r[[col]][r$Condition == cond]
    }
    intra <- g("Full data", "Intra_uniqueness")
    data.frame(
      Slide  = sl,
      Method = unname(METHODS[mt]),
      Axis   = unname(AXIS[g("Full data", "Source_axis")]),
      Full   = cell(g("Full data", "Matched_stability"), g("Full data", "SD_matched"), star = TRUE),
      Units  = cell(g("95% units", "Matched_stability"), g("95% units", "SD_matched")),
      Genes  = cell(g("95% genes", "Matched_stability"), g("95% genes", "SD_matched")),
      Jacc   = cell(g("Full data", "TopGene_Jaccard"),   g("Full data", "SD_TopGene"), star = TRUE),
      Intra  = paste0(f3(intra), if (round(intra, 3) == 1) "$^{\\dagger}$" else ""),
      Moran  = f3(g("Full data", "Moran_I_median")),
      stringsAsFactors = FALSE)
  }))
}))

# a slide, configuration or condition missing from S3 would leave a cell empty and shorten
# the table rather than stop, so the count is asserted before anything is written
stopifnot(nrow(rows) == length(SLIDES) * length(METHODS))

write.csv(rows, file.path(OUT_DIR, "Table_comparison_main.csv"), row.names = FALSE)

# the slide name once per block, the method always bold, and the rule that separates blocks
tex <- unlist(lapply(seq_len(nrow(rows)), function(i) {
  first <- i == 1 || rows$Slide[i] != rows$Slide[i - 1]
  last <- i == nrow(rows) || rows$Slide[i] != rows$Slide[i + 1]
  cells <- c(if (first) sprintf("\\textbf{%s}", rows$Slide[i]) else "",
             sprintf("\\textbf{%s}", rows$Method[i]),
             unlist(rows[i, c("Axis", "Full", "Units", "Genes", "Jacc", "Intra", "Moran")]))
  line <- paste0(paste(cells, collapse = " & "), " \\\\")
  if (last && i != nrow(rows))
    c(line, "\\addlinespace[2pt]", "\\midrule", "\\addlinespace[1pt]") else line
}))
writeLines(tex, file.path(OUT_DIR, "Table_comparison_main.tex"))

cat(sprintf("Table 2: %d rows, %d slides x %d configurations\n",
            nrow(rows), length(SLIDES), length(METHODS)))
cat(sprintf("  starred cells: %d | daggered: %d\n",
            sum(grepl("\\^\\{\\*\\}", unlist(rows))), sum(grepl("dagger", rows$Intra))))
print(rows[, c("Slide", "Method", "Axis", "Full", "Intra", "Moran")], row.names = FALSE)
