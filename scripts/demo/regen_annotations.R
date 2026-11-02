# Regenerate the demo annotation tables from the saved consICA objects, which already
# include $GO. No consICA re-run and no getGO: readRDS plus write_annotation.
# Produces annotation_ (genes) and annotation_GO_ (long, both tails) per slide, which
# tableS_components_go.R turns into Tables S1 and S2.

source("scripts/demo/annotation_export.R")
output_folder <- "output/demo/"
heap <- "output/heap/demo/"

DEMOS <- list(
  list(rds = "cica8_SCT_TENX117_382x1005.rds",   slide = "TENX117", postfix = "SCT", K = 8L),
  list(rds = "cica20_SCT_MEND40_17260x2793.rds",  slide = "MEND40",  postfix = "SCT", K = 20L),
  list(rds = "cica20_SCT_NCBI776_18045x4992.rds", slide = "NCBI776", postfix = "SCT", K = 20L)
)

for (d in DEMOS) {
  f <- paste0(heap, d$rds)
  if (!file.exists(f)) { message("MISSING object: ", f); next }
  cica <- readRDS(f)
  if (is.null(cica$GO$GOBP)) { message(d$slide, ": no $GO$GOBP in object — skipped"); next }
  write_annotation(cica, d$slide, d$postfix, d$K, output_folder)
  cat(sprintf("  %s: annotation_ + annotation_GO_ (k%d) written\n", d$slide, d$K))
}
