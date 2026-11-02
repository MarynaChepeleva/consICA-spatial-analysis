# Table 1 of the main text, the slides analysed. Static: it describes the slides, not a
# result, so it is rebuilt only when a slide is added or its filtering changes.
#
# Descriptive fields and the pre-filter dimensions come from the HEST per-slide JSON in
# data/metadata; the post-filter ones are the verified decomposition sizes. Tissue names
# are set in TUMOUR below rather than read from the metrics CSVs.

suppressPackageStartupMessages(library(jsonlite))

# The HEST JSON of a slide is located by search rather than by a fixed folder: this
# repository keeps them flat under data/metadata, and download_hest.R nests them by
# organ under data/HEST.
meta_file <- function(slide) {
  roots <- unique(c(Sys.getenv("HEST_ROOT", "data/HEST"), "data"))
  f <- unlist(lapply(roots, function(r)
       list.files(r, pattern = paste0("^", slide, "\\.json$"),
                  recursive = TRUE, full.names = TRUE)))
  if (!length(f)) stop("metadata not found for ", slide, "; run download_hest.R")
  f[1]
}
out_dir <- "output/tables/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

SLIDES <- c("TENX115", "TENX117", "MEND40", "NCBI776", "NCBI600")

TUMOUR <- c(TENX115 = "cutaneous melanoma",
            TENX117 = "cutaneous melanoma",
            MEND40 = "cutaneous squamous cell carcinoma",
            NCBI776 = "breast invasive ductal carcinoma",
            NCBI600 = "muscle-invasive bladder carcinoma")
ANALYSED <- c(TENX115 = "282 × 1568", TENX117 = "382 × 1005", MEND40 = "17260 × 2793",
              NCBI776 = "18045 × 4992", NCBI600 = "14437 × 832") # genes × units, post-QC

# HEST JSONs use bare NaN (invalid JSON) for missing fields -> jsonlite rejects it;
# read as text, turn NaN into null so missing fields arrive as NULL.
read_hest <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = " ")
  fromJSON(gsub("\\bNaN\\b", "null", txt))
}
nz <- function(x) { # NULL / NA / blank -> NA_character_
  if (is.null(x) || length(x) != 1) return(NA_character_)
  x <- trimws(as.character(x))
  if (is.na(x) || x == "" || x == "NaN") NA_character_ else x
}
# write UTF-8 with a BOM so Excel opens ×, µ and — correctly (not CP1252 mojibake)
write_csv_bom <- function(df, path) {
  tmp <- tempfile()
  write.csv(df, tmp, row.names = FALSE, fileEncoding = "UTF-8")
  raw <- readBin(tmp, "raw", file.info(tmp)$size)
  unlink(tmp)
  con <- file(path, "wb")
  writeBin(c(as.raw(c(0xEF, 0xBB, 0xBF)), raw), con)
  close(con)
}

tab1 <- do.call(rbind, lapply(SLIDES, function(s) {
  j <- read_hest(meta_file(s))
  is_xen <- identical(nz(j$st_technology), "Xenium")
  platform <- nz(j$st_technology)
  if (grepl("CytAssist", paste(nz(j$subseries), nz(j$dataset_title)), ignore.case = TRUE)){
    platform <- "Visium CytAssist"
  }
  # the deposit HEST took the slide from
  study <- if (is_xen) nz(j$download_page_link1) else nz(j$study_link)
  data.frame(
    Slide = s,
    Platform = platform,
    Tumour_type = unname(TUMOUR[s]),
    As_deposited = sprintf("%d × %d", as.integer(j$nb_genes), as.integer(j$spots_under_tissue)),
    Analysed = unname(ANALYSED[s]),
    Unit_size = "100 µm",
    Segmented_cells = if (is_xen && !is.na(nz(j$num_cells)))
                        format(as.integer(j$num_cells), big.mark = " ") else "—",
    Source_link = if (is.na(study)) "—" else study,
    stringsAsFactors = FALSE, row.names = NULL)
}))
tab1 <- tab1[match(SLIDES, tab1$Slide), ]
write_csv_bom(tab1, file.path(out_dir, "Table_slides.csv"))
cat("\nTable 1 (", nrow(tab1), " rows):\n", sep = "")
print(tab1)
