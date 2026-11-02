# Assemble the four supplementary tables into one submission workbook,
# output/tables/Supplementary_Tables.xlsx. Run after table_slides.R, table_comparison_main.R,
# tableS_comparison_full.R and tableS_components_go.R, when the numbering is final.
# Nothing is computed here: the source files stay content-named and this script alone
# knows the S-numbers. Formatting is display only and carries no information.

if (!requireNamespace("openxlsx", quietly = TRUE))
  stop("assemble_supplement.R needs the openxlsx package for styled multi-sheet output.")
suppressPackageStartupMessages(library(openxlsx))

out_dir <- "output/tables/"
contents_csv <- "data/supplement_contents.csv"   # the workbook's first sheet
wb_path <- file.path(out_dir, "Supplementary_Tables.xlsx")
xeno_dir <- "output/ST_paper/"   # xenograft deliverables

# a source file lives in out_dir, or (xenograft) in xeno_dir
resolve <- function(fname) {
  for (p in c(file.path(out_dir, fname), file.path(xeno_dir, fname))) {
    if (file.exists(p)) return(p)
  }
  file.path(out_dir, fname) # default path for the skip message
}

# sheet display name -> source file, row unit, columns to wrap. Order is sheet order.
# The names here are the ones the tabs carry; Excel caps a name at 31 characters.
SPEC <- list(
  list(sheet = "Table S1. Components", file = "TableS_components.xlsx",
       unit = "component", wrap = c("top_pos", "GO_pos", "top_neg", "GO_neg")),
  list(sheet = "Table S2. GO", file = "TableS_go_all.xlsx",
       unit = "term", wrap = c("Term")),
  list(sheet = "Table S3. Comparison", file = "TableS_comparison_full.csv",
       unit = "run of a configuration", wrap = character(0)),
  list(sheet = "Table S4. Xenograft", file = "TableS_xenograft.xlsx",
       unit = "tail of a component", wrap = c("top_genes", "GO_human", "GO_mouse")))

read_any <- function(path) {
  if (grepl("\\.xlsx$", path)) openxlsx::read.xlsx(path, sheet = 1)
  else read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

HEADER <- createStyle(textDecoration = "bold", fgFill = "#DDEBF7", halign = "left",
                      valign = "center", border = "bottom", borderColour = "#8EA9DB")
BAND <- createStyle(fgFill = "#F3F7FC")
WRAP <- createStyle(wrapText = TRUE, valign = "top")
TITLE <- createStyle(textDecoration = "bold", valign = "top", wrapText = TRUE)
BOLD <- createStyle(textDecoration = "bold")
MIN_W <- 6 # narrow columns stay readable: a header of one letter still needs a filter arrow
CAP_W <- 35 # the ceiling Excel's autofit lacks

# One data frame as a styled sheet. A title, when given, sits above the table rather than in
# it: headers, filter and freeze move down two rows, so the sentence cannot be sorted or
# filtered away with the data. bold_rows are row numbers within df, not within the sheet.
style_sheet <- function(wb, sheet, df, wrap, title = NULL, bold_rows = integer(0)) {
  addWorksheet(wb, sheet)
  # A wrapped column keeps a fixed 55 and lays its text over several lines; the rest are sized
  # to their longest value, which here never exceeds 33 characters, so CAP_W guards only
  # against a long value appearing in an unwrapped column later.
  fit <- function(x) {
    n <- suppressWarnings(max(nchar(as.character(x)), na.rm = TRUE))
    if (!is.finite(n)) n <- 0L
    return(n)
  }
  widths <- ifelse(names(df) %in% wrap, 55,
                   pmin(pmax(pmax(nchar(names(df)), vapply(df, fit, numeric(1))) + 2, MIN_W), CAP_W))
  setColWidths(wb, sheet, cols = seq_len(ncol(df)), widths = widths)

  off <- if (is.null(title)) 0L else 2L # rows the table is pushed down by
  if (!is.null(title)) {
    writeData(wb, sheet, title, startRow = 1, startCol = 1)
    mergeCells(wb, sheet, cols = seq_len(ncol(df)), rows = 1)
    addStyle(wb, sheet, TITLE, rows = 1, cols = seq_len(ncol(df)), gridExpand = TRUE)
    # a merged cell does not grow to fit its text, so the height is set by hand
    setRowHeights(wb, sheet, rows = 1,
                  heights = 15 * max(1, ceiling(nchar(title) / (0.95 * sum(widths)))) + 6)
  }
  writeData(wb, sheet, df, startRow = 1 + off, headerStyle = HEADER)
  freezePane(wb, sheet, firstActiveRow = 2 + off)
  addFilter(wb, sheet, rows = 1 + off, cols = seq_len(ncol(df)))
  if (nrow(df)) {
    if ("Slide" %in% names(df)) { # band whole rows by Slide block
      blk <- cumsum(c(TRUE, df$Slide[-1] != df$Slide[-nrow(df)]))
      shaded <- which(blk %% 2 == 0)
      if (length(shaded))
        addStyle(wb, sheet, BAND, rows = shaded + 1 + off, cols = seq_len(ncol(df)),
                 gridExpand = TRUE, stack = TRUE)
    }
    wc <- which(names(df) %in% wrap)
    if (length(wc))
      addStyle(wb, sheet, WRAP, rows = 1 + off + seq_len(nrow(df)), cols = wc,
               gridExpand = TRUE, stack = TRUE)
    if (length(bold_rows))
      addStyle(wb, sheet, BOLD, rows = bold_rows + 1 + off, cols = seq_len(ncol(df)),
               gridExpand = TRUE, stack = TRUE)
  }
}

wb <- createWorkbook()
report <- character(0)

# ---- Contents sheet (first): expand data/supplement_contents.csv ----
if (file.exists(contents_csv)) {
  cc <- read.csv(contents_csv, stringsAsFactors = FALSE, check.names = FALSE)  # sheet, column, definition
  # the contents file and SPEC may spell a sheet name differently, so they are matched on the
  # S-number they share; an exact match on the string would fail silently.
  snum <- function(x) sub("^.*?S\\s*([0-9]+).*$", "\\1", x)
  # a row naming no sheet is the legend's own opening, a sentence about the workbook. It is
  # held out of the matching, raises no warning, and is placed first whatever its position.
  head_row <- is.na(cc$sheet) | trimws(cc$sheet) == ""
  ord <- c(which(head_row), which(!head_row))
  cc <- cc[ord, , drop = FALSE]
  head_row <- head_row[ord]
  units <- setNames(vapply(SPEC, `[[`, character(1), "unit"),
                    snum(vapply(SPEC, `[[`, character(1), "sheet")))
  cc$unit <- unname(units[snum(cc$sheet)])
  cc$unit[is.na(cc$unit)] <- ""
  cc$unit[head_row] <- ""
  if (all(cc$unit[!head_row] == ""))
    warning("Contents: no sheet name matched SPEC, the Unit column will be empty", call. = FALSE)
  # the unit belongs to the table, so it is written once and left blank on the rows below
  cc$unit[!is.na(cc$column) & cc$column != ""] <- ""
  # the name shown is SPEC's, so renaming a tab cannot leave the legend pointing at a name no
  # tab carries; the file's own spelling serves only to say which sheet a row belongs to
  labels <- setNames(vapply(SPEC, `[[`, character(1), "sheet"),
                     snum(vapply(SPEC, `[[`, character(1), "sheet")))
  shown <- unname(labels[snum(cc$sheet)])
  if (any(is.na(shown) & !head_row))
    warning("Contents: rows naming a sheet SPEC does not have, kept as written: ",
            paste(unique(cc$sheet[is.na(shown) & !head_row]), collapse = ", "), call. = FALSE)
  cc$sheet <- ifelse(is.na(shown), cc$sheet, shown)
  cc$sheet[head_row] <- ""
  # the opening sentence becomes the title and leaves the table
  title <- if (any(head_row)) paste(cc$definition[head_row], collapse = " ") else NULL
  cc <- cc[!head_row, , drop = FALSE]
  contents <- cc[, c("sheet", "unit", "column", "definition")]
  names(contents) <- c("Sheet", "Unit", "Column", "Definition")
  # a row that introduces a table names no column; those are bold, so the tables stand apart
  # from the column definitions under them
  tbl_rows <- which(is.na(contents$Column) | trimws(contents$Column) == "")
  style_sheet(wb, "Contents", contents, wrap = c("Definition"),
              title = title, bold_rows = tbl_rows)
  report <- c(report, sprintf("Contents: title + %d rows from %s (%d table headings in bold)",
                              nrow(contents), contents_csv, length(tbl_rows)))
} else {
  report <- c(report, sprintf("Contents sheet SKIPPED — %s absent", contents_csv))
}

# ---- data sheets ----
for (s in SPEC) {
  path <- resolve(s$file)
  if (!file.exists(path)) {
    report <- c(report, sprintf("%s SKIPPED — %s absent", s$sheet, s$file))
    next
  }
  df <- read_any(path)
  style_sheet(wb, s$sheet, df, s$wrap)
  report <- c(report, sprintf("%s: %d rows x %d cols  (%s)", s$sheet, nrow(df), ncol(df), s$file))
}

saveWorkbook(wb, wb_path, overwrite = TRUE)
cat("Wrote", normalizePath(wb_path, mustWork = FALSE), "\n")
for (r in report) cat("  -", r, "\n")
