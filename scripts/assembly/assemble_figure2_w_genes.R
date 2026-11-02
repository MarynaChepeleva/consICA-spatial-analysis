#
# Figure 2 of the manuscript: sixteen panels, four rows of four. For each of the two Visium
# slides a histology panel and three marker-expression maps, then four component weight maps.
#   MEND40   a H&E | b KRT10 | c MKI67 | d CD74 || e IC11 | f IC5 | g IC14 | h IC9
#   NCBI776  i H&E | j SCGB2A2 | k CPB1 | l CDH2 || m IC18 | n IC15 | o IC14 | p IC10
# The supplementary grids S1 and S2 are built by assemble_figure_S1_2_grids.R.
#
# Expression is rebuilt from the HEST h5ad, there being no saved Seurat object, with the QC
# the demonstration used: every spot kept, control-matching and all-zero genes dropped. The
# value plotted is ln(CP10k + 1), not SCTransform scale.data: those are Pearson residuals,
# capped, and a marker drawn from them shows a flat plateau over the very region the map
# exists to display. One pass per slide, then cached.
#
# Geometry comes from figobj, the same object the component maps were drawn from, so frame
# and spots cannot drift from them. The component panels are those saved ggplots reused
# unchanged; only the colour scale is replaced, capped as the caption states.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(zellkonverter)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})
options(future.globals.maxSize = 4 * 1024^3)

OUT_DIR <- "output/demo"
WSI_DIR <- "data/wsi_export"
HEST_DIR <- Sys.getenv("HEST_ROOT", "data/HEST")
CACHE_DIR <- "output/heap/demo"
FIG_DIR <- "output/demo" # the figures land beside the other demo output
PREVIEW <- FALSE # optional sheets with the wider gene and component lists; off by default
for (d in c(OUT_DIR, CACHE_DIR, FIG_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

CTRL <- "^(BLANK|NegControlProbe|NegControlCodeword|antisense|DeprecatedCodeword|UnassignedCodeword|Intergenic)"
SCALE_NAME <- "ln(CP10k+1)" # printed on every expression legend; the actual transform

# cfg$pt is the point size of the preview sheets only: a preview panel is larger than a
# panel of the final figure, and the final sizes would leave gaps between spots here.
SLIDES <- list(
  MEND40 = list(
    h5ad = file.path(HEST_DIR, "skin/st/MEND40.h5ad"),
    figobj = file.path(OUT_DIR, "figobj_MEND40_SCT_k20.rds"),
    genes = c("KRT14", "KRT10", "MKI67", "CD74", "COL1A1"),
    comps = c("IC14", "IC9", "IC11"), pt = 0.9, ncol = 4,
    final = c("KRT10", "MKI67", "CD74"),
    # four per tissue, so that a row holds four panels and the histology is an ordinary cell
    final_comps = c("IC11", "IC5", "IC14", "IC9"),
    he = "MEND40_o24.tif", pixel_um = 0.1719394671378331, ds_he = 24, bar_um = 500,
    trim_top = 0, # fraction of the blank slide above the tissue to drop
    # pitch_frac is the median nearest-neighbour distance as a fraction of the bounding box
    # width, here 581 px over 32862, so the printed pitch follows whatever width the map is
    # given. fpt is the dot size that suited that pitch at 28 mm and is rescaled from it.
    pitch_frac = 581 / 32862, aspect = 0.866, fpt = 0.6,
    crop_q = 0.002,                    # drops the units stranded above the tissue
    he_mar = c(0.00, 0.04, 0.04, -0.01)),  # left, right, top, bottom: cut in on the two
                                            # sides where the glass edge and its dots intrude
  NCBI776 = list(
    h5ad = file.path(HEST_DIR, "breast/st/NCBI776.h5ad"),
    figobj = file.path(OUT_DIR, "figobj_NCBI776_SCT_k20.rds"),
    genes = c("SCGB2A2", "CPB1", "CDH2", "TACSTD2", "CD3D", "POSTN", "CD68"),
    # comps is drawn only when PREVIEW is on; the figure itself uses final_comps
    comps = c("IC18", "IC15", "IC10", "IC16"), pt = 0.75, ncol = 4,
    final = c("SCGB2A2", "CPB1", "CDH2"),
    final_comps = c("IC18", "IC15", "IC14", "IC10"),
    he = "NCBI776_o24.tif", pixel_um = 0.5475531951223623, ds_he = 24, bar_um = 500,
    pitch_frac = 182 / 11688, aspect = 0.948, fpt = 0.5,
    crop_q = 0)) # nothing stranded here, 0.1% of the frame

# mirrors the type sizes of figures_pre.R
TYPE <- list(tick = 6.5, legend = 6.5, legtitle = 7, letter = 9, title = 8)
ptt <- function(p) p / 2.845

# ---- expression: rebuilt once per slide, then cached ----
build_expression <- function(cfg, slide, genes) {
  cache <- file.path(CACHE_DIR, sprintf("markerexpr_%s_lnCP10k.rds", slide))
  if (file.exists(cache)) {
    e <- readRDS(cache)
    if (all(genes %in% rownames(e$expr))) {
      cat(sprintf("  %s: cache hit, %d genes x %d spots\n", slide,
                  nrow(e$expr), ncol(e$expr)))
      return(e)
    }
    cat(sprintf("  %s: cache lacks a requested gene, rebuilding\n", slide))
  }
  if (!file.exists(cfg$h5ad)) stop("h5ad not reachable: ", cfg$h5ad)

  cat(sprintf("  %s: reading %s\n", slide, basename(cfg$h5ad)))
  sce <- readH5AD(cfg$h5ad)
  srt <- as.Seurat(sce, counts = "X", data = NULL)
  raw <- dim(srt)

  # the demonstration ran both Visium slides with filter_spots = FALSE, so every spot is
  # kept and only genes are dropped; reproduced here rather than restated
  n_ctrl <- sum(grepl(CTRL, rownames(srt)))
  srt <- srt[!grepl(CTRL, rownames(srt)), ]
  srt <- srt[rowSums(GetAssayData(srt, assay = "originalexp", layer = "counts")) > 0, ]
  cat(sprintf("  %s: %d x %d raw -> %d x %d after QC (%d controls removed)\n",
              slide, raw[1], raw[2], nrow(srt), ncol(srt), n_ctrl))

  srt <- NormalizeData(srt, assay = "originalexp", normalization.method = "LogNormalize",
                       scale.factor = 1e4, verbose = FALSE)
  full <- GetAssayData(srt, assay = "originalexp", layer = "data")

  present <- intersect(genes, rownames(full))
  e <- list(expr = as.matrix(full[present, , drop = FALSE]),
            absent = setdiff(genes, rownames(full)),
            n_genes = nrow(full), n_spots = ncol(full),
            transform = SCALE_NAME, source = cfg$h5ad, built = Sys.time())
  rm(sce, srt, full)
  gc(verbose = FALSE)
  saveRDS(e, cache)
  cat(sprintf("  %s: cached %s\n", slide, basename(cache)))
  e
}

# ---- geometry, taken from the component maps rather than rebuilt ----
# The frame may be cropped: above the tissue of MEND40 a handful of units sit far from the
# rest, and the bounding box stretches to reach them, which leaves an empty band at the top
# of every panel of that slide. crop_q drops the topmost units by quantile. The histology is
# cropped to the same box further down, so the two keep showing one field.
spot_frame <- function(fig, crop_q = 0) {
  d <- fig$maps$data
  one <- d[d$IC == d$IC[1], c("x", "y")] # one component's worth: the spots themselves
  one$barcode <- rownames(d)[d$IC == d$IC[1]]
  stopifnot(!anyDuplicated(one$barcode))
  if (crop_q > 0) {
    cut <- stats::quantile(one$y, crop_q, names = FALSE) # y grows downwards on the image
    one <- one[one$y >= cut, ]
  }
  one
}

# ---- panels ----
gene_panel <- function(sf, values, gene, pt) {
  df <- data.frame(sf, value = unname(values[sf$barcode]))
  ggplot(df, aes(x, y, colour = value)) +
    geom_point(size = pt, stroke = 0) +
    # magma reversed: cream at zero through orange to near-black at the maximum. Zero reads
    # as empty, the way white does on the weight maps, which an ordinary viridis ramp does
    # not: there zero is near-black and a sparsely detected gene fills the frame with ink.
    scale_colour_viridis_c(name = SCALE_NAME, option = "magma", direction = -1,
                           breaks = scales::breaks_pretty(n = 3)) +
    coord_fixed() + scale_y_reverse() + theme_void() +
    ggtitle(gene) +
    theme(plot.title = element_text(size = TYPE$title, face = "italic", hjust = 0.5),
          legend.position = "right", legend.justification = "top",
          legend.key.width = unit(0.25, "lines"),
          legend.key.height = unit(1.1, "lines"),
          legend.title = element_text(size = TYPE$legtitle),
          legend.text = element_text(size = TYPE$legend),
          legend.box.spacing = unit(2, "pt"),
          plot.margin = margin(4, 3, 3, 4))
}

# the component map as it was saved; only the point size is matched to the preview
comp_panel <- function(fig, ic, pt) {
  p <- fig[[ic]]
  p$layers[[1]]$aes_params$size <- pt
  p$layers[[1]]$aes_params$stroke <- 0
  p + ggtitle(paste(ic, "weights")) +
    theme(plot.title = element_text(size = TYPE$title, hjust = 0.5),
          legend.position = "right", legend.justification = "top",
          legend.key.width = unit(0.25, "lines"),
          legend.key.height = unit(1.1, "lines"),
          legend.title = element_text(size = TYPE$legtitle),
          legend.text = element_text(size = TYPE$legend),
          legend.box.spacing = unit(2, "pt"),
          plot.margin = margin(4, 3, 3, 4))
}

# ---- run ----
report <- list()
for (slide in names(SLIDES)) {
  cfg <- SLIDES[[slide]]
  cat(sprintf("\n=== %s ===\n", slide))

  e <- build_expression(cfg, slide, cfg$genes)
  fig <- readRDS(cfg$figobj)
  sf <- spot_frame(fig)

  # barcodes must correspond one to one; a silent intersection would quietly draw the genes
  # on a different set of spots from the components beside them
  bc_e <- colnames(e$expr)
  bc_f <- sf$barcode
  only_e <- setdiff(bc_e, bc_f)
  only_f <- setdiff(bc_f, bc_e)
  cat(sprintf("  barcodes: %d in the matrix, %d in figobj, %d shared\n",
              length(bc_e), length(bc_f), length(intersect(bc_e, bc_f))))
  if (length(only_e) || length(only_f)) {
    cat("  stopping: the two spot sets differ.\n")
    cat("    only in the matrix:", length(only_e),
        if (length(only_e)) paste0("e.g. ", paste(head(only_e, 3), collapse = ", ")) else "", "\n")
    cat("    only in figobj:    ", length(only_f),
        if (length(only_f)) paste0("e.g. ", paste(head(only_f, 3), collapse = ", ")) else "", "\n")
    stop("spot sets differ on ", slide, " - resolve before drawing")
  }
  cat("  barcodes correspond one to one\n")

  # detection, reported for every requested gene and never silently patched
  for (g in cfg$genes) {
    if (g %in% e$absent) {
      report[[length(report) + 1]] <- data.frame(
        Slide = slide, gene = g, present = FALSE, detected_pct = NA_real_,
        min = NA_real_, max = NA_real_, median_nonzero = NA_real_)
      next
    }
    v <- e$expr[g, sf$barcode]
    report[[length(report) + 1]] <- data.frame(
      Slide = slide, gene = g, present = TRUE,
      detected_pct = round(100 * mean(v > 0), 1),
      min = round(min(v), 3), max = round(max(v), 3),
      median_nonzero = round(stats::median(v[v > 0]), 3))
  }

  # The candidate sheets are kept behind a switch rather than deleted, since rebuilding
  # them is the slow part.
  if (PREVIEW) {
    drawn <- setdiff(cfg$genes, e$absent)
    panels <- c(lapply(drawn, function(g) gene_panel(sf, e$expr[g, ], g, cfg$pt)),
                lapply(cfg$comps, function(ic) comp_panel(fig, ic, cfg$pt)))
    sheet <- wrap_plots(panels, ncol = cfg$ncol) +
      plot_annotation(title = sprintf("%s - marker expression (%s) and component weights",
                                      slide, SCALE_NAME),
                      theme = theme(plot.title = element_text(size = 10)))
    rows <- ceiling(length(panels) / cfg$ncol)
    f <- file.path(FIG_DIR, sprintf("preview_markers_%s.pdf", slide))
    ggsave(f, sheet, width = 24, height = 7 * rows, units = "cm", device = cairo_pdf,
           limitsize = FALSE)
    cat(sprintf("  wrote %s  (%d panels, %d rows)\n", basename(f), length(panels), rows))
  }
}

rep <- do.call(rbind, report)
write.csv(rep, file.path(OUT_DIR, "marker_maps__detection.csv"), row.names = FALSE)
cat("\n=== detection, per gene ===\n")
print(rep, row.names = FALSE)
cat("\nexpression source:", SCALE_NAME, "rebuilt from the HEST h5ad; not cica$X\n")

# ---- final layout ----

source("scripts/assembly/figures_pre.R")   # pt, TYPE, letter_theme, raster_panel, map_panel, assemble_figure

# Between a third and two thirds of the units read exactly zero for these genes, so the
# bottom of the scale is "not detected" rather than "low" and is given a very pale yellow:
# white would lose the yellow end of the ramp, and magma's own lightest colour is too
# saturated to read as empty tissue. The rest of the ramp is magma reversed.
EXPR_RAMP <- c("#fff9f5", viridisLite::magma(32, direction = -1))

# The bar size is set on the guide, not in the panel theme: through the theme it is laid out
# against the whole plot and overflows a panel whose drawing area has been pinned, running
# across the panel beside it. The bar is centred under its map, as the title is centred over it.
BAR_IC <- theme(legend.key.width = unit(16, "mm"), legend.key.height = unit(1.8, "mm"),
                   legend.text = element_text(size = 6),
                   legend.title = element_text(size = TYPE$legend, hjust = 0, vjust = 1),
                   legend.title.position = "left",
                   legend.ticks = element_line(colour = "white", linewidth = 0.2),
                   legend.ticks.length = unit(0.05, "cm"))
BAR_GENE <- BAR_IC + theme(legend.title = element_text(size = TYPE$legend, hjust = 0,
                                                        vjust = 1, face = "italic"))

# panel margins: tight on the sides so neighbouring maps sit close, and just enough on top
# for the letter. These are what separates the maps, so they are the knob for "closer together".
PANEL_M <- margin(t = 9, r = 1, b = 1, l = 4)

DRAFT <- FALSE   # TRUE renders png in seconds; FALSE writes the vector pdf
# Left to itself each panel keeps whatever width is left after its own legend, and the legends
# differ: a gene bar reads 0 2 4 6 and a weight bar -0.1 0.0 0.1 0.2, which gave maps of 37 mm
# in one row and 28 mm in the next. fix_size() pins the drawing area and the legend column to
# exact millimetres on the assembled gtable, so every map is identical by construction.
HE_MAR <- 0.05                       # crop the histology this much wider than the spot frame
MAP_MM <- 34                         # width of the drawing area in every panel
WIDTH_CM <- 4 * (MAP_MM + 4.6) / 10  # follows MAP_MM: four columns plus the side margins
# The panel letter is bold 9 pt against a 7 pt title, so aligned at the top their baselines
# differ and the letter reads high. TAG_DROP lowers it onto the title's line.
CLIP_Q <- 0.01                       # the capped panel ends at this percentile, mirrored
CAP_SLIDE <- "NCBI776"               # the one panel whose scale is capped; the others
CAP_IC <- "IC14"                     # keep the full range of their component
TAG_DROP <- 0                        # points the panel letter is pushed down by
tag_fix <- theme(plot.tag = element_text(size = TYPE$letter, face = "bold", hjust = 0,
                                          vjust = 1, margin = margin(t = TAG_DROP)))
HE_SCALE <- 1.00                       # 1 = histology at full panel size; below 1 it is padded
                                       # with white, but the scale bar is placed as a fraction
                                       # of image width and would move into that padding
LEGH_MM <- 6                           # height reserved for the bar under the map
TITLE_MM <- 4                          # height of the line holding the name
# The gap between two rows is the bottom margin of one panel plus the top margin of the next.
# The top one cannot go below the height of the panel letter, which is 9 pt and anchored to
# the top of the plot, so 10 pt is the floor; the bottom one carries nothing and goes to 1.
PAN_TOP <- 10 # points above the panel, the band the letter sits in
PAN_BOT <- 1 # points below it
PAN_M2 <- margin(t = PAN_TOP, r = 3, b = PAN_BOT, l = 10)
TOP_MM <- 4.0 # top margin, the band the panel letter sits in
MARG_MM <- (PAN_TOP + PAN_BOT) / 2.845   # follows the panel margins, so trimming them
                                         # actually shortens the row instead of adding white

find_obj <- function(rds) {
  for (d in c("output/demo/fig_obj/", "output/demo/"))
    if (file.exists(paste0(d, rds))) return(paste0(d, rds))
  stop("figobj not found: ", rds)
}

# the histology panel: crop the downsampled WSI to the spot frame, then hand it to raster_panel
he_panel <- function(letter, cfg, rds, sf = NULL) {
  path <- file.path(WSI_DIR, cfg$he)
  um_per_px <- cfg$pixel_um * cfg$ds_he
  if (!file.exists(path)) {
    message(sprintf("  %s: H&E %s absent -> grey placeholder", letter, cfg$he))
    return(ggplot() +
             annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1 / cfg$aspect,
                      fill = "grey92", colour = NA) +
             annotate("text", x = 0.5, y = 0.5 / cfg$aspect, label = "H&E (to add)",
                      size = pt(7), colour = "grey45") +
             coord_fixed(xlim = c(0, 1), ylim = c(0, 1 / cfg$aspect), expand = FALSE) +
             theme_void() + labs(tag = tolower(letter)) + letter_theme)
  }
  img <- tiff::readTIFF(path)
  H <- nrow(img)
  W <- ncol(img)
  # the frame comes from the cropped spot set when one is given, so that the stain shows the
  # same field as the maps rather than the uncropped bounding box
  dat <- if (is.null(sf)) readRDS(rds)$maps$data else sf
  xr <- range(dat$x, na.rm = TRUE)
  yr <- range(dat$y, na.rm = TRUE)
  # The crop is widened beyond the spot frame side by side, as fractions of its longer side,
  # in the order left, right, top, bottom. A positive value takes in more slide and so draws
  # the tissue smaller; a negative one cuts inside the frame, which is how the edge of the
  # glass and its row of fiducial dots are kept out on the sides where they intrude.
  # The panel keeps its size either way, so nothing in the grid moves: a crop whose
  # proportions no longer match the panel is simply centred in it with white to two sides.
  m <- if (is.null(cfg$he_mar)) rep(HE_MAR, 4) else cfg$he_mar
  ref <- max(diff(xr), diff(yr))
  c0 <- max(1L, floor((xr[1] - m[1] * ref) / cfg$ds_he))
  c1 <- min(W,  ceiling((xr[2] + m[2] * ref) / cfg$ds_he))
  r0 <- max(1L, floor((yr[1] - m[3] * ref) / cfg$ds_he))
  r1 <- min(H,  ceiling((yr[2] + m[4] * ref) / cfg$ds_he))
  crop <- if (length(dim(img)) == 3) img[r0:r1, c0:c1, , drop = FALSE] else img[r0:r1, c0:c1]

  # trim_top: drop part of the blank slide above the tissue. The first tissue row is found by
  # colour saturation, since H&E stain is coloured and bare slide is grey, and the fraction
  # given is how much of the empty band above it to remove. The crop then covers slightly less
  # than the spot bounding box, so the histology shows a little less field than the maps do;
  # the spots in that band are still drawn on them.
  if (!is.null(cfg$trim_top) && cfg$trim_top > 0 && length(dim(crop)) == 3) {
    sat <- apply(crop[, , 1:3, drop = FALSE], 1,
                  function(r) mean(apply(matrix(r, ncol = 3), 1, function(p) max(p) - min(p))))
    first <- which(sat > 0.35 * max(sat))[1]
    if (!is.na(first) && first > 1) {
      cut <- round(cfg$trim_top * (first - 1))
      message(sprintf("  %s: tissue starts at row %d of %d, trimming %d rows from the top",
                      letter, first, nrow(crop), cut))
      if (cut >= 1) crop <- crop[(cut + 1):nrow(crop), , , drop = FALSE]
    }
  }

  # HE_SCALE shrinks the picture, not the panel: white rows and columns are added around the
  # crop itself, so the gtable keeps its full size and the panel letter stays on the line it
  # shares with the letters beside it. Shrinking the panel instead centres it in the cell and
  # drags the letter down with it.
  if (HE_SCALE < 1 && length(dim(crop)) == 3) {
    padr <- round(nrow(crop) * (1 / HE_SCALE - 1) / 2)
    padc <- round(ncol(crop) * (1 / HE_SCALE - 1) / 2)
    big <- array(1, dim = c(nrow(crop) + 2 * padr, ncol(crop) + 2 * padc, dim(crop)[3]))
    big[padr + seq_len(nrow(crop)), padc + seq_len(ncol(crop)), ] <- crop
    crop <- big
    message(sprintf("  %s: padded to %d x %d px so the picture is %.0f%% of the panel",
                    letter, ncol(crop), nrow(crop), 100 * HE_SCALE))
  }
  message(sprintf("  %s: crop %d x %d px, %g um bar", letter, ncol(crop), nrow(crop), cfg$bar_um))
  raster_panel(crop, ds = um_per_px / PX_UM, letter, bar_um = cfg$bar_um,
               bar_txt = 5, bar_x = 0.08, bar_lwd = 0.5)
}

fix_size <- function(p, w, h, legw) {
  gt <- ggplotGrob(p)
  pan <- gt$layout$name == "panel"
  gt$widths[unique(gt$layout$l[pan])] <- unit(w, "mm")
  gt$heights[unique(gt$layout$t[pan])] <- unit(h, "mm")
  gb <- grep("guide-box-bottom", gt$layout$name)
  if (length(gb)) gt$heights[gt$layout$t[gb]] <- unit(legw, "mm")
  patchwork::wrap_elements(full = gt)
}

cat("\n=== final layout: the published grid, four rows ===\n")

gene_map <- function(sf, values, gene, letter, pt_size) {
  df <- data.frame(sf, value = unname(values[sf$barcode]))
  ggplot(df, aes(x, y, colour = value)) +
    geom_point(size = pt_size, stroke = 0) +
    scale_colour_gradientn(name = NULL, colours = EXPR_RAMP,
                           breaks = scales::breaks_pretty(n = 3),
                           guide = guide_colourbar(theme = BAR_GENE)) +
    coord_fixed() + scale_y_reverse() +
    ggtitle(gene) + labs(tag = letter) +
    theme_void(base_size = 8) +
    theme(plot.title = element_text(size = TYPE$legtitle, face = "italic", hjust = 0.5,
                                    margin = margin(t = 3, b = 1))) +
    theme(legend.position = "bottom", legend.justification = "center",
          legend.box.spacing = unit(1, "pt"),
          legend.margin = margin(0, 0, 0, 0)) +
    letter_theme + tag_fix + theme(plot.margin = PAN_M2)
}

rows <- list()
letters_used <- character(0)
i <- 0
row_mm <- numeric(0)
for (slide in names(SLIDES)) {
  cfg <- SLIDES[[slide]]
  rds <- find_obj(basename(cfg$figobj))
  fig <- readRDS(rds)
  sf <- spot_frame(fig, crop_q = cfg$crop_q)
  # the proportions follow the cropped frame, not a constant: dropping the band changes them
  asp <- diff(range(sf$x)) / diff(range(sf$y))
  if (!is.null(cfg$crop_q) && cfg$crop_q > 0) {
    n0 <- nrow(fig$maps$data) / 20
    cat(sprintf("  %-8s frame cropped: %d of %.0f units dropped, aspect %.3f -> %.3f\n",
                slide, round(n0) - nrow(sf), n0, cfg$aspect, asp))
  }
  e <- readRDS(file.path(CACHE_DIR, sprintf("markerexpr_%s_lnCP10k.rds", slide)))

  ltr <- letters[(i + 1):(i + 8)]
  i <- i + 8
  letters_used <- c(letters_used,
                    sprintf("%s: %s=H&E %s=%s | %s=%s", slide, ltr[1],
                            paste(ltr[2:4], collapse = ""), paste(cfg$final, collapse = ","),
                            paste(ltr[5:8], collapse = ""),
                            paste(cfg$final_comps, collapse = ",")))

  # a blank title line, so the letter stands above the image as it stands above a map title
  # and the top of the histology lines up with the top of the maps beside it
  he <- he_panel(ltr[1], cfg, rds, sf) +
    ggtitle(" ") +
    theme(plot.title = element_text(size = TYPE$legtitle, margin = margin(t = 3, b = 1))) +
    theme(plot.margin = PAN_M2) + tag_fix
  # The dot has to keep its ratio to the printed spot pitch, or a change of MAP_MM silently
  # turns a surface into a scatter. fpt in the config is the size that looked right at 28 mm,
  # so the ratio is read back from it and reapplied at whatever width is in force.
  dot <- cfg$fpt / (28 * cfg$pitch_frac) * MAP_MM * cfg$pitch_frac
  gs <- Map(function(g, l) gene_map(sf, e$expr[g, ], g, l, dot),
            cfg$final, ltr[2:4])
  cs <- Map(function(ic, l) {
    p <- fig[[ic]]
    # the same crop as the frame, or the panel would still reach to the dropped units
    p$data <- p$data[rownames(p$data) %in% sf$barcode, , drop = FALSE]
    p$layers[[1]]$aes_params$stroke <- 0
    # Only one panel is capped, and only because its scale is unreadable without it: the
    # maximum of NCBI776 IC14 is 6.4 times its 99th percentile, 0.999 against 0.157, so 84%
    # of the palette serves 50 units of 4992 and the component reads as empty although 37.9%
    # of its units exceed 0.1. The other seven run from 1.4 to 2.5 and read as they are. What
    # is not readable is corrected; what is merely not optimal is left alone.
    #
    # The limits are c(min, q) with q the 99th percentile: the bar covers the colours the
    # map actually carries, which is the rule the other seven follow. Symmetric limits
    # c(-q, q) were tried first and they colour the map identically, because a diverging
    # scale lays its colour out over 2 * max(|low|, |high|) whatever the limits say, so
    # the half-extent is q either way. What they change is the bar: it then reaches the
    # blue end of the palette, printing a navy no unit of this component attains, -0.157
    # against a minimum of -0.116. Only the red end departs from the other panels, and it
    # departs because that is the cap, which the caption states.
    v <- p$data$activity
    qq <- stats::quantile(v, 1 - CLIP_Q, names = FALSE)
    cap <- identical(slide, CAP_SLIDE) && identical(ic, CAP_IC)
    lims <- if (cap) c(min(v), qq) else range(v)
    if (cap)
      cat(sprintf("  %-8s %-5s scale %.3f .. %.3f (range was %.3f .. %.3f); %d of %d units held at the red end, %d at the blue\n",
                  slide, ic, lims[1], lims[2], min(v), max(v),
                  sum(v > qq), length(v), sum(v < -qq)))
    map_panel(p, l, NULL, pt_size = dot, leg_pos = "bottom") +
      # a fresh scale rather than an edit of the saved one: in this ggplot the fields of an
      # assembled scale no longer take an assignment, and the limits were silently ignored.
      # Colours are the demonstration's own, so only the ends of the range change.
      scale_colour_gradient2(low = "#3c5488", mid = "gray92", high = "#9e2b32", midpoint = 0,
                             limits = lims, oob = scales::squish, name = NULL,
                             # the same round ticks as the gene bars, and as the two-row
                             # figure carried before them. The default rule puts five to
                             # seven labels on a 16 mm bar and they overlap; a horizontal
                             # bar has room for three or four and no more.
                             breaks = scales::breaks_pretty(n = 3),
                             guide = guide_colourbar(theme = BAR_IC)) +
      ggtitle(ic) +
      theme(plot.title = element_text(size = TYPE$legtitle, hjust = 0.5,
                                      margin = margin(t = 3, b = 1))) +
      theme(legend.justification = "center",
            legend.box.spacing = unit(1, "pt"),
            legend.margin = margin(0, 0, 0, 0),
            plot.margin = PAN_M2) + tag_fix
  }, cfg$final_comps, ltr[5:8])
  # every panel of this tissue gets exactly the same drawing area, histology included: its
  # crop is the spot bounding box plus 3%, so it carries the same proportions as the maps
  map_h <- MAP_MM / asp
  row <- map_h + TITLE_MM + LEGH_MM + MARG_MM
  row_mm <- c(row_mm, row, row)
  cat(sprintf("  %-8s panel %.1f x %.1f mm, row height %.1f mm, both rows\n",
              slide, MAP_MM, map_h, row))
  # the maps all take the same forced drawing area; the histology takes a smaller one, so it
  # sits in the same cell with white around it and nothing else on the row moves
  sized <- function(ps) lapply(ps, function(p) fix_size(p, MAP_MM, map_h, LEGH_MM))
  he_s <- fix_size(he, MAP_MM, map_h, LEGH_MM)
  rows <- c(rows, list(c(list(he_s), sized(gs))), list(sized(cs)))
  cat(sprintf("  %-8s dot %.2f mm at a %.0f mm map (%.2f was the size at 28 mm)\n",
              slide, dot, MAP_MM, cfg$fpt))
}
for (s in letters_used) cat("  ", s, "\n")
out <- file.path(FIG_DIR, if (DRAFT) "Figure2_draft.png" else "Figure2_final_markers.pdf")
# assemble_figure() is hardwired to cairo_pdf, so the composition it would build is made here
# and the device is chosen, which is what lets a draft render as png in seconds
row_plots <- lapply(rows, function(ps) wrap_plots(ps, nrow = 1, widths = rep(1, 4)))
fig_all <- wrap_plots(row_plots, ncol = 1, heights = row_mm)
total_mm <- sum(row_mm)
cat(sprintf("  figure %.0f x %.0f mm | the journal page leaves about 200 mm\n", WIDTH_CM * 10, total_mm))
ggsave(out, fig_all, width = WIDTH_CM, height = total_mm / 10, units = "cm",
       device = if (DRAFT) ragg::agg_png else grDevices::cairo_pdf,
       dpi = if (DRAFT) 200 else 300, limitsize = FALSE)
cat("wrote", normalizePath(out, mustWork = FALSE), "\n")
