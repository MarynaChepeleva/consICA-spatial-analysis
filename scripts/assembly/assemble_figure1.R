# Assemble manuscript Figure 1 from the functions in figures_pre.R.
#
# Panels A-D are raster H&E from data/wsi_export (TENX117). The inset rectangles come
# from QuPath and are part of the exported image; only the scale bar is drawn here.
# Panels E-G are the saved ggplot objects of output/demo/fig_obj, composed directly and
# left vector; nothing is recomputed.
#
# E and F are framed to the same full-resolution extent as the overview, so at equal
# panel width a weight map overlays the H&E pixel for pixel. Without the H&E exports
# the panels fall back to grey and the maps keep their own extent.

source("scripts/assembly/figures_pre.R")

wsi_dir <- "data/wsi_export"
fig_dir <- "output/demo"
obj_rds <- file.path(fig_dir, "fig_obj", "figobj_TENX117_SCT_k8.rds")
out_pdf <- file.path(fig_dir, "Figure1.pdf")

# inset border colours; the QuPath rectangles on the overview carry the same three
INSET_COL <- c(B = "#33a457", C = "#c8641e", D = "#3d9bd6")

read_he <- function(name, letter) {
  path <- file.path(wsi_dir, name)
  if (file.exists(path)) return(tiff::readTIFF(path, native = TRUE))
  message(sprintf("  %s: H&E %s absent -> grey placeholder", letter, name))
  matrix("grey85", nrow = 600, ncol = 600)
}

obj <- readRDS(obj_rds)

rects <- data.frame(x = c(14744, 14673, 5035), y = c(2795, 16402, 27233),
                    w = 698, h = 1099, lab = c("B", "C", "D"), stringsAsFactors = FALSE)

# overview once; its full-res extent is what E and F span to overlay A. Without the
# export there is no extent to match, and the maps keep their own.
he_ok <- file.exists(file.path(wsi_dir, "TENX117_col_ann_t30.tif"))
ov_img <- read_he("TENX117_col_ann_t30.tif", "A")
W_full <- ncol(ov_img) * DS_OV
H_full <- nrow(ov_img) * DS_OV
message(sprintf("overview %d x %d px -> full-res extent %d x %d", ncol(ov_img), nrow(ov_img), W_full, H_full))

p_A <- raster_panel(ov_img, DS_OV, "A", bar_um = 500, rects = rects)
p_B <- raster_panel(read_he("TENX117 (1, x=14744, y=2795, w=698, h=1099).tif", "B"),
                    DS_IN, "B", bar_um = 20, border = INSET_COL["B"])
p_C <- raster_panel(read_he("TENX117 (1, x=14673, y=16402, w=698, h=1099).tif", "C"),
                    DS_IN, "C", bar_um = 20, border = INSET_COL["C"])
p_D <- raster_panel(read_he("TENX117 (1, x=5035, y=27233, w=698, h=1099).tif", "D"),
                    DS_IN, "D", bar_um = 20, border = INSET_COL["D"])
p_E <- map_panel(obj$IC7, "E", "IC7",
                 xlim = if (he_ok) c(0, W_full), ylim = if (he_ok) c(0, H_full), shift = 0.16,
                 leg_pos = c(1.15, 0.5)) # E legend pushed right, off the map
p_F <- map_panel(obj$IC6, "F", "IC6",
                 xlim = if (he_ok) c(0, W_full), ylim = if (he_ok) c(0, H_full),
                 leg_pos = c(0.98, 0.5), lpad = 18) # F pushed right -> room for E's legend
p_G <- moran_panel(obj$moran, "G", lpad = 18) # G pushed right too

# Relative widths. wA is the overview and map cell. A smaller wEF pulls the legend up
# against the map and removes the whitespace coord_fixed leaves, at the cost of shrinking the map itself.
wA <- 3.7; wEF <- wA; wI <- 2.4; wG <- 5.0 # wEF = wA so E,F align under A
assemble_figure(
  rows = list(list(p_A, p_B, p_C, p_D), list(p_E, p_F, p_G)),
  widths = list(c(wA, wI, wI, wI), c(wEF, wEF, wG)),
  heights = c(1, 1),
  file = out_pdf, width_cm = 17.5, height_cm = 15)
