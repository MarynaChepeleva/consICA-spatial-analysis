# Run before Figures 1 and 2
# Figures are built at final size, 17.5 cm wide.
#
# Scale bars are computed from the pixel size. PX_UM is the full-resolution pixel of
# TENX117; a caller on another slide passes ds so that PX_UM * ds is its own um per pixel.

library(ggplot2)
library(patchwork)
library(tiff)
library(grid)

PX_UM <- 0.13685  # full-res um per pixel
DS_OV <- 8        # overview downsample
DS_IN <- 1        # inset downsample

pt <- function(p) p / 2.845 # pt -> ggplot text size (mm)
TYPE <- list(tick = 6.5, axis = 7.5, legend = 6.5, legtitle = 7, letter = 9, bar = 6.5)

DRAW_LETTERS <- FALSE # TRUE labels the rectangles on A instead of colour-coding them
PT_SIZE <- 0.8   # weight-map point size at final panel size

# panel letter in the top-left margin, above the panel, never on the image
letter_theme <- theme(
  plot.tag = element_text(size = TYPE$letter, face = "bold", hjust = 0, vjust = 1),
  plot.tag.position = c(0, 1),
  plot.margin = margin(t = 12, r = 3, b = 3, l = 10))

# ---- raster panel: image on its own pixel axes, optional scale bar + inset border ----
raster_panel <- function(img, ds, letter, bar_um = NA_real_, border = NA, rects = NULL,
                         bar_txt = TYPE$bar, bar_x = 0.05, bar_lwd = 0.7) {
  H <- nrow(img)
  W <- ncol(img)
  p <- ggplot() +
    annotation_custom(grid::rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"),
                                       interpolate = TRUE), xmin = 0, xmax = W, ymin = 0, ymax = H) +
    coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE) +
    theme_void() + labs(tag = tolower(letter)) + letter_theme

  if (!is.na(border))
    p <- p + annotate("rect", xmin = 0, xmax = W, ymin = 0, ymax = H,
                      fill = NA, colour = border, linewidth = 0.9)

  if (!is.na(bar_um)) { # scale bar, bottom-left, in image px
    bar_px <- bar_um / (PX_UM * ds)
    x0 <- bar_x * W
    y0 <- 0.06 * H
    p <- p +
      annotate("segment", x = x0, xend = x0 + bar_px, y = y0, yend = y0,
               linewidth = bar_lwd, colour = "black") +
      annotate("text", x = x0 + bar_px / 2, y = y0 + 0.025 * H,
               label = paste0(bar_um, " µm"), size = pt(bar_txt), vjust = 0)
    message(sprintf("  %s: %g um bar = %.1f px (panel is %d px wide)", letter, bar_um, bar_px, W))
  }

  # optional small letters to the left of each rectangle, on the overview only
  if (DRAW_LETTERS && !is.null(rects)) {
    r <- rects
    r$cy <- H - (r$y + r$h / 2) / ds
    r$lx <- r$x / ds - 0.02 * W
    p <- p + annotate("text", x = r$lx, y = r$cy, label = r$lab, colour = "black",
                      size = pt(TYPE$bar), fontface = "bold", hjust = 1)
  }
  p
}

# ---- weight-map panel: shrink points, drop axes, title the legend by component ----
map_panel <- function(p, letter, leg_title, pt_size = PT_SIZE, xlim = NULL, ylim = NULL, shift = 0,
                      leg_pos = c(1, 0.5), lpad = 0) {
  # a ggplot is an environment, so these two assignments change p itself and not a copy:
  # the object the caller passed in is modified as well. Each figobj is read fresh and
  # used once, so nothing downstream sees the edit. labs() cannot do the rename, because
  # a scale that carries a name overrides it.
  if (length(p$layers)) p$layers[[1]]$aes_params$size <- pt_size
  for (s in p$scales$scales)
    if (any(c("colour", "color") %in% s$aesthetics)) s$name <- leg_title
  p <- p + ggtitle(NULL) + labs(tag = tolower(letter))
  # shift, as a fraction of width, slides the map; clip off keeps the edge points whole
  if (!is.null(xlim))
    p <- p + coord_fixed(xlim = xlim - shift * diff(xlim), ylim = ylim, expand = FALSE, clip = "off")
  p + theme_void(base_size = 8) +
    theme(legend.position = leg_pos,                    # per-panel (E, F independent); inside the panel ->
          legend.justification = c(0, 0.5),             # panel = full width, so the map aligns under A
          legend.text = element_text(size = TYPE$legend),
          legend.title = element_text(size = TYPE$legtitle),
          legend.key.height = unit(3.2, "mm"), legend.key.width = unit(2.4, "mm"),
          legend.ticks = element_line(colour = "white", linewidth = 0.2),
          legend.ticks.length = unit(0.06, "cm")) +
    letter_theme +
    theme(plot.margin = margin(t = 12, r = 3, b = 3, l = 10 + lpad))   # lpad pushes this panel right
}

# ---- Moran panel: small type, no title. ratio is panel height over width, 0.5 giving
# a 2:1 wide-short panel; pt_size and seg_lwd shrink the lollipops; lpad pushes it right ----
moran_panel <- function(p, letter, ratio = 0.75, pt_size = 1.7, seg_lwd = 0.4, lpad = 0) {
  for (l in p$layers) {
    g <- class(l$geom)[1]
    if (g == "GeomPoint") l$aes_params$size <- pt_size
    if (g == "GeomSegment") l$aes_params$linewidth <- seg_lwd
  }
  p + ggtitle(NULL) + labs(tag = tolower(letter)) +
    theme(plot.title = element_blank(),
          axis.text = element_text(size = TYPE$tick, colour = "black"),
          axis.title = element_text(size = TYPE$axis, colour = "black"),
          aspect.ratio = ratio) +
    letter_theme +
    theme(plot.margin = margin(t = 12, r = 3, b = 3, l = 10 + lpad))
}

# ---- composer: rows is a list of panel lists, one per row ----
assemble_figure <- function(rows, widths, heights, file, width_cm, height_cm) {
  row_plots <- Map(function(ps, w) wrap_plots(ps, nrow = 1, widths = w), rows, widths)
  fig <- wrap_plots(row_plots, ncol = 1, heights = heights)
  ggsave(file, fig, width = width_cm, height = height_cm, units = "cm",
         device = grDevices::cairo_pdf, limitsize = FALSE)
  cat("wrote", normalizePath(file, mustWork = FALSE), "\n")
  invisible(fig)
}
