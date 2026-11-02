#
# Supplementary component grids S1 and S2, as single vector pdfs 17.5 cm wide at final size.
# Reuses the functions of figures_pre.R.
#   S1  Xenium TENX117 at k = 8, one grid, no panel letter
#   S2  Visium MEND40 and NCBI776 at k = 20, one grid per slide, letters a and b
# Both grids are drawn from $maps of the saved figobj, which holds per-component z-scores
# winsorised at |z| = 3. The share of values that fall beyond the limit is reported at the
# end and quoted in the figure captions.

source("scripts/assembly/figures_pre.R")

fig_dir <- "output/demo"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

out_pdf_s1 <- file.path(fig_dir, "FigureS1.pdf") # Xenium (TENX117) component grid
out_pdf_s2 <- file.path(fig_dir, "FigureS2.pdf") # Visium component grids (MEND40, NCBI776)
S2_PT <- 0.1 # dot size in the S2 (Visium) grids (mm; stroke=0); minimal
S1_PT <- 0.7 # dot size in the S1 (Xenium) grid — fewer components, larger facets

# figobj rds saved by the demos; the Xenium demo writes into fig_obj, the Visium one beside it
find_obj <- function(rds) {
  for (d in c(file.path(fig_dir, "fig_obj"), fig_dir)) {
    p <- file.path(d, rds)
    if (file.exists(p)) return(p)
  }
  stop("figobj not found: ", rds, " - run the demo to generate it")
}

# One grid from the $maps of a saved figobj. The panel letter is optional, and omitted on S1.
grid_panel <- function(rds, letter = NULL, pt_size = S2_PT, leg_pull = 0) {
  g <- readRDS(find_obj(rds))$maps
  g$layers[[1]]$aes_params$size <- pt_size # the maps carry a single geom_point layer
  g$layers[[1]]$aes_params$stroke <- 0 # drop the stroke floor so pt_size controls the dot
  message(sprintf("  grid %s: %d layer(s), layer1 = %s, size -> %g",
                  rds, length(g$layers), class(g$layers[[1]]$geom)[1], pt_size))
  g <- g + ggtitle(NULL) +
    # the legend takes a right column of its own and sits at the top of it, clear of the
    # facets; a negative left margin drags it back into the whitespace beside them
    theme(plot.margin = margin(12, 3, 3, 10),
          legend.position = "right",
          legend.justification = "top",
          legend.box.spacing = unit(2, "pt"),
          legend.margin = margin(0, 0, 0, -leg_pull))
  if (!is.null(letter)){ # panel tag only for the multi-grid figure (S2 = a, b)
    g <- g + labs(tag = tolower(letter)) +
      theme(plot.tag = element_text(size = TYPE$letter, face = "bold"), plot.tag.position = c(0, 1))
  }
  return(g)
}

# S1 — Xenium (TENX117, k=8): a single grid, no panel letter
assemble_figure(
  rows = list(list(grid_panel("figobj_TENX117_SCT_k8.rds",
                              pt_size = S1_PT, leg_pull = 40))),
  widths = list(1),
  heights = 1,
  file = out_pdf_s1, width_cm = 17.5, height_cm = 10)

# S2 — Visium (MEND40, NCBI776): one grid per row, letters a and b
assemble_figure(
  rows = list(list(grid_panel("figobj_MEND40_SCT_k20.rds", "A")),
              list(grid_panel("figobj_NCBI776_SCT_k20.rds", "B"))),
  widths = list(1, 1),
  heights = c(1, 1),
  file = out_pdf_s2, width_cm = 17.5, height_cm = 24)

# Share of values drawn at the limit, quoted in the captions of S1 and S2.
clipped <- do.call(rbind, lapply(
  c("figobj_TENX117_SCT_k8.rds", "figobj_MEND40_SCT_k20.rds", "figobj_NCBI776_SCT_k20.rds"),
  function(f) {
    d <- readRDS(find_obj(f))$maps$data
    data.frame(figobj = f, pct_clipped = round(100 * mean(abs(d$activity) > 3, na.rm = TRUE), 2))
  }))
print(clipped)
write.csv(clipped, file.path(fig_dir, "clipped_at_z3.csv"), row.names = FALSE)
