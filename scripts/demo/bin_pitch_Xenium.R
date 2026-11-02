# Xenium pooling step (HEST): the grid step measured directly from the coordinates.
# Settles 55 um against 100 um for the Xenium unit. HEST leaves the step unrecorded,
# spot_diameter and inter_spot_dist being NaN on Xenium, so it is measured from the grid.
# All units are used, with no QC filter: filtering punches holes in the grid and inflates
# the nearest-neighbour distance.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(zellkonverter)
})

HEST_DIR <- Sys.getenv("HEST_ROOT", "data/HEST")

# per-slide constants from data/metadata/<slide>.json (full precision)
SLIDES <- list(
  TENX117 = list(file = file.path(HEST_DIR, "skin/st/TENX117.h5ad"),
                 pixel_um = 0.1368527629551209, fw = 22209, fh = 43287, n_expect = 1830),
  TENX115 = list(file = file.path(HEST_DIR, "skin/st/TENX115.h5ad"),
                 pixel_um = 0.2736897890251921, fw = 20342, fh = 22673, n_expect = 3886))

# modal positive consecutive difference of sorted unique values = one grid step.
# (gaps between occupied grid lines are multiples of the step; the mode is the step.)
modal_step <- function(v) {
  u <- sort(unique(v))
  d <- diff(u); d <- d[d > 1e-6]          # positive diffs only
  dr <- round(d, 1)                         # collapse float noise
  tab <- sort(table(dr), decreasing = TRUE)
  list(mode = as.numeric(names(tab)[1]), mode_n = as.integer(tab[1]),
       min = min(d), n_unique = length(u), top = head(tab, 5))
}

measured <- list()

for (sl in names(SLIDES)) {
  s <- SLIDES[[sl]]
  sce <- readH5AD(s$file)
  cat(sprintf("\n=====================  %s  =====================\n", sl))
  cat(sprintf("units (UNFILTERED): %d  | metadata spots_under_tissue = %d\n", ncol(sce), s$n_expect))

  # coordinate space: raw spatial coords vs fullres pixel dimensions
  stopifnot("spatial" %in% reducedDimNames(sce))
  xy <- as.data.frame(reducedDim(sce, "spatial"))[, 1:2]
  colnames(xy) <- c("x", "y")
  cat(sprintf("max(x) = %.1f  vs fullres_px_width  = %d\n", max(xy$x), s$fw))
  cat(sprintf("max(y) = %.1f  vs fullres_px_height = %d\n", max(xy$y), s$fh))
  cat("  -> if max(x/y) match fullres px, coordinates are in full-resolution pixels.\n")

  # modal grid step per axis, in pixels
  sx <- modal_step(xy$x); sy <- modal_step(xy$y)
  cat(sprintf("modal step  x = %.1f px (mode n=%d, min=%.1f, unique=%d)\n", sx$mode, sx$mode_n, sx$min, sx$n_unique))
  cat(sprintf("modal step  y = %.1f px (mode n=%d, min=%.1f, unique=%d)\n", sy$mode, sy$mode_n, sy$min, sy$n_unique))
  cat("top-5 x diffs (px : count):\n"); print(sx$top)
  cat("top-5 y diffs (px : count):\n"); print(sy$top)

  # convert to um, then the verdict
  ux <- sx$mode * s$pixel_um; uy <- sy$mode * s$pixel_um
  near <- function(u) c(55, 100)[which.min(abs(u - c(55, 100)))]
  cat(sprintf("step  x = %.2f um  (nearest of {55,100} = %d)\n", ux, near(ux)))
  cat(sprintf("step  y = %.2f um  (nearest of {55,100} = %d)\n", uy, near(uy)))
  cat(sprintf("axes agree (px): %s  | verdict: %d um pooling\n",
              isTRUE(all.equal(sx$mode, sy$mode, tolerance = 0.02)), near((ux + uy) / 2)))

  # look for a pooling field in the object / uns, and any HEST version tag
  md <- metadata(sce)
  cat("metadata()/uns names:\n"); print(names(md))
  hit <- grep("spot|pool|bin|size|pixel|hest|version|resolution", names(md), ignore.case = TRUE, value = TRUE)
  if (length(hit)) { cat("matching uns fields:\n"); print(md[hit]) } else cat("no spot_size/pooling/version field in uns.\n")
  cd_hit <- grep("spot|pool|bin|size", colnames(colData(sce)), ignore.case = TRUE, value = TRUE)
  if (length(cd_hit)) { cat("matching colData columns:\n"); print(cd_hit) }
  measured[[sl]] <- data.frame(slide = sl, n_units = ncol(sce),
                               step_px_x = sx$mode, step_px_y = sy$mode,
                               step_um_x = round(ux, 2), step_um_y = round(uy, 2),
                               pitch_um = round((ux + uy) / 2, 2))
  rm(sce); gc()
}

cat(sprintf("\nzellkonverter version: %s\n", as.character(packageVersion("zellkonverter"))))
cat("The HEST version is not stored in these objects; report it from the download environment.\n")

output_file <- "output/demo/bin_pitch_Xenium.csv"
if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)
write.csv(do.call(rbind, measured), output_file, row.names = FALSE)
message("written: ", output_file)
