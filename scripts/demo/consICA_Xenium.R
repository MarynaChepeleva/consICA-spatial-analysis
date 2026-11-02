# Xenium demo: load, preprocess, apply consICA, report and spatial component maps.
# Writes the annotation tables, ics_stat, the figures and the saved plot objects
# that the assembly scripts read.
# One slide per run: set slide_id under Slide params.

library(SingleCellExperiment)
library(zellkonverter)
library(Seurat)
library(consICA)
library(ggplot2)
library(spdep)

## Slide params ----

HEST_DIR <- Sys.getenv("HEST_ROOT", "data/HEST")
SLIDES <- list(TENX115 = "skin/st/TENX115.h5ad",
               TENX117 = "skin/st/TENX117.h5ad")

slide_id <- "TENX117"   # <- set the slide to run
input_file <- file.path(HEST_DIR, SLIDES[[slide_id]])

ncomp <- 8
nc <- 1                  # consICA cores; 1 is fully reproducible, and from 2.11.2 a seed is available in parallel computations

norm_method <- "SCT"     # input to consICA: "SCT" | "LogNorm" | "log2"
filter_spots <- TRUE     # count-based per-spot QC; FALSE = keep all spots
postfix <- norm_method   # tag for output filenames

# Xenium technical controls
ctrl_pattern <- "^(BLANK|NegControlProbe|NegControlCodeword|antisense|DeprecatedCodeword|UnassignedCodeword|Intergenic)"


## Output folders and functions ----
output_folder <- "output/demo/"
heap_folder <- "output/heap/demo/"
for (d in c(output_folder, paste0(output_folder, "fig_obj/"), heap_folder))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

source("scripts/demo/functions.R")


## Data load ----
sce <- readH5AD(input_file)
srt <- as.Seurat(sce, counts = "X", data = NULL)

# spatial coordinates
cd <- as.data.frame(colData(sce))
if ("spatial" %in% reducedDimNames(sce)) {
  coord <- as.data.frame(reducedDim(sce, "spatial"))[, 1:2]
} else {
  xc <- intersect(c("pxl_col_in_fullres", "array_col", "imagecol"), colnames(cd))[1]
  yc <- intersect(c("pxl_row_in_fullres", "array_row", "imagerow"), colnames(cd))[1]
  coord <- cd[, c(xc, yc)]
}
colnames(coord) <- c("x", "y")
rownames(coord) <- colnames(sce)

cat(sprintf("Slide_id: %s | raw: %d genes x %d spots\n", slide_id, nrow(srt), ncol(srt)))

## Preprocess ----
counts <- GetAssayData(srt, assay = "originalexp", layer = "counts")

# Count-based per-spot QC
if (filter_spots) {
  min_genes <- auto_spot_threshold(counts)
  keep_spots <- Matrix::colSums(counts > 0) >= min_genes
} else {
  min_genes <- 0L
  keep_spots <- rep(TRUE, ncol(srt))
}

# verify the cut removed background, not tissue (auto can over-cut full-tissue slides)
print(
  ggplot(data.frame(coord[colnames(srt), ],
                    spot = ifelse(keep_spots, "tissue", "removed")),
         aes(x, y, colour = spot)) +
    geom_point(size = 1.5) +
    scale_colour_manual(values = c(tissue = "#2c7fb8", removed = "grey80")) +
    coord_fixed() + scale_y_reverse() + theme_void() +
    ggtitle(sprintf("%s — QC: %d kept / %d removed (min_genes=%d)",
                    slide_id, sum(keep_spots), sum(!keep_spots), min_genes))
)

srt <- srt[, keep_spots]
n_ctrl <- sum(grepl(ctrl_pattern, rownames(srt)))
srt <- srt[!grepl(ctrl_pattern, rownames(srt)), ]
srt <- srt[rowSums(GetAssayData(srt, assay = "originalexp", layer = "counts")) > 0, ]
cat(sprintf("\nQC: kept %d / %d spots (min_genes=%d) | removed %d controls; %d genes x %d spots\n",
            sum(keep_spots), length(keep_spots), min_genes, n_ctrl, nrow(srt), ncol(srt)))

# Overview
tc <- colSums(GetAssayData(srt, assay = "originalexp", layer = "counts"))
df_tc <- data.frame(coord[colnames(srt), ], total_counts = tc)
print(
  ggplot(df_tc, aes(x, y, colour = total_counts)) +
    geom_point(size = 2) + scale_colour_viridis_c(trans = "log10") +
    coord_fixed() + scale_y_reverse() + theme_void() +
    ggtitle(sprintf("%s — spots by total counts (n=%d)", slide_id, ncol(srt)))
)

# Normalization
counts <- GetAssayData(srt, assay = "originalexp", layer = "counts")   # final (post-QC) counts
if (norm_method == "SCT") {
  srt <- suppressWarnings(SCTransform(srt, assay = "originalexp",
                                      return.only.var.genes = FALSE, verbose = FALSE))
  X <- as.matrix(srt[["SCT"]]@scale.data)                             # Pearson residuals
} else if (norm_method == "LogNorm") {
  srt <- NormalizeData(srt, assay = "originalexp",
                       normalization.method = "LogNormalize", scale.factor = 1e4)
  X <- as.matrix(GetAssayData(srt, assay = "originalexp", layer = "data"))  # ln(CP10k+1)
} else if (norm_method == "log2") {
  X <- as.matrix(log2(counts + 1))                                    # no depth normalization
} else stop("Unknown norm_method: ", norm_method)
dim(X)

## Signal deconvolution with consICA ----
set.seed(2026)
cica <- consICA(X = X, ncomp = ncomp, ntry = 20, ncores = nc, verbose = F)
cica <- getGO(cica, db = c("BP", "CC", "MF"), ncores = 4)
saveRDS(cica, paste0(heap_folder, "cica", ncomp, "_", postfix, "_", slide_id, "_",
                     nrow(X), "x", ncol(X), ".rds"))

saveReport(cica, file = paste0(output_folder, "report_cica", ncomp, "_", postfix, "_",
                               slide_id, "_", nrow(X), "x", ncol(X), ".pdf"))

# cica <- readRDS(paste0(heap_folder, "cica", ncomp, "_", postfix, "_", slide_id, "_", nrow(X), "x", ncol(X), ".rds"))

spots <- intersect(colnames(cica$M), rownames(coord))

## Diagnostics: depth QC + spatial autocorrelation ----
# depth_cor: cor with total counts
# moran_I  : global Moran's I, kNN (k=8) row-standardised weights (spdep)
coords <- as.matrix(coord[spots, c("x", "y")])
lw <- nb2listw(knn2nb(knearneigh(coords, k = 8)), style = "W", zero.policy = TRUE)
tc <- colSums(counts)[spots]

# perm_p: permutation null p per component (shuffle activity, recompute Moran, n=500)
n_perm <- 500L; S0 <- Szero(lw); n0 <- length(spots)
perm_p <- function(ic) {
  obs <- moran(cica$M[ic, spots], lw, n0, S0, zero.policy = TRUE)$I
  set.seed(2026)
  null_I <- replicate(n_perm, moran(sample(cica$M[ic, spots]), lw, n0, S0, zero.policy = TRUE)$I)
  (sum(null_I >= obs) + 1) / (n_perm + 1)
}

moran_res <- lapply(seq_len(ncomp), function(ic){moran.test(cica$M[ic, spots], lw, zero.policy = TRUE)})
ics_stat <- data.frame(IC = paste0("IC", seq_len(ncomp)),
                          depth_cor = round(vapply(seq_len(ncomp), function(ic) cor(cica$M[ic, spots], tc), numeric(1)), 2),
                          moran_I   = round(vapply(moran_res, function(m) m$estimate[[1]], numeric(1)), 3),
                          moran_p   = signif(vapply(moran_res, function(m) m$p.value, numeric(1)), 3),
                          perm_p    = signif(vapply(seq_len(ncomp), perm_p, numeric(1)), 3))
ics_stat <- ics_stat[order(-ics_stat$moran_I), ]
print(ics_stat)
write.csv(ics_stat, paste0(output_folder, "ics_stat_", slide_id, "_", postfix, "_k", ncomp, ".csv"), row.names = FALSE)

# Moran's I per component
p_moran <- ggplot(ics_stat, aes(reorder(IC, as.integer(sub("IC", "", IC))), moran_I)) +
  geom_segment(aes(xend = IC, y = 0, yend = moran_I), colour = "#522081") +
  geom_point(colour = "#522081", size = 3) +
  labs(x = NULL, y = "Moran's I",
       title = sprintf("%s: spatial autocorrelation per component", slide_id)) +
  theme_minimal(base_size = 12)

## Annotation export: genes (both tails) + GO:BP long (both tails, classicFisher+FDR) ----
source("scripts/demo/annotation_export.R")
# the export sizes are set in annotation_export.R, not here
write_annotation(cica, slide_id, postfix, ncomp, output_folder)

## Spatial component activity maps ----
# biological labels from the annotation CSV (fill `label`, then re-read here)
ann <- read.csv(paste0(output_folder, "annotation_", slide_id, "_", postfix, "_k", ncomp, ".csv"),
                stringsAsFactors = FALSE)
has_label <- !is.na(ann$label) & ann$label != ""
ic_lab <- setNames(ifelse(has_label, paste0(ann$IC, ": ", ann$label), ann$IC), ann$IC)

# Plot all components
winsorize <- T   # winsorised: +/-3 SD
df <- do.call(rbind, lapply(seq_len(ncomp), function(ci)
  data.frame(coord[spots, ],
             IC       = factor(paste0("IC", ci), paste0("IC", seq_len(ncomp))),
             # activity = cica$M[ci, spots],                    # raw weights, one shared scale
             activity = as.vector(scale(cica$M[ci, spots])))))  # per-component z-score
message(sprintf("%s: %.2f%% of unit-values clipped at |z|>3", slide_id,
                100 * mean(abs(df$activity) > 3, na.rm = TRUE)))
p_maps <- ggplot(df, aes(x, y, colour = activity)) +
  geom_point(size = 1) +
  scale_colour_gradient2(low = "#3c5488", mid = "grey92", high = "#9e2b32", midpoint = 0,
                         limits = if (winsorize) c(-3, 3) else NULL,
                         oob = scales::squish, name = "z-score") +
  coord_fixed() + scale_y_reverse() +
  facet_wrap(~ IC, nrow = 2) + theme_void() +
  theme(legend.box.spacing = grid::unit(25, "pt"),
        panel.spacing      = grid::unit(10, "pt")) +
  ggtitle(sprintf("%s: consICA component activity (ncomp=%d, %s)", slide_id, ncomp, postfix))

# Plot single component
plot_component <- function(ic, flip = 1, size = 2.5) {
  ggplot(data.frame(coord[spots, ], activity = flip * cica$M[ic, spots]),
         aes(x, y, colour = activity)) +
    geom_point(size = size) +
    scale_colour_gradient2(low = "#3c5488", mid = "grey92", high = "#9e2b32", midpoint = 0, name = "weight") +
    coord_fixed() + scale_y_reverse() + theme_void() +
    ggtitle(sprintf("%s: %s", slide_id, ic_lab[paste0("IC", ic)]))
}
plot_component(ic = 6)

## Export figures ----
fig <- function(tag) paste0(output_folder, "fig_", tag, "_", slide_id, "_", postfix, "_k", ncomp, ".pdf")
ggsave(fig(if (winsorize) "maps_winz" else "maps"), p_maps, width = 200, height = 200, units = "mm", dpi = 300)
ggsave(fig("moran"), p_moran,           width = 120, height = 100,  units = "mm", dpi = 300)
ggsave(fig("IC6"),   plot_component(6), width = 200, height = 200, units = "mm", dpi = 300)
ggsave(fig("IC7"),   plot_component(7), width = 200, height = 200,  units = "mm", dpi = 300)
ggsave(fig("IC3"),   plot_component(3), width = 200, height = 200,  units = "mm", dpi = 300)

## Save plot objects ----
saveRDS(c(list(maps = p_maps, moran = p_moran),
          setNames(lapply(seq_len(ncomp), plot_component), paste0("IC", seq_len(ncomp)))),
        paste0(output_folder, "fig_obj/figobj_", slide_id, "_", postfix, "_k", ncomp, ".rds"))

