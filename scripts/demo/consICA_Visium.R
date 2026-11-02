# Visium demo: load, preprocess, apply consICA, report and spatial component maps.
# Writes the annotation tables, ics_stat, the figures and the saved plot objects
# that the assembly scripts read.
# One slide per run: set slide_id under Slide params.

library(SingleCellExperiment)
library(zellkonverter)
library(Seurat)
library(consICA)
library(ggplot2)
library(spdep)

options(future.globals.maxSize = 4 * 1024^3)   # for SCTransform on full transcriptome

## Slide params ----

HEST_DIR <- Sys.getenv("HEST_ROOT", "data/HEST")
SLIDES <- list(MEND40  = "skin/st/MEND40.h5ad",
               NCBI776 = "breast/st/NCBI776.h5ad")

slide_id <- "NCBI776"   # <- set the slide to run
input_file <- file.path(HEST_DIR, SLIDES[[slide_id]])

ncomp <- 20              # full transcriptome supports more components than a panel
nc <- 1                  # consICA cores; 1 = reproducible + avoids OOM on large Visium
getgo_ncores <- 4        # getGO cores; NCBI776 can OOM at 4, drop to 1 if it does

norm_method <- "SCT"     # input to consICA: "SCT" | "LogNorm" | "log2"
filter_spots <- FALSE    # count-based per-spot QC; FALSE = keep all spots
postfix <- norm_method   # tag for output filenames

# Xenium control pattern — harmless on Visium (matches no real gene symbols)
ctrl_pattern <- "^(BLANK|NegControlProbe|NegControlCodeword|antisense|DeprecatedCodeword|UnassignedCodeword|Intergenic)"


## Output folders and functions ----
output_folder <- "output/demo/"
heap_folder <- "output/heap/demo/"
for (d in c(output_folder, heap_folder)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

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

# Count-based per-spot QC (full-tissue Visium is usually unimodal -> floor cut only)
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
    geom_point(size = 1.2) +
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

tc <- colSums(GetAssayData(srt, assay = "originalexp", layer = "counts"))
df_tc <- data.frame(coord[colnames(srt), ], total_counts = tc)
print(
  ggplot(df_tc, aes(x, y, colour = total_counts)) +
    geom_point(size = 1.2) + scale_colour_viridis_c(trans = "log10") +
    coord_fixed() + scale_y_reverse() + theme_void() +
    ggtitle(sprintf("%s — spots by total counts (n=%d)", slide_id, ncol(srt)))
)

## Normalization ----
counts <- GetAssayData(srt, assay = "originalexp", layer = "counts")   # final (post-QC) counts
if (norm_method == "SCT") { # Pearson residuals
  srt <- suppressWarnings(SCTransform(srt, assay = "originalexp",
                                      return.only.var.genes = FALSE, verbose = FALSE))
  X <- as.matrix(srt[["SCT"]]@scale.data)
} else if (norm_method == "LogNorm") { # ln(CP10k+1)
  srt <- NormalizeData(srt, assay = "originalexp",
                       normalization.method = "LogNormalize", scale.factor = 1e4)
  X <- as.matrix(GetAssayData(srt, assay = "originalexp", layer = "data"))
} else if (norm_method == "log2") { # no depth normalization
  X <- as.matrix(log2(counts + 1))
} else stop("Unknown norm_method: ", norm_method)
dim(X)


## Signal deconvolution with consICA ----
set.seed(2026)
cica <- consICA(X = X, ncomp = ncomp, ntry = 20, ncores = nc, verbose = F)
cica <- getGO(cica, db = c("BP", "CC", "MF"), ncores = getgo_ncores)
saveRDS(cica, paste0(heap_folder, "cica", ncomp, "_", postfix, "_", slide_id, "_",
                     nrow(X), "x", ncol(X), ".rds"))

saveReport(cica, file = paste0(output_folder, "report_cica", ncomp, "_", postfix, "_",
                               slide_id, "_", nrow(X), "x", ncol(X), ".pdf"))

cica <- readRDS(paste0(heap_folder, "cica", ncomp, "_", postfix, "_", slide_id, "_",
                       nrow(X), "x", ncol(X), ".rds"))

spots <- intersect(colnames(cica$M), rownames(coord))

## Diagnostics: depth QC + spatial autocorrelation ----
# depth_cor: cor with total counts
# moran_I  : global Moran's I, kNN (k=6, hex Visium) row-standardised weights (spdep)
# moran_p  : analytic p ; perm_p: permutation null p (shuffle activity, n=500)
coords <- as.matrix(coord[spots, c("x", "y")])
lw <- nb2listw(knn2nb(knearneigh(coords, k = 6)), style = "W", zero.policy = TRUE)
tc <- colSums(counts)[spots]

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
  geom_point(colour = "#522081", size = 2) +
  labs(x = NULL, y = "Moran's I",
       title = sprintf("%s: spatial autocorrelation per component", slide_id)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

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

# Plot all components. Colour = per-component z-score (ICA scale is arbitrary, so raw
# magnitudes are not comparable across panels); winsorize clips the DISPLAY at |z|=3.
winsorize <- TRUE   # FALSE = full z-range
df <- do.call(rbind, lapply(seq_len(ncomp), function(ci)
  data.frame(coord[spots, ],
             IC       = factor(paste0("IC", ci), paste0("IC", seq_len(ncomp))),
             # activity = cica$M[ci, spots],                    # OLD: raw weights, one shared scale
             activity = as.vector(scale(cica$M[ci, spots])))))  # per-component z-score
message(sprintf("%s: %.2f%% of unit-values clipped at |z|>3", slide_id,
                100 * mean(abs(df$activity) > 3, na.rm = TRUE)))
p_maps <- ggplot(df, aes(x, y, colour = activity)) +
  geom_point(size = 0.35) +
  scale_colour_gradient2(low = "#3c5488", mid = "gray92", high = "#9e2b32", midpoint = 0,
                         limits = if (winsorize) c(-3, 3) else NULL,
                         oob = scales::squish, name = "z-score") +
  coord_fixed() + scale_y_reverse() +
  facet_wrap(~ IC) + theme_void() +
  theme(legend.box.spacing = grid::unit(16, "pt"),   # gap between panels and the legend
        panel.spacing      = grid::unit(8, "pt")) +   # gap between panels
  ggtitle(sprintf("%s: consICA component activity (ncomp=%d, %s)", slide_id, ncomp, postfix))

# Plot single component (own colour scale; flip = -1 reverses sign for display)
plot_component <- function(ic, flip = 1, size = 2.5) {
  ggplot(data.frame(coord[spots, ], activity = flip * cica$M[ic, spots]),
         aes(x, y, colour = activity)) +
    geom_point(size = size) +
    scale_colour_gradient2(low = "#3c5488", mid = "gray92", high = "#9e2b32", midpoint = 0, name = "weight") +
    coord_fixed() + scale_y_reverse() + theme_void() +
    ggtitle(sprintf("%s: %s", slide_id, ic_lab[paste0("IC", ic)]))
}
plot_component(ic = 10)

## Export figures for the paper ----
# vector PDF; for raster change extension to .tiff (ggsave picks device from it)
fig <- function(tag) paste0(output_folder, "fig_", tag, "_", slide_id, "_", postfix, "_k", ncomp, ".pdf")
ggsave(fig(if (winsorize) "maps_winz" else "maps"),  p_maps,            width = 300, height = 300, units = "mm", dpi = 300)
ggsave(fig("moran"), p_moran,width = 120, height = 100,  units = "mm", dpi = 300)

# single-component PDFs
ggsave(fig("IC15"), plot_component(15), width = 200, height = 200, units = "mm", dpi = 300)
ggsave(fig("IC18"),  plot_component(18),  width = 200, height = 200, units = "mm", dpi = 300)
ggsave(fig("IC10"), plot_component(10), width = 200, height = 200, units = "mm", dpi = 300)
ggsave(fig("IC7"),  plot_component(7),  width = 200, height = 200, units = "mm", dpi = 300)

ggsave(fig("IC14"),  plot_component(14),  width = 200, height = 200, units = "mm", dpi = 300)
ggsave(fig("IC9"),  plot_component(9),  width = 200, height = 200, units = "mm", dpi = 300)
ggsave(fig("IC11"),  plot_component(11),  width = 200, height = 200, units = "mm", dpi = 300)

ggsave(fig("IC13"), plot_component(13), width = 200, height = 200, units = "mm", dpi = 300)

## Save plot objects for the assembly scripts ----
saveRDS(c(list(maps = p_maps, moran = p_moran),
          setNames(lapply(seq_len(ncomp), plot_component), paste0("IC", seq_len(ncomp)))),
        paste0(output_folder, "figobj_", slide_id, "_", postfix, "_k", ncomp, ".rds"))
