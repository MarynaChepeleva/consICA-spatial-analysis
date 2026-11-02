# Figures for the xenograft section, consICA at k = 12
# main: composition of the component tails, then the two weight maps
# supplementary 1: human share of every tail, then all twelve components mapped
# supplementary 2: the marker-carrying component of every section

library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)
library(mgcv)

# IRanges and AnnotationDbi, attached by consICA, mask the dplyr select() and
# slice() generics, so the dplyr versions are bound here.
select <- dplyr::select
slice <- dplyr::slice
rename <- dplyr::rename


# ---- settings ----

cica_dir <- "heap/cica_ntry20"
f_counts <- "heap/counts_mat.rds"
f_seurat <- "heap/merged_seurat.rds"
f_species <- "heap/gene_species.rds"
f_tumour <- "data/tumor_spots.csv"
out_dir <- "output/ST_paper"

ncomp <- 12
ref_slide <- 6
HAVE_SLIDE1_HM <- FALSE   # have human alignment for slide 1
rep_slides <- if (HAVE_SLIDE1_HM) c(6, 5, 7, 8, 1) else c(6, 5, 7, 8)   # 1 is the tumour-free control

z_thresh <- 3     # |z| on the loadings that defines a component side
min_side <- 10    # smallest side kept
rot_maps <- 90    # rotation of the maps, degrees clockwise
mirror <- TRUE  # mirror the maps about the vertical axis to match the H&E

# Dot size in millimetres, one per figure, set against the panel width each is
# printed at rather than against the canvas it is composed on. The five replicate
# panels are the smallest cells and take the smallest dot.
pt_main <- 0.4
pt_all <- 0.25
pt_rep <- 0.4

# which component plays which role is decided in xeno_numbers.R by behaviour and
# read from there, so a rerun that renumbers the components cannot mislabel them
f_roles <- file.path(out_dir, "components.csv")

# the melanocytic set already used in the manuscript; the same list is used here
# so the two sections are comparable
markers <- c("PMEL", "TYR", "TYRP1", "MLANA", "SOX10", "TRPM1", "DCT")

# The diverging scale is the one the demonstration figures use, from
# consICA_Visium.R and consICA_Xenium.R, so a weight map reads the same
# everywhere in the paper.
pal <- c(prog = "#754E97", abun = "#C46F45", shared = "#5B9A7A",
         other = "#9AA6B0", tumour = "#C85C6C", rest = "#9AA6B0",
         neg = "#3c5488", zero = "grey92", pos = "#9e2b32",
         # the two tails of the bar panel: the same two hues as the weight
         # scale, at lightness 48 against 38 at the poles
         bar_neg = "#64739D", bar_pos = "#BC6864")

# Type sizes, panel letters and the black-text rule follow those of figures_pre.R,
# so both halves of the paper set their figures alike.
# Type size is carried by base_size on each theme: BASE sets the axis and legend
# titles, and the ticks and legend labels come out at 0.8 of it.
pt <- function(p) p / 2.845                     # pt -> ggplot text size (mm)
BASE <- 8
TYPE <- list(letter = 9, title = 7.5, sub = 6.5)

# The right margin is the same on every panel and wide enough for a colour bar.
# The bars are drawn into it as legends past the panel edge, so they are not
# taken out of the width the panel is given.
LEG_PAD <- 34            # points reserved to the right of every panel
SCAT_W <- c(1.3, 1)     # the composition panel with its legend, against the space left beside it
COL_W <- c(1, 1)         # relative widths of the two maps
ROW_H <- c(0.8, 1)         # relative heights of the two rows

# COL_W stays equal: the two maps share one colour scale. All three panels are of
# fixed aspect, so grid locks width to height and resizing one panel resizes the
# whole figure. The levers are FIG_W and these three vectors.

# The letter goes in the top margin, not in the plot area, where it would land on
# the axis labels. The position has to be named rather than given as coordinates,
# numeric coordinates being read relative to the plot.
letter_theme <- theme(
  plot.tag = element_text(size = TYPE$letter, face = "bold", hjust = 0, vjust = 1),
  plot.tag.position = "topleft",
  plot.tag.location = "margin",
  plot.margin = margin(t = 12, r = LEG_PAD, b = 3, l = 10))

# Of everything theme_minimal draws, only the tick labels and the facet strips
# come out grey. Naming an axis element is why this cannot go on the maps, which
# are drawn on theme_void.
black_text <- theme(axis.text = element_text(colour = "black"),
                    strip.text = element_text(colour = "black"))

# short white ticks at the ends of a colour bar rather than one solid rule
legend_ticks <- theme(
  legend.ticks = element_line(colour = "white", linewidth = 0.2),
  legend.ticks.length = unit(0.06, "cm"))


# ---- functions ----

# the two sets of extreme loadings of a component, by sign
get_sides <- function(S, ic) {
  s <- as.vector(scale(S[, ic]))
  list("+" = rownames(S)[s > z_thresh],
       "-" = rownames(S)[s < -z_thresh])
}

# naive orthologue map, GAPDH to Gapdh; renamed pairs such as HBB and Hbb-bs
# are missed, so the paired share is a slight underestimate
to_mouse <- function(g) paste0(substr(g, 1, 1), tolower(substr(g, 2, nchar(g))))

pct_human <- function(genes) {
  if (!length(genes)) return(NA)
  100 * mean(genes %in% human_genes)
}

# of the human genes of a side, the share whose mouse counterpart is in the same
# side; the counterpart has to be a gene that really came from the mouse matrix
pct_paired <- function(genes) {
  h <- intersect(genes, human_genes)
  if (!length(h)) return(NA)
  partner <- to_mouse(h)
  100 * mean(partner %in% intersect(genes, mouse_genes))
}

# composition of every side of every component
side_table <- function(S) {
  lapply(colnames(S), function(ic) {
    sides <- get_sides(S, ic)
    data.frame(ic = ic,
               side = names(sides),
               n = lengths(sides),
               human = sapply(sides, pct_human),
               paired = sapply(sides, pct_paired),
               markers = sapply(sides, function(g) sum(markers %in% g)),
               row.names = NULL,
               stringsAsFactors = FALSE)
  }) %>%
    bind_rows() %>%
    filter(n >= min_side)
}

# The same rule as in xeno_numbers.R: among the tails holding at least one gene
# of the melanocytic set, the most human one. Ordering by marker count instead
# picks the largest tail, which is the abundance axis.
pick_programme <- function(tab) {
  hits <- tab %>% filter(markers > 0)
  if (!nrow(hits)) hits <- tab
  # dplyr:: qualified because IRanges, pulled in by consICA and the annotation
  # databases, defines its own slice generic and masks this one
  hits %>% arrange(desc(human)) %>% dplyr::slice(1)
}

cica_path <- function(slide, k = ncomp) {
  hits <- list.files(cica_dir, full.names = TRUE,
                     pattern = sprintf("^cica%d_slide_%d_[0-9]+x[0-9]+\\.rds$",
                                       k, slide))
  if (!length(hits)) stop("no decomposition for slide ", slide)
  hits[1]
}

slide_image <- function(slide) {
  if (slide == 1) "slice1" else paste0("slice1.", slide)
}

# Spot coordinates placed the same way round as the stored H&E. Indexed by column
# rather than by [[ ]]: on a matrix [[1]] is the first element, not the first
# column, which silently collapsed every spot onto one point.
place_xy <- function(xy) {
  stopifnot(ncol(xy) == 2, nrow(xy) > 1)
  x <- xy[, 1]
  y <- xy[, 2]
  out <- switch(as.character(rot_maps %% 360),
                "0"   = cbind(x, y),
                "90"  = cbind(-y, x),
                "180" = cbind(-x, -y),
                "270" = cbind(y, -x))
  if (mirror) out[, 1] <- -out[, 1]
  out
}

plot_he <- function(img, title = NULL) {
  ggplot() +
    annotation_raster(img, 0, ncol(img), 0, nrow(img)) +
    scale_x_continuous(limits = c(0, ncol(img)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, nrow(img)), expand = c(0, 0)) +
    coord_fixed() +
    labs(title = title) +
    theme_void(base_size = BASE) +
    theme(plot.title = element_text(size = TYPE$title, hjust = 0.5))
}

# leg_title names the colour bar, so the panel title is free to say what the
# component is. shape 19 carries a border, which puts a floor of about a third of
# a millimetre under the drawn dot, and shape 16 has none. breaks are left to
# breaks_pretty on the raw scales and named explicitly on a standardised one.
plot_map <- function(weight, xy, title, lim, pt_size = 0.7, tumour = NULL,
                     leg_title = "weight", shape = 19,
                     breaks = scales::breaks_pretty(n = 4), subtitle = NULL) {
  df <- data.frame(x = xy[, 1], y = xy[, 2], weight = weight)
  p <- ggplot(df, aes(x, y, colour = weight)) +
    geom_point(size = pt_size, shape = shape)
  # black rings mark the spots annotated as tumour, an
  # annotation the decomposition never saw
  if (!is.null(tumour)) {
    p <- p + geom_point(data = df[tumour, ], shape = 21, size = pt_size * 1.4,
                        colour = "grey10", stroke = 0.2, fill = NA)
  }
  p +
    scale_colour_gradient2(low = pal[["neg"]], mid = pal[["zero"]],
                           high = pal[["pos"]], midpoint = 0,
                           limits = c(-lim, lim), oob = scales::squish,
                           breaks = breaks, name = leg_title) +
    scale_y_reverse() +
    coord_fixed(clip = "off") +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_size = BASE) +
    theme(plot.title = element_text(size = TYPE$title, hjust = 0.5),
          # a second line of its own element rather than of the title, so that it
          # can be set smaller: a title carries one size for all its lines
          plot.subtitle = element_text(size = TYPE$sub, hjust = 0.5,
                                       colour = "black"),
          legend.key.width = unit(0.25, "cm")) +
    legend_ticks
}


# ---- data ----

roles <- read.csv(f_roles, stringsAsFactors = FALSE)
ic_prog <- roles$ic[roles$role == "programme"]
ic_abun <- roles$ic[roles$role == "abundance"]
ic_shared <- roles$ic[roles$role == "orthologue pairs"]
cat("roles from", f_roles, ":", ic_prog, ic_abun, ic_shared, "\n")

# Reading the count matrix and the Seurat object is the slow part of this script.
# A second source() in the same session reuses them, and rm(he_list) forces a
# fresh read. How the coordinates are placed is settled further down, so a change
# to rot_maps or mirror takes effect without reloading.
if (exists("he_list")) {
  cat("reusing the matrices and images already in the session\n")
} else {

cica <- readRDS(cica_path(ref_slide))
S <- cica$S
M <- cica$M
spots <- colnames(M)

# species comes from which of the two Space Ranger matrices a row was merged
# from, recorded gene by gene in the lookup rather than guessed from the symbol
gene_species <- readRDS(f_species)[rownames(S)]
stopifnot(!anyNA(gene_species))
human_genes <- names(gene_species)[gene_species == "human"]
mouse_genes <- names(gene_species)[gene_species == "mouse"]

counts_slide <- readRDS(f_counts)[, spots]
human_frac <- Matrix::colSums(counts_slide[human_genes, ]) /
  Matrix::colSums(counts_slide)
rm(counts_slide)
gc()

tumour_barcodes <- as.character(read.csv(f_tumour)[[1]])
is_tumour <- spots %in% tumour_barcodes

# The stored object predates the installed Seurat, so its slots are migrated on
# load. The migration report and its warnings are silenced here and nowhere else.
merged_seurat <- suppressWarnings(
  suppressMessages(UpdateSeuratObject(readRDS(f_seurat))))
coords_list <- lapply(rep_slides, function(s) {
  xy <- GetTissueCoordinates(merged_seurat, image = slide_image(s))
  cols <- if (all(c("imagecol", "imagerow") %in% colnames(xy))) {
    c("imagecol", "imagerow")
  } else {
    c("y", "x")
  }
  xy[, cols]
})
he_list <- lapply(rep_slides, function(s) {
  as.raster(merged_seurat@images[[slide_image(s)]]@image)
})
names(coords_list) <- as.character(rep_slides)
names(he_list) <- as.character(rep_slides)
rm(merged_seurat)
gc()

}

# genes absent from the matrix could never be found in any tail, so they are not
# counted against; TYRP1 and TRPM1 fall out at the detection filter
markers <- intersect(markers, rownames(S))
cat("marker set in the matrix:", paste(markers, collapse = " "), "\n")

coords_ref <- place_xy(coords_list[[as.character(ref_slide)]][spots, ])
he_ref <- he_list[[as.character(ref_slide)]]


# ---- how a component, a tail and a section are named ----
# One form throughout, the manuscript's own: IC12, the tail as a typographic
# minus or plus, and XG6 for the section. The internal names stay as they are and
# are converted only where something is drawn.

ic_name <- function(ic) paste0("IC", sub("^ic\\.", "", ic))
tail_name <- function(side) ifelse(side == "-", "neg", "pos")
xg_name <- function(slide) paste0("XG", slide)


# ---- composition of the component sides, reference slide ----

comp_sides <- side_table(S) %>%
  mutate(label = paste(ic_name(ic), tail_name(side)),
         role = case_when(ic == ic_prog & human > 50 ~ "prog",
                          ic == ic_abun & human > 50 ~ "abun",
                          ic == ic_shared & human > 40 ~ "shared",
                          TRUE ~ "other"))

# which tail of a component holds its human genes; the sign is arbitrary and
# differs between runs, so panel titles take it from here rather than assume it
top_side <- comp_sides %>%
  group_by(ic) %>%
  slice_max(human, n = 1) %>%
  ungroup()
human_side <- setNames(top_side$side, top_side$ic)


# ---- main figure panels ----

# built but not placed in the figure
p_he <- plot_he(he_ref)

# clip = "off" lets a circle sitting against the top of the panel be drawn whole.
# The size key is drawn in the colour of the unlabelled tails, a black key having
# outweighed the grey points it stands for.
p_comp <- ggplot(comp_sides, aes(human, paired)) +
  geom_point(aes(size = n, colour = role), alpha = 0.85) +
  geom_text(data = filter(comp_sides, role != "other"),
            aes(label = label, colour = role), size = pt(7),
            hjust = 0.5, vjust = 0, nudge_y = 6, show.legend = FALSE) +
  scale_colour_manual(values = pal[c("prog", "abun", "shared", "other")],
                      breaks = c("prog", "abun", "shared", "other"),
                      guide = "none") +
  scale_size_area(max_size = 4.5, name = "genes in tail",
                  guide = guide_legend(
                    override.aes = list(colour = pal[["other"]]))) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  coord_fixed(ratio = 1, clip = "off") +
  labs(x = "human genes in the tail (%)",
       y = "human genes whose mouse\ncounterpart is in the same tail (%)") +
  theme_minimal(base_size = BASE) +
  black_text +
  theme(panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.4),
        # outside the panel, to its right, the column aligned to the top. Drawn
        # outside it is taken out of the width of the cell, and the panel is
        # square, so SCAT_W widens that cell by about the width of the legend to
        # keep the square the size it was.
        legend.position = "right",
        legend.justification = "top")

# the per-spot human read fraction locates the tumour without any decomposition,
# so the panel asks how much of a component that statistic already accounts for
# the same order as the maps above, programme first
redundancy <- rbind(
  data.frame(ic = ic_prog, weight = M[ic_prog, spots],
             frac = human_frac, tumour = is_tumour),
  data.frame(ic = ic_abun, weight = M[ic_abun, spots],
             frac = human_frac, tumour = is_tumour)) %>%
  mutate(ic = factor(ic, levels = c(ic_prog, ic_abun),
                     labels = ic_name(c(ic_prog, ic_abun))),
         region = ifelse(tumour, "annotated tumour", "rest of the section"))

# deviance explained by a smooth allows any shape; no line is fitted and no
# slope is interpreted
gam_fit <- redundancy %>%
  group_by(ic) %>%
  summarise(dev = summary(gam(weight ~ s(frac, bs = "cs"),
                              method = "REML"))$dev.expl,
            r2 = cor(weight, frac)^2,
            rho = cor(weight, frac, method = "spearman"),
            .groups = "drop") %>%
  mutate(label = sprintf("deviance explained = %.1f%%", 100 * dev))

# One component per panel rather than two facets of one: a facetted panel divides
# a single cell of the layout, so its split cannot line up with the columns of
# the grid. The legend goes beneath each panel, where it is taken out of the
# height and not the width.
plot_redundancy <- function(lab) {
  ggplot(redundancy[redundancy$ic == lab, ], aes(frac, weight)) +
    # shape 16 rather than the default 19: 19 draws a border round the dot, which
    # at this size is most of what the dot is
    geom_point(aes(colour = region), size = 0.5, alpha = 0.5, shape = 16) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                method.args = list(method = "REML"),
                colour = "black", fill = "grey60", alpha = 0.3,
                linewidth = 0.4) +
    geom_text(data = gam_fit[gam_fit$ic == lab, ], aes(label = label),
              x = -Inf, y = Inf, hjust = -0.08, vjust = 1.1, size = pt(7),
              inherit.aes = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
    scale_colour_manual(values = c("annotated tumour" = pal[["tumour"]],
                                   "rest of the section" = pal[["rest"]]),
                        name = NULL) +
    # the key is enlarged but not made opaque: the annotated spots are few and do
    # not overlap, so at full opacity the key showed a rose that appears nowhere in
    # the panel, while the grey, which is thousands of overlapping points, does
    # stack to near solid
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 0.7))) +
    labs(x = "human read fraction of the spot", y = paste(lab, "weight")) +
    theme_minimal(base_size = BASE) +
    black_text +
    # inside the panel, anchored by its top right corner and dropped a little from
    # the top edge; the corner is empty on both components, and a legend drawn
    # inside costs the layout nothing
    theme(panel.grid = element_blank(),
          axis.line = element_line(linewidth = 0.4),
          legend.position = "inside",
          legend.position.inside = c(1.12, 0.90),
          legend.justification = c(1, 1),
          legend.direction = "vertical",
          legend.key.spacing.y = unit(0, "pt"),
          legend.key.height = unit(10, "pt"),
          legend.background = element_rect(fill = "white", colour = NA))
}


# ---- main figure ----

lim_main <- as.numeric(quantile(abs(M[c(ic_prog, ic_abun), spots]), 0.99))

# The map carries no title: the colour bar names the component, and what the
# component is belongs in the caption. The bar sits past the right edge of the
# panel, in the margin every panel reserves, so the map keeps its full width.
map_grid <- function(ic) {
  plot_map(M[ic, spots], coords_ref, NULL, lim_main, pt_size = pt_main,
           tumour = is_tumour, leg_title = ic_name(ic)) +
    # For a colour bar ggplot reads legend.key.height as the length of one
    # division and multiplies it by about five, so the value here is a fifth of
    # the bar wanted: 0.3 cm gives a bar of roughly 1.5 cm against a map of four.
    theme(legend.position = "inside",
          legend.position.inside = c(1.02, 0.5),
          legend.justification = c(0, 0.5),
          legend.key.height = unit(0.3, "cm"))
}

# One flat grid of six cells: the composition panel alone on the first row, the
# two maps side by side on the second. The stain and the two redundancy panels
# are built above but not placed. Letters are set panel by panel, because
# tag_levels would letter the row rather than the panel.
fig_main <- ((p_comp + labs(tag = "a") | plot_spacer()) +
               plot_layout(widths = SCAT_W)) /
  ((map_grid(ic_prog) + labs(tag = "b")) |
     (map_grid(ic_abun) + labs(tag = "c"))) +
  plot_layout(heights = ROW_H) &
  letter_theme

# All three panels carry a fixed aspect, so how much of its cell each one fills
# depends on it. Printed because the width of the figure and the three numbers of
# COL_W are what trade their whitespace against one another.
cat(sprintf("aspect (width/height): weight map %.2f, composition panel 1.00\n",
            diff(range(coords_ref[, 1])) / diff(range(coords_ref[, 2]))))


# ---- supplementary 1: the composition of every tail, and every component mapped ----

BAR_W <- c(2.4, 1)   # the bar panel against the space left empty beside it

p_bars <- ggplot(comp_sides,
                 aes(human, reorder(ic_name(ic), human, FUN = max),
                     fill = side)) +
  # dodge rather than dodge2: dodge2 packs the bars into the space a category has,
  # and at the first and last category that space is cut by the expansion of the
  # discrete axis, which pushed those two pairs apart. preserve keeps the single
  # bar of the one component with only one tail the same width as the rest.
  geom_col(position = position_dodge(width = 0.7, preserve = "single"),
           width = 0.7) +
  scale_fill_manual(values = c("+" = pal[["bar_pos"]], "-" = pal[["bar_neg"]]),
                    labels = tail_name, name = "tail") +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = "human genes in the tail (%)", y = NULL) +
  theme_minimal(base_size = BASE) +
  black_text +
  theme(panel.grid = element_blank(),
        axis.line.x = element_line(linewidth = 0.4),
        legend.position = "right")

# Standardised per component and clipped at three standard deviations, as in the
# other supplementary grids, so the twelve are comparable and a few extreme units
# cannot set the scale. The main figure stays on raw weights.
Z_CLIP <- 3

z_all <- lapply(rownames(M), function(ic) as.vector(scale(M[ic, spots])))
names(z_all) <- rownames(M)
cat(sprintf("component grid: %.2f%% of weights clipped at |z| > %d\n",
            100 * mean(abs(unlist(z_all)) > Z_CLIP), Z_CLIP))

# One facetted plot rather than twelve assembled ones. Assembled, each map was a
# leaf of the layout and took a letter of its own, so the twelve ran from c to n;
# facetted they are one panel with one letter and one colour bar, and the component
# names move into the strips.
maps_all <- bind_rows(lapply(rownames(M), function(ic) {
  data.frame(ic = ic_name(ic), x = coords_ref[, 1], y = coords_ref[, 2],
             z = pmax(pmin(z_all[[ic]], Z_CLIP), -Z_CLIP), row.names = NULL)
})) %>%
  mutate(ic = factor(ic, levels = ic_name(rownames(M))))

p_all <- ggplot(maps_all, aes(x, y, colour = z)) +
  # shape 16 has no border, so pt_all is the whole of the dot; the default 19 puts
  # a floor of about a third of a millimetre under it
  geom_point(size = pt_all, shape = 16) +
  facet_wrap(~ ic, ncol = 4) +
  scale_colour_gradient2(low = pal[["neg"]], mid = pal[["zero"]],
                         high = pal[["pos"]], midpoint = 0,
                         limits = c(-Z_CLIP, Z_CLIP),
                         name = "z-score") +
  scale_y_reverse() +
  coord_fixed() +
  theme_void(base_size = BASE) +
  theme(strip.text = element_text(size = TYPE$title, colour = "black"),
        legend.key.width = unit(0.25, "cm")) +
  legend_ticks

# The bar panel is given part of the row and the rest left empty, which keeps it
# from stretching over the whole width. Letters are set by hand because the first
# row is a nested layout, where tag_levels would letter the row and not the panel.
fig_supp <- (((p_bars + labs(tag = "a")) | plot_spacer()) +
               plot_layout(widths = BAR_W)) /
  (p_all + labs(tag = "b")) +
  plot_layout(heights = c(1, 2.4)) &
  letter_theme


# ---- supplementary 2: the marker-carrying component of every section ----
# On each section the component is chosen by the rule in xeno_numbers.R, without
# reference to which component it was on the reference section. Slide 1 carries
# no tumour cells, so the rule returns its most human tail.

# the five decompositions are read again here, so this block is reused on a
# second source() as well; the coordinates are kept unplaced for the same reason
# as above
if (exists("rep_data")) {
  cat("reusing the five decompositions already in the session\n")
} else {
rep_data <- lapply(rep_slides, function(s) {
  obj <- readRDS(cica_path(s))
  tab <- pick_programme(side_table(obj$S))
  out <- list(
    slide = s, tab = tab,
    xy = coords_list[[as.character(s)]][colnames(obj$M), ],
    weight = obj$M[tab$ic, ],
    tumour = colnames(obj$M) %in% tumour_barcodes,
    # the genes themselves, not a count: a larger tail collects more of any
    # list, so a count makes the purest tail look the weakest
    found = intersect(markers, get_sides(obj$S, tab$ic)[[tab$side]]))
  rm(obj)
  gc()
  out
})
}

rep_summary <- bind_rows(lapply(rep_data, function(d) cbind(slide = d$slide, d$tab)))

# Standardised within each component and clipped, not put on one raw scale. A
# component of an ICA decomposition is determined only up to a scalar, S times c
# and M over c giving the same product, so the size of a weight carries no fixed
# meaning and weights from five separate decompositions cannot share a colour bar.
z_rep <- lapply(rep_data, function(d) as.vector(scale(d$weight)))
cat(sprintf("replicate maps: %.2f%% of weights clipped at |z| > %d\n",
            100 * mean(abs(unlist(z_rep)) > Z_CLIP), Z_CLIP))

# The title carries the section, the component with its tail, and the melanocytic
# genes the tail holds. The percentages of human and of paired genes are in Table 3
# of the manuscript and are not repeated. Where the tail holds no marker, a dash
# rather than nothing: an empty line reads as something missing, not as a count of
# zero, and on the control section that count is the finding.
panels_rep <- lapply(seq_along(rep_data), function(i) {
  d <- rep_data[[i]]
  plot_map(pmax(pmin(z_rep[[i]], Z_CLIP), -Z_CLIP), place_xy(d$xy),
           sprintf("%s, %s %s", xg_name(d$slide), ic_name(d$tab$ic),
                   tail_name(d$tab$side)),
           Z_CLIP, pt_size = pt_rep, tumour = d$tumour, leg_title = "z-score",
           shape = 16, breaks = waiver(),
           subtitle = if (length(d$found)) paste(d$found, collapse = " ") else "—")
})

# four sections that received tumour cells on the first row, the control alone on
# the second, in the order of rep_slides; the empty cells hold the control to the
# width of the panels above it
fig_rep <- wrap_plots(panels_rep, guides = "collect",
                      design = "ABCD\nE###")


# ---- save ----

FIG_W <- 17.5   # cm

save_fig <- function(p, file, height_cm) {
  ggsave(file.path(out_dir, file), p, width = FIG_W, height = height_cm,
         units = "cm", device = grDevices::cairo_pdf, limitsize = FALSE)
  cat(sprintf("%-32s %.1f x %.1f cm\n", file, FIG_W, height_cm))
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
save_fig(fig_main, "Figure3.pdf", 13.0)
save_fig(fig_supp, "FigureS3.pdf", 18.3)
save_fig(fig_rep, "FigureS4.pdf", 10.5)

cat("figures written to", out_dir, "\n")


# ---- values quoted in the text ----

cat("annotated tumour spots:", sum(is_tumour), "of", length(spots), "\n")
cat("human read fraction: median", round(median(human_frac), 3),
    "max", round(max(human_frac), 3), "\n")
cat("composition panel:", nrow(comp_sides), "sides of", 2 * nrow(M), "\n\n")
print(as.data.frame(arrange(comp_sides, desc(human))), digits = 3,
      row.names = FALSE)
cat("\n")
print(as.data.frame(gam_fit[, c("ic", "dev", "r2", "rho")]), digits = 4,
      row.names = FALSE)
cat("\nmarker-carrying component per section\n")
print(rep_summary, digits = 3, row.names = FALSE)
