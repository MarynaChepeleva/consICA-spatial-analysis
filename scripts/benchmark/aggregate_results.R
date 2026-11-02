
# Aggregate all benchmark + GO results into tidy master tables and draw core figures.
# Read-only over output/; writes ALL_*.csv summaries and returns ggplot objects.

suppressPackageStartupMessages({
  library(ggplot2)
})

output_dir <- "output/ica_perturbations"
go_dir <- "output/go_enrichment"

# Per-slide labels from config/runs.csv; real POST-FILTER dims from each run's metadata.
runs_cfg <- read.csv("config/runs.csv", stringsAsFactors = FALSE)
slide_meta <- unique(runs_cfg[, c("slide", "platform", "tissue")])
names(slide_meta) <- c("Slide", "Platform", "Tissue")

# dims recorded AFTER spot/gene filtering (metadata$n_genes / n_spots of a pure_stability run)
.slide_dims <- function(output_dir, slide) {
  f <- list.files(output_dir, sprintf("^runs__%s__.*pure_stability.*\\.rds$", slide),
                  full.names = TRUE)
  if (!length(f)) return(c(n_genes = NA_real_, n_spots = NA_real_))
  md <- readRDS(f[1])$metadata
  c(n_genes = md$n_genes, n_spots = md$n_spots)
}
dims <- t(vapply(slide_meta$Slide, function(s) .slide_dims(output_dir, s), numeric(2)))
slide_meta$n_genes <- dims[, "n_genes"]
slide_meta$n_spots <- dims[, "n_spots"]
slide_meta$spots_per_gene <- slide_meta$n_spots / slide_meta$n_genes

# parse "SLIDE__kK__..." tag -> data.frame(Slide, K)
.parse_tag <- function(basename_no_prefix) {
  parts <- strsplit(sub("\\.csv$", "", basename_no_prefix), "__")[[1]]
  k <- suppressWarnings(as.integer(sub("^k", "", parts[grep("^k[0-9]+$", parts)][1])))
  list(Slide = parts[1], K = k)
}

# 1. AGGREGATE ----

aggregate_benchmark <- function(output_dir) {
  files <- list.files(output_dir, "^metrics_FINAL__.*\\.csv$", full.names = TRUE)
  df <- do.call(rbind, lapply(files, function(f) {
    d <- read.csv(f, check.names = FALSE)
    t <- .parse_tag(sub("^metrics_FINAL__", "", basename(f)))
    d$Slide <- t$Slide; d$K <- t$K
    d
  }))
  merge(slide_meta, df, by = "Slide", all.y = TRUE)   # keep rows even for unlisted slides
}

aggregate_go <- function(go_dir) {
  files <- list.files(go_dir, "^go_semantic_uniqueness__.*\\.csv$", full.names = TRUE)
  if (!length(files)) return(NULL)
  df <- do.call(rbind, lapply(files, function(f) {
    d <- read.csv(f, check.names = FALSE)
    t <- .parse_tag(sub("^go_semantic_uniqueness__", "", basename(f)))
    if (is.null(d$Slide)) d$Slide <- t$Slide
    d$K <- t$K
    d
  }))
  merge(slide_meta[, c("Slide", "Platform", "Tissue")], df, by = "Slide")
}

# 2. CORE FIGURES ----

# (a) Diversity headline: Intra_Uniqueness per method, by slide (K=20, pure_stability)
fig_intra_uniqueness <- function(bench, k = 20L) {
  d <- bench[bench$K == k & bench$Mode == "pure_stability", ]
  ggplot(d, aes(reorder(Method, Intra_Uniqueness), Intra_Uniqueness, fill = Method)) +
    geom_col() + coord_flip() +
    facet_wrap(~ paste0(Tissue, " (", Platform, ")")) +
    labs(title = sprintf("Loading diversity (Intra_Uniqueness), K=%d", k),
         x = NULL, y = "Intra_Uniqueness") +
    theme_minimal() + theme(legend.position = "none")
}

# (b) Stability caveat: Mean_Stability per method x mode, by slide (K=20), with SD
fig_stability <- function(bench, k = 20L) {
  d <- bench[bench$K == k, ]
  ggplot(d, aes(Mode, Mean_Stability, fill = Method)) +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_errorbar(aes(ymin = Mean_Stability - SD_Stability,
                      ymax = Mean_Stability + SD_Stability),
                  position = position_dodge(0.8), width = 0.25) +
    facet_wrap(~ paste0(Tissue, " (", Platform, ")")) +
    labs(title = sprintf("Cross-run stability, K=%d", k), x = NULL, y = "Mean_Stability") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# (c) K-invariance: Intra_Uniqueness vs K per method (slides with both K)
fig_uniqueness_vs_k <- function(bench) {
  d <- bench[bench$Mode == "pure_stability", ]
  d <- d[d$Slide %in% names(which(tapply(d$K, d$Slide, function(x) length(unique(x)) > 1))), ]
  ggplot(d, aes(factor(K), Intra_Uniqueness, colour = Method, group = Method)) +
    geom_line() + geom_point() +
    facet_wrap(~ paste0(Tissue, " (", Platform, ")")) +
    labs(title = "Intra_Uniqueness vs number of components (K)",
         x = "K", y = "Intra_Uniqueness") +
    theme_minimal()
}

# (d) Underdetermination: consICA diversity gap vs spots/genes ratio (K=20, pure)
fig_gap_vs_ratio <- function(bench, k = 20L) {
  d <- bench[bench$K == k & bench$Mode == "pure_stability", ]
  gap <- do.call(rbind, by(d, d$Slide, function(s) {
    cc <- s$Intra_Uniqueness[s$Method == "consICA"]
    oth <- mean(s$Intra_Uniqueness[s$Method != "consICA"])
    data.frame(Slide = s$Slide[1], Tissue = s$Tissue[1],
               spots_per_gene = s$spots_per_gene[1], gap = cc - oth)
  }))
  ggplot(gap, aes(spots_per_gene, gap, label = Slide, colour = Tissue)) +
    geom_point(size = 3) + geom_text(vjust = -0.8, size = 3) +
    scale_x_log10() +
    labs(title = "consICA diversity advantage vs ICA determinedness",
         x = "spots / genes (log scale)", y = "Intra_Uniqueness gap (consICA - others)") +
    theme_minimal() + theme(legend.position = "none")
}

# 3. RUN ----

bench_all <- aggregate_benchmark(output_dir)
go_all <- aggregate_go(go_dir)

write.csv(bench_all, file.path(output_dir, "ALL_benchmark_metrics.csv"), row.names = FALSE)
if (!is.null(go_all)) write.csv(go_all, file.path(go_dir, "ALL_go_semantic.csv"), row.names = FALSE)

cat("Benchmark rows:", nrow(bench_all), "| slides:",
    paste(unique(bench_all$Slide), collapse = ", "), "\n")

# draw (comment out any you don't need)
print(fig_intra_uniqueness(bench_all))
print(fig_stability(bench_all))
print(fig_uniqueness_vs_k(bench_all))
print(fig_gap_vs_ratio(bench_all))
