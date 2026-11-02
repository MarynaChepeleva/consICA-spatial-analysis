# Functions shared by the two demonstration scripts.
# bench_ica_perturbations.R keeps its own copy of auto_spot_threshold, with a
# different seed.

#' @description Auto per-spot QC threshold (min genes/spot): valley between
#' background and tissue in log10(genes/spot). If not clearly bimodal
#' (full-tissue slide, no background) -> returns the floor, no aggressive cut.
auto_spot_threshold <- function(counts, floor = 50L) {
  set.seed(2026)
  ng <- colSums(counts > 0)
  lg <- log10(ng[ng > 0] + 1)
  km <- stats::kmeans(lg, centers = 2, nstart = 10)
  if (abs(diff(km$centers)) < 0.5) {
    return(floor)
  }
  return(max(floor, round(10^max(lg[km$cluster == which.min(km$centers)]) - 1)))
}
