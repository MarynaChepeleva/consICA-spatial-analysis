# ICASpatial of CellsFromSpace, by Corentin Thuilliez, copied from
# github.com/gustaveroussy/CFS at commit 3926422e7b5264e92a75d665aa9ef447be88c273.
#
# This file is covered by that package's GPL, version 2 or any later version,
# see LICENSE.GPL-2. It is not covered by the MIT licence of the rest of this
# repository.
#
# The seed is the only change: CFS fixes it at 42, here it is a parameter the
# benchmark sets. It makes no difference to the result, since icafast is
# initialized from a fixed rotation matrix (Rmat = diag(nc)) and draws no
# random numbers. correct_sign is called from the package itself.

custom_ICASpatial <- function(data = NULL, nics = 100, maxit = 600,
                              method = "icafast", kurtosis = 3,
                              sd = 3, seed.use = 42, ...) {
  set.seed(seed.use)
  object <- data@assays$SCT@scale.data
  nics <- min(nics, ncol(x = data))
  ica.results <- ica::icafast(t(x = as.matrix(object)), nc = nics, maxit = maxit)

  cell.embeddings <- ica.results$S
  feature.loadings <- as.matrix(object) %*% as.matrix(cell.embeddings)
  colnames(feature.loadings) <- paste0("IC_", seq_len(ncol(feature.loadings)))
  colnames(cell.embeddings) <- paste0("IC_", seq_len(ncol(cell.embeddings)))
  rownames(cell.embeddings) <- Cells(data)

  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = "SCT", key = "IC_")
  data@reductions$ica <- reduction.data
  data <- CFS:::correct_sign(data)

  kurt_cob <- apply(data@reductions$ica@feature.loadings, 2,
                         function(x) e1071::kurtosis(x))
  Contrib_logic <- apply(data@reductions$ica@feature.loadings, 2,
                         function(x) abs((x - mean(x)) / sd(x)) > sd)
  Contrib_gene <- purrr::map2(
    as.data.frame(Contrib_logic),
    colnames(as.data.frame(Contrib_logic)),
    function(x, y) data.frame(
      gene = rownames(data@reductions$ica@feature.loadings)[x],
      Sig = as.data.frame(data@reductions$ica@feature.loadings)[x, y]))
  list_all_gene <- purrr::map(
    Contrib_gene[names(which(kurt_cob > kurtosis))],
    ~.x$gene) %>% unlist() %>% unique()
  data@misc$GeneAndStat <- list(
    Kurtosis_ICs = kurt_cob,
    Contrib_gene = Contrib_gene,
    All_Contrib_Gene = list_all_gene,
    kurtosis_value = kurtosis)
  return(data)
}
