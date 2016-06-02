#' Run DSA
#'
#' Runs DSA with provided clusters as putative signatures
#'
#' @param dataset gene expression matrix
#' @param clustering numeric vector, clustering of the rows
#' @param clusters numeric vector, which clusters use as putative signatures
#'
#' @import CellMix
#'
#' @return NMF object, deconvolution results
runDSA <- function(dataset, clustering, clusters) {
    genes <- lapply(clusters, function(i) rownames(dataset[clustering == i, ]))
    ged(dataset, MarkerList(genes), method="DSA")
}
