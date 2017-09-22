#' Title
#'
#' @param dds distances
#' @param ngenes number of genes
#'
#' @return
#' @export
#'
#' @examples
markersFromDist <- function(dds, ngenes=200) {
    wholeClustering <- apply(dds, 1, which.min)
    clusters <- lapply(1:ncol(dds), function(i) {
        head(rownames(dds)[order(dds[, i])], ngenes)
    })
    allGenes <- unlist(clusters)
    repeatedGenes <- unique(allGenes[duplicated(allGenes)])
    clusters <- lapply(clusters, function(cluster) cluster[!cluster %in% repeatedGenes])
    return(list(clusters=clusters, whole=wholeClustering))
}
