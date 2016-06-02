#' Cosine simillarity clustering
#'
#' Clusters given dataset using cosine simillarity as distance
#' Uses Kmeans from amap
#'
#' @param ge matrix, given gene expression dataset
#' @param k numeric, number of clusters
#'
#' @return matrix, clustered dataset
#' @import amap
clusterCosine <- function(ge, k) {
    # cosine simillarity between two vectors is proportional to
    # euclidean distance between two vectors normalized by their length
    geNorm <- t(apply(ge, 1, function(r) r / sqrt(sum(r^2))))
    clustering <- Kmeans(x=geNorm, centers=k, method="euclidean", iter.max=20000)
    clustered <- cbind(clustering$cluster, ge)
    colnames(clustered)[1] <- paste0("X", k, ".Clusters")
    rownames(clustered) <- rownames(ge)
    return(clustered)
}
