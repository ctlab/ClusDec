#' Accuracy for all combination of clusters
#'
#' @param dataset matrix, clustered dataset
#' @param cellTypesNumber numeric, expected number of cell types
#' @param cores numeric, how many cores to use
#'
#' @return matrix, matrix of estimated accuracy (log Frobenius norm) for every combination
#' @export
clusdecAccuracy <- function(dataset, cellTypesNumber, cores=1) {
    clusterNumber <- max(dataset[, 1])
    datasetGE <- dataset[, 2:ncol(dataset)]
    clustering <- as.numeric(dataset[, 1])

    iteration <- function(cls) {
        evalClusters(datasetGE, clustering, cls,  clusterNumber, cellTypesNumber)
    }

    cmb <- combn(1:clusterNumber, cellTypesNumber)
    lcmb <- lapply(seq_len(ncol(cmb)), function(i) c(cmb[, i]))

    results <- do.call(rbind, mclapply(lcmb, iteration, mc.cores=cores))
    colnames(results)<- c(paste0("Cluster ", 1:cellTypesNumber), "sumToOneError", "LogFrobNorm")
    results
}

#' Evaluate accuracy of combination
#'
#' Uses given combination as putative signatures for deconvolution
#' Runs DSA algorithm
#'
#' Then calculates Frobenius norm between logarithms of
#' observed gene expression matrix
#' and estimated gene expression matrix
#'
#' @param dataset matrix, given dataset to perform deconvolution
#' @param clustering numeric vector, cluster of every row in dataset
#' @param cls numeric vector, given cluster to evaluate accuracy
#' @param clusterNumber numeric, total number of clusters in dataset
#' @param cellTypesNumber numeric, expected number of cell types
#'
#' @return vector of cluster, sum-to-one error and log frobenius norm
evalClusters <- function(dataset, clustering, cls, clusterNumber, cellTypesNumber) {
    dsaResults <- runDSA(dataset, clustering, cls)
    proportions <- coef(dsaResults)
    sumToOneError <- norm(as.matrix(apply(proportions, 2, function(x) (1 - sum(x)))), "F")

    evaluated <- (basis(dsaResults) %*% coef(dsaResults))[rownames(dataset), ]

    logEval = log2(evaluated + 1)
    logDataset = log2(dataset + 1)

    logFrobNorm <- norm(logDataset - logEval, "F")

    result <- c(cls, sumToOneError, logFrobNorm)
    message(paste0("Deconvolving by clusters  ", paste(cls, collapse = " "), ". Accuracy is ", logFrobNorm))
    return(result)
}


#' Deconvolution by given set of clusters
#'
#' Runs DSA with  given clusters as putative signatures
#'
#' @param dataset matrix, given dataset to perform deconvolution
#' @param clusters numeric vector, clusters to use as putative signatures
#'
#' @return NMF object, deconvolution results
#' @import NMF
#' @import CellMix
#' @export
deconvolveClusters <- function(dataset, clusters) {
    message(paste0("Deconvolving dataset using clusters: ", paste(clusters, collapse = " ")))
    datasetGE <- dataset[, 2:ncol(dataset)]
    clustering <- as.numeric(dataset[, 1])
    results <- runDSA(datasetGE, clustering, clusters)
    rownames(coef(results)) <- paste0("Cluster ", clusters)
    results
}

#' Choose best combination of clusters and deconvolve
#'
#' Chooses best combination of clusters in accuracy table and then
#' runs DSA using this combination of clusters as putative signatures
#'
#' @param dataset matrix, given dataset to perform deconvolution
#' @param accuracy matrix, given accuracy table, result of clusdecAccuracy function
#'
#' @return NMF object, deconvolution results
#' @export
chooseBest <- function(dataset, accuracy) {
    logFrobNorm <- ncol(accuracy)
    cellTypes <- ncol(accuracy) - 2
    accuracy <- accuracy[order(accuracy[, logFrobNorm]), ]
    deconvolveClusters(dataset, accuracy[1, 1:cellTypes])
}
