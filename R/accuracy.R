#' Accuracy for all combination of clusters
#'
#' @param dataset matrix, clustered dataset
#' @param cellTypesNumber numeric, expected number of cell types
#' @param cores numeric, how many cores to use
#'
#' @return matrix, matrix of estimated accuracy (log Frobenius norm) for every combination
#'
#' @examples
#' data('datasetLiverBrainLung')
#' preprocessed <- preprocessedrocessDataset(datasetLiverBrainLung, k=5) # 5 clusters
#' accuracy <- clusdecAccuracy(preprocessed, 3) # assuming 3 cell types
#'
#' \dontrun{
#' accuracy <- clusdecAccuracy(preprocessed, 3, cores=2) # using 2 CPU cores
#' }
#'
#' @export
clusdecAccuracy <- function(dataset, cellTypesNumber, cores = 1) {
    clusterNumber <- max(dataset[, 1])
    datasetGE <- dataset[, 2:ncol(dataset)]
    clustering <- as.numeric(dataset[, 1])
    
    iteration <- function(cls) {
        evalClusters(datasetGE, clustering, cls, clusterNumber, cellTypesNumber)
    }
    
    cmb <- combn(1:clusterNumber, cellTypesNumber)
    lcmb <- lapply(seq_len(ncol(cmb)), function(i) c(cmb[, i]))
    
    results <- do.call(rbind, mclapply(lcmb, iteration, mc.cores = cores))
    colnames(results) <- c(paste0("Cluster ", 1:cellTypesNumber), "sumToOneError", 
        "LogFrobNorm")
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
    proportions <- dsaResults$H
    sumToOneError <- norm(as.matrix(apply(proportions, 2, function(x) (1 - sum(x)))), 
        "F")
    
    evaluated <- (dsaResults$W %*% dsaResults$H)[rownames(dataset), ]
    
    logFrobNorm <- logFrobNormAccuracy(dataset, evaluated)
    
    result <- c(cls, sumToOneError, logFrobNorm)
    message(paste0("Deconvolving by clusters  ", paste(cls, collapse = " "), ". Accuracy is ", 
        logFrobNorm))
    return(result)
}


#' Deconvolution by given set of clusters
#'
#' Runs DSA with  given clusters as putative signatures
#'
#' @param dataset matrix, given dataset to perform deconvolution
#' @param clusters numeric vector, clusters to use as putative signatures
#'
#' @return list with matrices W and H
#'
#' @examples
#' data('datasetLiverBrainLung')
#' preprocessed <- preprocessedrocessDataset(datasetLiverBrainLung, k=5) # 5 clusters
#' results <- deconvolveClusters(preprocessed, c(1, 2, 4)) # use clusters 1, 2 and 4 as putative signatures
#'
#' @export
deconvolveClusters <- function(dataset, clusters) {
    message(paste0("Deconvolving dataset using clusters: ", paste(clusters, collapse = " ")))
    datasetGE <- dataset[, 2:ncol(dataset)]
    clustering <- as.numeric(dataset[, 1])
    results <- runDSA(datasetGE, clustering, clusters)
    rownames(results$H) <- paste0("Cluster ", clusters)
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
#' @return list with matrices W and H
#'
#' @examples
#' data('datasetLiverBrainLung')
#' preprocessed <- preprocessedrocessDataset(datasetLiverBrainLung, k=5) # 5 clusters
#' accuracy <- clusdecAccuracy(preprocessed, 3) # assuming 3 cell types
#' results <- chooseBest(preprocessed, accuracy) # choose best combination of clusters as putative sigantures
#'
#'
#' @export
chooseBest <- function(dataset, accuracy) {
    logFrobNorm <- ncol(accuracy)
    cellTypes <- ncol(accuracy) - 2
    accuracy <- accuracy[order(accuracy[, logFrobNorm]), ]
    deconvolveClusters(dataset, accuracy[1, 1:cellTypes])
}

logFrobNormAccuracy <- function(data1, data2) {
    data1 <- logDataset(data1)
    data2 <- logDataset(data2)
    norm(data1 - data2, "F")
}
