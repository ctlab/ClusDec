# new preprocessing

#' Preprocess Dataset
#'
#' Preprocesses given dataset. Preprocessing consists of 3 major steps:
#' 1) If needed, probes corresponding to the same genes are collapsed, only most expressed probe is taken for further analysis.
#'    It's common technique in microarray data analysis.
#' 2) If needed, only highly expressed genes are taken for further analysis. (Say hello to noize reduction)
#' 3) All genes are clustered with Kmeans using cosine simillarity as distance.
#'
#' @param dataset matrix, data.frame, path to file or GSE accession with expression data
#' @param annotaion dataframe, matrix, named vector with annotation to probes
#' @param k number of clusters to perform Kmeans, default value is 10
#' @param geneSymbol column from annotation to collapse the genes, deafult value is "Gene Symbol"
#' @param samples character vector of samples. If column were not in samples, it would be excluded from analysis.
#' Default value is NULL, which takes every sample from dataset
#' @param topGenes integer How many genes include in analysis. We suppose to include only expressed genes. Default value is 10000
#' @param ...
#'
#' @return clustered dataset, matrix, first column identifies cluster of the row
#' @import CellMix
#' @export
setGeneric("preprocessDataset", function(dataset, annotaion, ...) {
    standardGeneric("preprocessDataset")
})
setMethod("preprocessDataset", c("matrix", "missing"),
          function(dataset, annotation, k=10, geneSymbol=NULL, samples=NULL, topGenes=10000) {

              if (!is.null(samples)) {
                  dataset <- dataset[, samples]
              }

              # finding top genes in log scale
              dataset <- logDataset(dataset)
              topGenes <- min(topGenes, nrow(dataset))
              topRows <- order(rowSums(dataset), decreasing = TRUE)[1:topGenes]
              topDataset <- dataset[topRows, ]
              # clustering in linear space
              topDataset <- linearizeDataset(topDataset)
              clustered <- clusterCosine(topDataset, k)
              return(clustered)
          }
)
setMethod("preprocessDataset", c("data.frame", "missing"),
          function(dataset, annotation, ...) {
              dataset <- as.matrix(dataset)
              preprocessDataset(dataset, ...)
          }
)
setMethod("preprocessDataset", c("character", "missing"),
          function(dataset, annotation, geneSymbol="Gene Symbol", samples=NULL, ...) {
              if (file.exists(dataset)) {
                  message(paste0("Reading dataset from file ", dataset))
                  dataset <- read.table.mine(dataset)
                  preprocessDataset(dataset, geneSymbol="Gene Symbol", samples=NULL, ...)
              } else {
                  message(paste0("Trying to get dataset from public GEO dataset ", dataset))
                  gse <- getGSE(dataset, verbose=TRUE)

                  dataset <- exprs(gse)

                  if (!is.null(samples)) {
                      dataset <- dataset[, samples]
                  }

                  fdata <- fData(gse)[, geneSymbol, drop=FALSE]
                  dataset <- collapseGenes(dataset, fdata)

                  preprocessDataset(dataset, geneSymbol="Gene Symbol", samples=NULL, ...)
              }

          }
)
