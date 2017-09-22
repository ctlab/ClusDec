#' Get gene stats
#'
#' Gets gene statistics corresponding to found simplex ndpoints
#'
#' @param dataset gene expression dataset
#' @param endpoints signature gene expression patterns
#'
#' @return data.frame with different statistics for genes in dataset
#' @export
#'
#' @examples
geneMarkerStats <- function(dataset, endpoints) {
    Y <- dataset / rowSums(dataset)
    mm <- colMeans(Y)
    endshifted <- endpoints - mm

    dds <- apply(Y, 1, function(r) {
        shifted <- endpoints - r

        sr <- r - mm
        cosprod <- crossprod(endshifted, sr)
        cosine <- cosprod / sqrt(sum(sr^2))
        cosine <- cosine / sqrt(colSums(endshifted^2))
        c(sqrt(colSums(shifted^2)), cosine)
    })

    return(t(dds))

}
