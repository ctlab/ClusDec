#' Pairwise linear fits between genes
#'
#' Performs pairwise linear fit gene1 ~ k * gene2
#'
#' @param X matrix, given gene expression dataset, rows supposed to be genes
#'
#' @return list containing two N x N matrices: R2 and slopes
#' @export
pairwiseLinearFit <- function(X) {
    results <- .Call('clusdec_pairwiseLinearFit', t(X), PACKAGE="clusdec")
    rownames(results$slopes) <- rownames(X)
    colnames(results$slopes) <- rownames(X)
    rownames(results$R2) <- rownames(X)
    colnames(results$R2) <- rownames(X)
    return(results)
}


#' Pairwise Deming regression fits between genes
#'
#' Performs pairwise linear fit gene1 ~ k * gene2
#'
#' @param X matrix, given gene expression dataset, rows supposed to be genes
#'
#' @return list containing two N x N matrices: R2 and slopes
#' @export
pairwiseDemingRegression <- function(X) {
    results <- .Call('clusdec_pairwiseDemingRegression', t(X), PACKAGE="clusdec")
    rownames(results$slopes) <- rownames(X)
    colnames(results$slopes) <- rownames(X)
    rownames(results$R2) <- rownames(X)
    colnames(results$R2) <- rownames(X)
    return(results)
}
