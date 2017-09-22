

#' Generate mixed data with or without noise
#'
#' Generates mixed data with or without noise under assumption of linear model
#'
#' @param samples number of samples
#' @param genes number of genes
#' @param cellTypes number of cell types
#' @param bias if not null bias should be a number in between 0 and 1, one of cell types presented in mix will be more abundunt and one cell type will be less abundunt
#' @param pureGenes number of genes which are not noisy
#' @param noiseDeviation standart deviation of normally distributed noise value
#' @param removeAngles if TRUE there will be no signature genes in the mix (simplex corners will be removed)
#' @param cutCoeff corners are removed as points that lie outside of cellTypes-dimensional sphere with radius of cutCoef
#' @param removeBorders if TRUE where will be no genes close to simplex border
#' @param borderShift number in between 0 and 1, every value in basis is guaranteed to be at least borderShift
#'
#' @return
#' @export
#'
#' @examples
generateMixedData <- function(samples=40, genes=12000, cellTypes=3, bias = NULL,
                              pureGenes = 0, noiseDeviation = 0,
                              removeAngles = F, cutCoef = 0.85,
                              removeBorders = F, borderShift = 0.25) {
    proportions <- sampleFromSimplexUniformly(samples, cellTypes, 100000)
    colnames(proportions) <- paste0("Sample ", 1:samples)
    rownames(proportions) <- paste0("Cell type ", 1:cellTypes)

    if (!is.null(bias)) {
        bias <- proportions[1, ] * bias
        proportions[1, ] <- proportions[1, ] - bias
        proportions[cellTypes, ] <- proportions[cellTypes, ] + bias
    }

    basis <- sampleFromSimplexUniformly(genes, cellTypes, 100000)
    if (removeAngles) {
        while (!all(sqrt(colSums(basis^2)) < cutCoef)) {
            badIds <- sqrt(colSums(basis^2)) >= cutCoef
            basis[, badIds] <- sampleFromSimplexUniformly(sum(badIds), cellTypes, 100000)
        }
    }

    if (removeBorders) {
        shift <- (1 / borderShift - 1) / k
        basis <- basis + shift
        basis <- apply(basis, 2, function(x) x / sum(x))
    }


    rownames(basis) <- paste0("Cell type ", 1:cellTypes)
    colnames(basis) <- paste0("Gene ", 1:genes)
    basis <- t(basis)
    basis <- basis * runif(genes, min=200, max=10000)

    data <- basis %*% proportions

    if (noiseDeviation > 0) {
        noise <- matrix(rnorm(length(data), sd=noiseDeviation),
                        nrow = genes, ncol = samples)
        noised <- data + noise
        noised[noised < 0] <- 0

        if (pureGenes > 0) {
            pure <- sample(1:genes, pureGenes)
            noised[pure, ] <- data[pure, ]
        }

        data <- noised
    }

    return(list(
        data = data, proportions = proportions, basis = basis
    ))

}

#' Generation of points uniformly distributed on k-dimensional standard simplex
#'
#' @param n number of poitns
#' @param k dimensionality
#' @param M grid size
#'
#'
#' @return matrix where columns are points
#'
#' @import data.table
#'
#' @examples
sampleFromSimplexUniformly <- function(n, k=3, M=100000) {
    X <- matrix(0, nrow = k + 1, ncol=n)
    X[k + 1, ] <- M

    X[2:k, ] <- replicate(n, sample(1:(M-1), k - 1))
    X <- apply(X, 2, sort)
    Y <- (X - X[c(k + 1, 1:k), ])[2:(k + 1), ]
    return(Y / M)
}
