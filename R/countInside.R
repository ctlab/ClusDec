# check how mush points are inside


#' Check how many points are inside of given simplex
#'
#' Barycentric transformation can be applied to set of points
#' to check which of points are lie inside fiven simplex
#'
#' @param V matrix, given simplex, n \times (n + 1) matrix describing (n + 1) points in n-dimensional space
#' @param P matrix, given points to check, n \times N matrix describing N points in n-dimensional space
#'
#' @return logical vector wether point is inside of simplex
#' @export
#'
#' @examples
countInside <- function(V, P) {
    # dimensionality checks
    stopifnot(nrow(V) == nrow(P))
    stopifnot(nrow(V) + 1 == ncol(V))


    nd <- nrow(V)
    nd1 <- nd + 1
    baryTransform <- matrix(ncol=nd, nrow=nd)
    # print(baryTransform)
    for (i in 1:nd) {
        # print(V[, i] - V[, nd1])
        baryTransform[, i] = V[, i] - V[, nd1]
    }

    if (!(class(try(solve(baryTransform), silent=T)) == "matrix")) {
        return(0)
    }
    baryTransformInverse <- solve(baryTransform)
    # print(baryTransform)
    # print(baryTransformInverse)
    baryCoords <- baryTransformInverse %*% (P - V[, nd1])
    baryCoords <- rbind(baryCoords, 1 - colSums(baryCoords))
    # print(baryCoords)
    inside <- apply(baryCoords, 2, function(x) all(x >= 0))
    return(inside)
}

#' Compute the volume of the given simplex
#'
#' @param V matrix, given simplex, n \times (n + 1) matrix describing (n + 1) points in n-dimensional space
#'
#' @return float, hypervolume of given simplex
#' @export
#'
#' @examples
hyperVolume <- function(V) {
    stopifnot(nrow(V) + 1 == ncol(V))
    n <- nrow(V)
    dV <- matrix(ncol=n, nrow=n)
    dV <- (V - V[, 1])[, 2:(n+1)]
    return(abs(det(dV) / factorial(n)))
}

#' Check how many points are inside of given simplex
#'
#' Number of points inside of given simplex is devided by hypervolume of given simplex
#'
#' @param V matrix, given simplex, n \times (n + 1) matrix describing (n + 1) points in n-dimensional space
#' @param P matrix, given points to check, n \times N matrix describing N points in n-dimensional space
#'
#' @return float, how "dense" points are in given simplex
#' @export
#'
#' @examples
simplexDensity <- function(V, P) {
    if (checkVolume(V)) {
        k <- sum(countInside(V, P))
        return(k / hyperVolume(V))
    } else {
        return(0)
    }
}

checkVolume <- function(V) {
    n <- nrow(V)
    dV <- matrix(ncol=n, nrow=n)
    dV <- (V - V[, 1])[, 2:(n+1)]
    return(class(try(solve(dV), silent=T)) == "matrix")
}


#' Check how many points are inside of given simplex
#'
#' Number of points inside of given simplex is devided by hypervolume of given simplex
#'
#' @param V matrix, given simplex, n \times (n + 1) matrix describing (n + 1) points in n-dimensional space
#' @param P matrix, given points to check, n \times N matrix describing N points in n-dimensional space
#'
#' @return float, how "dense" points are in given simplex
#' @export
#'
#' @examples
simplexEDensity <- function(V, P, scale = 1) {
    if (checkVolume(V)) {
        k <- sum(countInside(V, P))
        return(k / exp(1 / (scale * hyperVolume(V))))
    } else {
        return(0)
    }
}

#' Root density
#'
#' Number of points inside of given simplex is devided by hypervolume of given simplex
#'
#' @param V matrix, given simplex, n \times (n + 1) matrix describing (n + 1) points in n-dimensional space
#' @param P matrix, given points to check, n \times N matrix describing N points in n-dimensional space
#'
#' @return float, how "dense" points are in given simplex
#' @export
#'
#' @examples
simplexRDensity <- function(V, P, scale = 1) {
    if (checkVolume(V)) {
        k <- sum(countInside(V, P))
        return(k / (scale * hyperVolume(V)) ^ (1 / nrow(V)))
    } else {
        return(0)
    }
}

# ## show case
#
# exampleEndpoints <- matrix(c(
#     c(10, 10), c(10, 14), c(15, 14)
# ), nrow=2)
#
# exampleEndpoints <- matrix(c(
#  c(10, 10), c(11, 12), c(15, 14)
# ), nrow=2)
#
# x_min <- min(exampleEndpoints[1, ])
# x_max <- max(exampleEndpoints[1, ])
# y_min <- min(exampleEndpoints[2, ])
# y_max <- max(exampleEndpoints[2, ])
#
# pointsToCheck <- matrix(c(
#     runif(12000, min = x_min, max = x_max),
#     runif(12000, min = y_min, max = y_max)
# ), ncol=2)
# pointsToCheck <- t(pointsToCheck)
#
# isInside <- countInside(exampleEndpoints, pointsToCheck)
# pointsToPlot <- as.data.frame(t(pointsToCheck))
# colnames(pointsToPlot)  <- c("x", "y")
# pointsToPlot$inside <- isInside
#
# library(ggplot2)
# toPlot <- as.data.frame(t(exampleEndpoints))
# colnames(toPlot) <- c("x", "y")
#
#
# ggplot(toPlot, aes(x=x, y=y)) + geom_point() + theme_bw() +
#     geom_point(data=pointsToPlot, aes(color=inside)) +
#     geom_polygon(color="grey", fill="grey", alpha=0.5)




