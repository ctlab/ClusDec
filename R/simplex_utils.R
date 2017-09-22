## the whole fucking bunch of current utilities

#' Selection of points at the ends of simplex using VCA algorithms several times
#'
#' @param reduced matrix or data.frame object representing objects in reduced space
#' @param dims positive number, dimensionality of reduced space
#' @param number_of_iterations number of iterations to perform VCA
#'
#' @return data.frame contains endpoints of reduced dataset and associated cluster number
#' @export
#'
#' @examples
point_selection <- function(reduced, dims, number_of_iterations=100) {
    results <- matrix(nrow=0, ncol=ncol(reduced))
    colnames(results) <- colnames(reduced)

    center <- colMeans(reduced)

    for (i in 1:number_of_iterations) {
        endm <- vca05(reduced, p=dims)
        results <- rbind(results, reduced[endm, ])
        reduced <- reduced[-endm, ]
    }

    directions <- t(apply(results, 1, function(r) r - center))
    directions <- t(apply(directions, 1, function(r) r / sqrt(sum(r^2))))
    clusters <- Kmeans(directions, directions[1:dims, ], method="euclidean", iter.max=20000)$cluster

    return(cbind(results, clusters))
}

point_selection_vca <- point_selection

getEndpoints <- function(x, k) {
    x_mean <- rowMeans(x)
    em <- matrix(x_mean, ncol=1)

    get_next <- function(em) {
        props <- pseudoinverse(em) %*% x
        reconstruct <- em %*% props
        errors <- x - reconstruct
        rmse <- apply(errors, 2, function(cc) mean(cc^2))
        which.max(rmse)
    }

    endpoints <- get_next(em)

    while (length(endpoints) < k) {
        em <- x[, endpoints, drop=F]
        endpoints <- c(endpoints, get_next(em))
    }

    return(endpoints)
}

#' Selection of points at the ends of simplex using IEA algorithms several times
#'
#' @param reduced matrix or data.frame object representing objects in reduced space
#' @param dims positive number, dimensionality of reduced space
#' @param number_of_iterations number of iterations to perform VCA
#'
#' @return data.frame contains endpoints of reduced dataset and associated cluster number
#' @export
#'
#' @examples
point_selection_iea <- function(x, dims, number_of_iterations=100, toNormalise=T) {
    results <- matrix(nrow=0, ncol=ncol(x))
    colnames(results) <- colnames(x)

    if (toNormalise) {
        x_norm <- x / rowSums(x)
    } else {
        x_norm <- x
    }
    center <- colMeans(x_norm)
    x_norm_t <- t(x_norm)

    for (i in 1:number_of_iterations) {
        endm <- getEndpoints(x_norm_t, dims)
        results <- rbind(results, x_norm[names(endm), ])
        x_norm_t <- x_norm_t[, -endm]
    }

    directions <- t(apply(results, 1, function(r) r - center))
    directions <- t(apply(directions, 1, function(r) r / sqrt(sum(r^2))))
    clusters <- Kmeans(directions, directions[1:dims, ], method="euclidean", iter.max=20000)$cluster

    return(cbind(results, clusters))
}





plot_rotation <- function(reduced, dims, step=5, colors=NULL, fileName="rotation.gif") {
    if (is.null(colors)) {
        for (i in seq(0, 180, step)) {
            png(sprintf("tmpRotation%03d.png", i), width=4, height=4, units="in", res=300)
            s3d <- scatterplot3d(as.matrix(reduced[, (dims-2):dims]),
                                 angle = i, pch = 16, scale.y = 0.8, grid=F, cex.symbols = 0.5)
            dev.off()
        }
    } else {
        for (i in seq(0, 180, step)) {
            png(sprintf("tmpRotation%03d.png", i), width=4, height=4, units="in", res=300)
            s3d <- scatterplot3d(as.matrix(reduced[, (dims-2):dims]),
                                 color=colors,
                                 angle = i, pch = 16, scale.y = 0.8, grid=F, cex.symbols = 0.5)
            dev.off()
        }
    }

    system(paste0("convert -delay 20 tmpRotation*.png ", fileName))
    system("rm tmpRotation*")
}


#' Dimensionality reconstruction
#'
#' Returns a matrix of given rank which approximates (in terms of SVD) given matrix m
#'
#' @param m given matrix to be approximated
#' @param dims dimensionality of return matrix
#'
#' @return approximation of m with matrix of rank dims
#' @export
#'
#' @examples
dimensionality_reconstruction <- function(m, dims) {
    dec <- svd(m, nu = dims, nv = dims)
    approx <- matrix(0, nrow=nrow(m), ncol=ncol(m))
    for (i in 1:dims) {
        approx <- approx + dec$u[, i, drop=F] %*% dec$v[, i] * dec$d[i]
    }
    return(approx)
}

#' SVD based dimensionality reduction
#'
#' Projection of given amtrix m to space of rank dims
#'
#' If needed given points are projected to hyperplane in resulting space
#'
#' @param m given matrix to be projected
#' @param dims dimensionality of resulting matrix
#' @param projection project resulting matrix to dims-1 dimensional hyperplane
#'
#' @return
#' @export
#'
#' @examples
dimensionality_reduction <- function(m, dims, projection=FALSE) {
    m <- t(m)
    um <- svd(m, nu = dims, nv = 0)$u
    x <- crossprod(um, m)
    x <- t(x)
    if (projection) {
        u <- colMeans(x)
        x <- x / as.vector(x %*% u)
    }
    return(x)
}

#' Out of space error
#'
#' Calculates an error between given matrix M and its approximation matrix M_red of lower rank for every gene
#'
#' @param m given matrix M
#' @param m_red approximation of M of lower rank
#'
#' @return sorted vector of errors
#' @export
#'
#' @examples
dim_error <- function(m, m_red) {
    diff <- m - m_red
    diffs <- diff^2
    return(sort(rowSums(diffs), decreasing = T))
}







