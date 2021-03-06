#' Title
#'
#' @param R matrix describing points (possibly lyinh in a simplex) in high dimensional space
#' @param p number endpoints to find
#' @param SNR signal to noise ratio, NULL by default
#' @param verbose verbosity, deafult value is FALSE
#'
#' @return matrix of columns from R which are considered to be endpoints
#' @export
#'
#' @examples
vca <- function(R, p, SNR=NULL, verbose=F) {
    L <- nrow(R)
    N <- ncol(R)

    if (p < 0 || p > L) {
        stop("p is out of range (negative or too big)")
    }

    if (is.null(SNR)) {
        rm <- apply(R, 1, mean)
        R0 <- R - rm
        svdObj <- svd(R0 %*% t(R0), nu = p, nv = p)
        Ud <- svdObj$u
        xp <- t(Ud) %*% R0
        SNR <- estimateSnr(R, rm, xp)
        if (verbose) message(sprintf("SNR estimated = %f", SNR))
    } else {
        if (verbose) message(sprintf("input SNR = %f"))
    }

    SNRth <- 15 + 10 * log10(p)

    if (!is.na(SNR) && (SNR < SNRth)) {
        if (verbose) message("Select the projective projection")
        d <- p - 1

        if (exists("xp")) {
            Ud <- Ud[, 1:d]
        } else {
            rm <- apply(R, 1, mean)
            R0 <- R - rm
            svdObj <- svd(R0 %*% t(R0), nu = p, nv = p)
            Ud <- svdObj$u
            xp <- t(Ud) %*% R0
        }

        Rp <- Ud %*% xp[1:d, ] + rm
        x <- xp[1:d, ]
        c <- sqrt(max(colSums(x^2)))
        y <- rbind(x, rep(c, N))
    } else {
        if (verbose) message("Select prjection to p-1")
        d <- p
        Ud <- svd(R %*% t(R) / N, nu = d, nv = d)$u

        xp <- t(Ud) %*% R
        Rp <- Ud %*% xp[1:d, ]

        x <- xp
        u <- rowMeans(x)
        y <- x / matrix(kronecker(colSums(x * u), rep(1, d)), nrow=d)
    }

    ## VCA itself

    indice <- rep(0, p)
    A <- matrix(0, nrow=p, ncol=p)
    A[p, 1] <- 1

    for (i in 1:p) {
        w <- matrix(runif(p), ncol=1)
        f <- w - A %*% pseudoinverse(A) %*% w;
        f <- f / sqrt(sum(f^2))

        v <- t(f) %*% y
        indice[i] <- which.max(abs(v))
        A[, i] <- y[, indice[i]]
    }
    Ae = Rp[, indice]
    return(Ae)
}


#' Signal To Noise estimation
#'
#' @param R matrix containing points in high dimensional space
#' @param rM vector of feature means
#' @param x projection of R (shifted to zero) to lower dimensional space produced by SVD
#'
#' @return numeric, signal to noise ratio
#' @export
#'
#' @examples
estimateSnr <- function(R, rM, x) {
    L <- nrow(R)
    N <- ncol(R)
    p <- nrow(x)

    py <- sum(R^2) / N
    px <- sum(x^2) / N + crossprod(rM)

    return(10 * log10( (px - p * py / L) / (py - px) ))
}
