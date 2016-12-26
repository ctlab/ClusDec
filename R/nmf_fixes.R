# .fccnls fixes
library(NMF)

.fcnnls <- function (x, y, verbose = FALSE, pseudo = FALSE, eps = 0, maxiter=100)
{
    if (any(dim(y) == 0L)) {
        stop("Empty target matrix 'y' [", paste(dim(y), collapse = " x "), "]")
    }
    if (any(dim(x) == 0L)) {
        stop("Empty regression variable matrix 'x' [",
             paste(dim(x), collapse = " x "), "]")
    }
    C <- x
    A <- y
    nObs = nrow(C)
    lVar = ncol(C)
    if (nrow(A) != nObs)
        stop("C and A have imcompatible sizes")
    pRHS = ncol(A)
    W = matrix(0, lVar, pRHS)
    iter = 0
    CtC = crossprod(C)
    CtA = crossprod(C, A)
    K = NMF:::.cssls(CtC, CtA, pseudo = pseudo)
    Pset = K > 0
    K[!Pset] = 0
    D = K
    Fset = which(colSums(Pset) != lVar)
    oitr = 0
    while (length(Fset) > 0) {
        oitr = oitr + 1
        if (verbose && oitr > 5)
            cat(sprintf("%d ", oitr))
        K[, Fset] = NMF:::.cssls(CtC, CtA[, Fset, drop = FALSE],
                           Pset[, Fset, drop = FALSE], pseudo = pseudo)
        Hset = Fset[colSums(K[, Fset, drop = FALSE] < eps) >
                        0]
        if (length(Hset) > 0) {
            nHset = length(Hset)
            alpha = matrix(0, lVar, nHset)
            while (nHset > 0 && (iter < maxiter)) {
                iter = iter + 1
                alpha[, 1:nHset] = Inf
                ij = which(Pset[, Hset, drop = FALSE] &
                               (K[, Hset, drop = FALSE] < eps), arr.ind = TRUE)
                i = ij[, 1]
                j = ij[, 2]
                if (length(i) == 0)
                    break
                hIdx = (j - 1) * lVar + i
                negIdx = (Hset[j] - 1) * lVar + i
                alpha[hIdx] = D[negIdx]/(D[negIdx] - K[negIdx])
                alpha.inf <- alpha[, 1:nHset, drop = FALSE]
                minIdx = max.col(-t(alpha.inf))
                alphaMin = alpha.inf[minIdx + (0:(nHset - 1) *
                                                   lVar)]
                alpha[, 1:nHset] = matrix(alphaMin, lVar, nHset,
                                          byrow = TRUE)
                D[, Hset] = D[, Hset, drop = FALSE] -
                    alpha[, 1:nHset, drop = FALSE] * (D[, Hset, drop = FALSE] - K[, Hset, drop = FALSE])
                idx2zero = (Hset - 1) * lVar + minIdx
                D[idx2zero] = 0
                Pset[idx2zero] = FALSE
                K[, Hset] = NMF:::.cssls(CtC, CtA[, Hset, drop = FALSE],
                                   Pset[, Hset, drop = FALSE], pseudo = pseudo)
                Hset = which(colSums(K < eps) > 0)
                nHset = length(Hset)
            }
        }
        W[, Fset] = CtA[, Fset, drop = FALSE] -
            CtC %*% K[, Fset, drop = FALSE]
        Jset = which(colSums((ifelse(!(Pset[, Fset, drop = FALSE]),
                                     1, 0) * W[, Fset, drop = FALSE]) > eps) == 0)
        Fset = setdiff(Fset, Fset[Jset])
        if (length(Fset) > 0) {
            mxidx = max.col(t(ifelse(!Pset[, Fset, drop = FALSE],
                                     1, 0) * W[, Fset, drop = FALSE]))
            Pset[(Fset - 1) * lVar + mxidx] = TRUE
            D[, Fset] = K[, Fset, drop = FALSE]
        }
    }
    list(coef = K, Pset = Pset)
}

assignInNamespace(".fcnnls", ".fcnnls", asNamespace("NMF"))
