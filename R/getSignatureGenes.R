## getting signature genes from dense endpoints

getSignatureGenes <- function(V, P, k=50) {
    n <- ncol(V)
    inside <- countInside(V, P)
    P <- P[, inside]
    sig <- list()
    for (i in 1:n) {
        diff <- P - V[, i]
        dists <- apply(diff, 2, function(x) sqrt(sum(x^2)))
        sig[[i]] <- names(sort(dists))[1:k]
    }
    return(sig)
}
