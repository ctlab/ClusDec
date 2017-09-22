## complicated analysis



#' getting simplex endpoints
#'
#' @param V
#'
#' @return
#' @export
#'
#' @examples
getDenseEndpoints <- function(P, V=NULL, ...) {
    # dimensionality checks
    stopifnot(ncol(P) > nrow(P))

    n <- nrow(P)


    lower <- apply(P, 1, function(x) min(x))
    upper <- apply(P, 1, function(x) max(x))

    lower <- rep(lower, n + 1) - 1e-10
    upper <- rep(upper, n + 1) + 1e-10


    scaling <- apply(P, 1, function(x) max(x) - min(x))
    scaling <- 10 / scaling
    scaling <- prod(scaling)

    getDensity <- function(V) {
        if (checkVolume(V)) {
            k <- sum(countInside(V, P))
            return(k / ((scaling * hyperVolume(V) + 10^n / n) ^ (1 / n)))
            # return(k / (scaling * hyperVolume(V) + 10^n))
        } else {
            return(0)
        }
    }

    score <- function(par) {
        V <- matrix(par, nrow=n)
        return(-getDensity(V))
    }

    if (is.null(V)) {
        opt <- GenSA(NULL, score, lower, upper, control=list(
            smooth=FALSE, verbose=T, ...
        ))
    } else {
        par <- as.numeric(V)
        opt <- GenSA(par, score, lower, upper, control=list(
            smooth=FALSE, verbose=T, ...
        ))
    }

    return(opt)
}
