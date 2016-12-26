# alpha-beta deconvolution
# library(nnls)

findMinumum <- function(clustered, clusters,
                        alphasRange=c(0, 1),
                        betasRange=c(0, 1),
                        maxIter=10) {
    mixedData <- clustered[, 2:ncol(clustered)]


    subsets <- lapply(clusters, function(i) {
        normSub <- t(apply(clustered[clustered[, 1] == i, 2:ncol(clustered)],
                           1, function(r) r / sqrt(sum(r^2))))
        colMeans(normSub)
    })
    gg <- do.call(rbind, subsets)
    ones <- matrix(1, nrow = ncol(gg), ncol=1)
    cc <- lsfit(t(gg), ones, intercept = F)$coefficients

    goodness <- function(alpha, beta) {
        coef1 = (beta - 1) / (alpha * beta - 1)
        coef2 = (alpha - 1) / (alpha * beta - 1)

        cc1 <- 1 / (cc / c(coef1, coef2))
        A <- cc1[1]
        B <- cc1[2]

        pseudoBasis <- matrix(c(A, beta * B, alpha * A, B), nrow=2, ncol=2)
        res <- Inf
        tryCatch({
            props <- fcnnls(pseudoBasis, gg)$x
            a <- t(props)
            b <- t(mixedData)

            basis_results <- t(fcnnls(a, b)$x)
            results <- basis_results %*% props
            res <- norm(log2(results + 1) - log2(mixedData + 1), "F")
        }, error=function(e) {
            print(alpha)
            print(beta)
            print(e)
        })
        return(res)

    }
    goodnessBeta <- function(beta, alpha) goodness(alpha, beta)

    x <- alphasRange[1]
    y <- betasRange[1]
    iter <- 0
    turn <- 1
    while (iter < maxIter) {
        print(paste0("Alpha is ", x, " ; beta is ", y))
        if (turn == 1) {
            x <- ternarySearch(goodness, alphasRange, beta=y)
        } else {
            y <- ternarySearch(goodnessBeta, alphasRange, alpha=x)
        }
        turn <- turn %% 2 + 1
        iter <- iter + 1
    }
    return(c(x, y))
}

abDeconvolution <- function(clustered, clusters, alpha, beta) {
    mixedData <- clustered[, 2:ncol(clustered)]
    subsets <- lapply(clusters, function(i) {
        normSub <- t(apply(clustered[clustered[, 1] == i, 2:ncol(clustered)],
                           1, function(r) r / sqrt(sum(r^2))))
        colMeans(normSub)
    })
    gg <- do.call(rbind, subsets)
    ones <- matrix(1, nrow = ncol(gg), ncol=1)
    cc <- lsfit(t(gg), ones, intercept = F)$coefficients

    coef1 = (beta - 1) / (alpha * beta - 1)
    coef2 = (alpha - 1) / (alpha * beta - 1)

    cc1 <- 1 / (cc / c(coef1, coef2))
    A <- cc1[1]
    B <- cc1[2]

    pseudoBasis <- matrix(c(A, beta * B, alpha * A, B), nrow=2, ncol=2)
    props <- fcnnls(pseudoBasis, gg)$x
    a <- t(props)
    a1 <- cbind(a, rep(1, nrow(a)))
    b <- t(mixedData)

    basis_results <- t(fcnnls(a, b)$x)
    print(a)
    basis_triple <- t(NMF:::fcnnls(a1, b, pseudo=T)$x)
    results <- basis_results %*% props
    return(list(W=basis_results, X=results, H=props, BT=basis_triple))
}
