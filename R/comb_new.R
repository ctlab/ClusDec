# new combinatorial approach
library(data.table)

demingRegression <- function(x, y) {
    n <- length(x)

    sxx <- crossprod(x) / (n - 1)
    syy <- crossprod(y) / (n - 1)
    sxy <- crossprod(x, y) / (n - 1)

    coef <- ((syy - sxx) + sqrt( (syy - sxx) * (syy - sxx) + 4 * sxy * sxy ) ) /
        (2 * sxy)

    return(coef)
}

linearityScore <- function(subset1, subset2) {
    subsetNorm1 <- t(apply(subset1, 1, function(r) r / max(r)))
    subsetNorm2 <- t(apply(subset2, 1, function(r) r / max(r)))
    s1 <- colMeans(subsetNorm1)
    s2 <- colMeans(subsetNorm2)

    # print(s1)
    # print(s2)

    coef <- demingRegression(s1, s2)

    s2_pred <- s1 * coef
    r2_1 <- 1 - sum(crossprod(s2 - s2_pred)) / sum(crossprod(s2 - mean(s2)))

    coef <- demingRegression(s2, s1)
    s1_pred <- s2 * coef
    r2_2 <- 1 - sum(crossprod(s1 - s1_pred)) / sum(crossprod(s1 - mean(s1)))

    # print(r2_1)
    # print(r2_2)
    return(mean(c(r2_1, r2_2)))
}


subsamplingScore <- function(subset1, subset2) {
    sampleSize <- 10
    samplingCount <- 20

    x <- sapply(1:samplingCount, function(i) {
        samp <- sample(1:nrow(subset1), sampleSize)
        norm(log2(subset1[samp, ] + 1) - log2(subset2[samp, ] + 1), "F")
    })

    return(mean(x))
}

checkColinear <- function(dataset, source, dest) {
    allClusters <- c(source, dest)
    subset <- as.matrix(dataset[dataset[, 1] %in% allClusters, ])
    clustering <- subset[, 1]
    subset <- subset[, 2:ncol(subset)]

    results <- clusdec:::runDSA(subset, clustering, source)

    subsetDest <- subset[clustering == dest, ]
    resultDest <- (results$W %*% results$H)[rownames(subsetDest), ]

    score_full <- norm(log2(subsetDest) - log2(resultDest), "F")
    score_sample <- subsamplingScore(subsetDest, resultDest)

    zSubsetDest <- t(apply(subsetDest, 1, zscore))
    zResultDest <- t(apply(resultDest, 1, zscore))

    cor_score <- cor(colMeans(zSubsetDest), colMeans(zResultDest))
    lin_score <- linearityScore(subsetDest, resultDest)

    return(c(score_full, score_sample, cor_score, lin_score))
}


checkCollinearFull <- function(dataset, source, dest, verbose=F) {
    ggs <- lapply(c(source, dest), function(i) {
        subset <- dataset[dataset[, 1] == i, 2:ncol(dataset)]
        subset <- t(apply(subset, 1, function(r) r / sqrt(sum(r^2))))
        return(colMeans(subset))
    })
    ggs <- do.call(rbind, ggs)
    ggs <- t(ggs)
    vars <- ggs[, 1:(ncol(ggs) - 1), drop=F]
    response <- ggs[, ncol(ggs), drop=F]
    fit <- lsfit(vars, response, intercept = F)$coef

    if (verbose) {
        print("FIT:")
        print(fit)

    }
    # print(fit)

    s1 <- response
    s2 <- vars %*% t(t(fit))

    # print(s1)
    # print(s2)
    coef <- demingRegression(s1, s2)

    # print(coef)

    s2_pred <- s1 * as.numeric(coef)
    r2_1 <- 1 - sum(crossprod(s2 - s2_pred)) / sum(crossprod(s2 - mean(s2)))

    coef <- demingRegression(s2, s1)
    s1_pred <- s2 * as.numeric(coef)
    r2_2 <- 1 - sum(crossprod(s1 - s1_pred)) / sum(crossprod(s1 - mean(s1)))
    return(mean(c(r2_1, r2_2)))

}

zscore <- function(x) (x - mean(x)) / sd(x)

allCombs <- function(set) {
    n <- length(set)
    args <- rep(list(0:1), n)
    do.call(expand.grid, args)

}

checkColinear <- function(dataset, source, dest) {
    allClusters <- c(source, dest)
    subset <- as.matrix(dataset[dataset[, 1] %in% allClusters, ])
    clustering <- subset[, 1]
    subset <- subset[, 2:ncol(subset)]

    results <- clusdec:::runDSA(subset, clustering, source)

    subsetDest <- subset[clustering == dest, ]
    resultDest <- (results$W %*% results$H)[rownames(subsetDest), ]

    score_full <- norm(log2(subsetDest + 1) - log2(resultDest + 1), "F")
    score_sample <- subsamplingScore(subsetDest, resultDest)

    zSubsetDest <- t(apply(subsetDest, 1, zscore))
    zResultDest <- t(apply(resultDest, 1, zscore))

    cor_score <- cor(colMeans(zSubsetDest), colMeans(zResultDest))
    lin_score <- linearityScore(subsetDest, resultDest)

    return(c(score_full, score_sample, cor_score, lin_score))
}

createBasis <- function(clustered, threshold=1) {
    clusters <- unique(clustered[, 1])
    clusters <- clusters[clusters > 0]
    # clusters <- sample(clusters)

    basis <- c(clusters[1])
    current <- 2
    while (current <= length(clusters)) {
        next_cluster <- clusters[current]
        message(" ")
        message("Current basis is ")
        message(paste(basis, collapse=" "))
        message("Trying to add ", next_cluster)


        n <- length(basis)
        x <- as.matrix(do.call(expand.grid, rep(list(c(FALSE, TRUE)), n)))
        x <- x[order(rowSums(x)), , drop=F]
        x <- x[2:nrow(x), , drop=F]

        scores <- apply(x,  1, function(r) {
            source <- basis[r]
            checkCollinearFull(clustered, source, next_cluster)
        })

        if (any(scores >= threshold)) {
            ind <- which(scores >= threshold)[1]
            source <- basis[x[ind, , drop=T]]
            message(
                sprintf(
                    "%s -> %s with score %s",
                    paste(source, collapse=" "),
                    paste(next_cluster),
                    paste(scores[ind])
                )
            )
            message("Can not add cluster ", next_cluster, " to basis")
            current <- current + 1
            next
        }

        message("Adding cluster ", next_cluster, " to basis")

        bad <- T
        while (bad) {
            bad <- F
            for (i in 1:length(basis)) {
                cl <- basis[i]
                x <- as.matrix(do.call(expand.grid, rep(list(c(FALSE, TRUE)), length(basis))))
                x <- x[order(rowSums(x)), , drop=F]
                x <- x[x[, i] == FALSE, , drop=F]

                scores <- apply(x,  1, function(r) {
                    source <- c(next_cluster, basis[r])
                    checkCollinearFull(clustered, source, cl)
                })

                if (any(scores >= threshold)) {
                    ind <- which(scores >= threshold)[1]
                    source <- c(next_cluster, basis[x[ind, , drop=T]])
                    message(
                        sprintf(
                            "%s -> %s with score %s",
                            paste(source, collapse=" "),
                            paste(cl),
                            paste(scores[ind])
                        )
                    )
                    message("removing cluster ", cl, " from basis")
                    basis <- basis[basis != cl]
                    bad <- T
                    break
                }

            }
        }

        basis <- c(basis, next_cluster)
        current <- current + 1
    }
    return(basis)
}

