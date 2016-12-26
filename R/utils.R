# utils

write.table.mine <- function(data, fn, ...) write.table(data, fn, sep = "\t", quote = FALSE,
    col.names = NA)
read.table.mine <- function(fn, ...) read.table(fn, sep = "\t", header = 1, row.names = 1,
    ...)
norm.length <- function(r) r/sqrt(sum(r^2))
norm.length.matrix <- function(m) t(apply(m, 1, norm.length))
norm.relative <- function(r) (r - min(r))/(max(r) - min(r))
norm.relative.matrix <- function(m) t(apply(m, 1, norm.relative))
similarity.cosine <- function(x, y) crossprod(x, y)/sqrt(crossprod(x) * crossprod(y))
space.metric <- function(x, y) 1 - similarity.cosine(x, y)

clusterStats <- function(dataset) {
    data.norm <- cbind(dataset[, 1], norm.length.matrix(dataset[, 2:ncol(dataset)]))

    clusters.counts <- max(data.norm[, 1])
    clusters.means <- matrix(ncol = (ncol(data.norm) - 1), nrow = clusters.counts)

    for (i in 1:clusters.counts) {
        subset <- data.norm[data.norm[, 1] == i, 2:ncol(data.norm)]
        means <- colMeans(subset)
        clusters.means[i, ] <- means/sqrt(sum(means^2))
    }

    cor.pairwise <- lapply(1:clusters.counts, function(i) {
        subset <- as.matrix(data.norm[data.norm[, 1] == i, 2:ncol(data.norm)])
        sapply(1:10000, function(i) {
            pair <- sample(1:nrow(subset), 2)
            cor(subset[pair[1], ], subset[pair[2], ])
        })
    })

    cos.pairwise <- lapply(1:clusters.counts, function(i) {
        subset <- as.matrix(data.norm[data.norm[, 1] == i, 2:ncol(data.norm)])
        sapply(1:10000, function(i) {
            pair <- sample(1:nrow(subset), 2)
            similarity.cosine(subset[pair[1], ], subset[pair[2], ])
        })
    })


    cos.pairwise.means <- sapply(cos.pairwise, mean)
    cos.pairwise.max <- sapply(cos.pairwise, max)

    cor.pairwise.means <- sapply(cor.pairwise, mean)

    result <- data.frame(pairwise_cor = cor.pairwise.means, pairwise_cos_mean = cos.pairwise.means)
    rownames(result) <- paste0("cluster", 1:clusters.counts)
    result <- result[order(-result[, 1]), ]
    return(result)
}

linearizeDataset <- function(ge) {
    if (is_logscale(ge))
        return(2^ge - 1)
    return(ge)
}

logDataset <- function(ge) {
    if (is_logscale(ge))
        return(ge)
    return(log2(ge + 1))
}

createReport <- function(dataset, accuracy, where = ".") {
    requireNamespace("ggplot2")
    dir.create(where, showWarnings = F)
    clustersDirName <- paste0(where, "/clusters")
    dir.create(clustersDirName, showWarnings = F)
    clusterNumber <- max(dataset[, 1])
    lapply(1:clusterNumber, function(i) {
        subset <- dataset[dataset[, 1] == i, ]
        write.table(rownames(subset), paste0(clustersDirName, "/cluster", i, ".txt"),
            sep = "\n", col.names = FALSE, row.names = FALSE, quote = FALSE)
    })
    stats <- clusterStats(dataset)
    write.table.mine(stats, paste0(where, "/stats.tsv"))
    write.table.mine(dataset, paste0(where, "/dataset.tsv"))
    results <- chooseBest(dataset, accuracy)
    plotProportions(results$H)
    ggplot2::ggsave(paste0(where, "/results.png"), height = 5, width = 5, units = "in")

    ggplot2::ggplot(as.data.frame(accuracy), ggplot2::aes(LogFrobNorm)) +
        ggplot2::geom_histogram() + ggplot2::theme_bw() +
        ggplot2::geom_vline(xintercept=min(accuracy[, ncol(accuracy)]),
                   linetype="dashed", color="green")
    ggplot2::ggsave(paste0(where, "/histo.png"), height = 5, width = 5, units = "in")

    write.table.mine(basis(results), paste0(where, "/basis.tsv"))
    write.table.mine(coef(results), paste0(where, "/coef.tsv"))
}

is_logscale <- function(x) {
    qx <- quantile(as.numeric(x), na.rm = T)
    if (qx[5] - qx[1] > 100 || qx[5] > 100) {
        return(FALSE)
    } else {
        return(TRUE)
    }
}


r2Profiler <- function(r2Table,
                       rCheck=c(seq(0, 0.5, 0.1), seq(0.55, 0.9, 0.05), seq(0.92, 0.98, 0.02), 0.99),
                       kCheck=c(1, 3, 5, 10)) {
        # rTemp <- (r2Table + t(r2Table)) / 2
        results <- do.call(cbind, lapply(rCheck, function(rTreshold) {
            rMask <- (r2Table >= rTreshold)
            rMaskSums <- rowSums(rMask)
            sapply(kCheck, function(kn) {
                sum(rMaskSums > kn)
            })
        }))
        colnames(results) <- rCheck
        rownames(results) <- kCheck
        return(results)
}

visualizeProfilerResults <- function(results) {
    requireNamespace("ggplot2")
    requireNamespace("reshape2")
    melted <- reshape2::melt(t(results))
    colnames(melted) <- c("R2", "K", "candidates")

    pl <- ggplot2::ggplot(data=melted, ggplot2::aes(x=R2, y=log10(candidates),
                                  group=as.factor(K), color=as.factor(K))) +
        ggplot2::geom_point() + ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = log10(500), color="green") +
        ggplot2::geom_hline(yintercept = log10(1000), color="red") +
        ggplot2::theme_bw()
    pl
}

r2Filtering <- function(r2Table, r2val, k) {
    if (r2val > 1) stop("Value of R^2 to filter should be less or equal")
    if (r2val < 0) warning("Negative values of R^2 might be not meaningful")
    if (k <= 0) stop("Number of neighbours should be positive")
    # rTemp <- (r2Table + t(r2Table)) / 2
    rSums <- rowSums(r2Table >= r2val)
    filteredGenes <- rownames(r2Table)[rSums > k]
    subset <- r2Table[filteredGenes, filteredGenes]
    return(list(filteredGenes=filteredGenes,
                r2Filtered=subset))
}

provideClusterInfo <- function(clusters, dataset) {
    cluster <- rep(0, nrow(dataset))
    dataset <- cbind(cluster, dataset)
    for (i in 1:length(clusters)) {
        dataset[clusters[[i]], 1] = i
    }
    return(dataset)
}

