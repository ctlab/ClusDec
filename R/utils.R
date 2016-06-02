#utils

write.table.mine <- function(data, fn, ...) write.table(data, fn, sep="\t", quote=F, col.names=NA)
read.table.mine <- function(fn, ...) read.table(fn, sep="\t", header=1, row.names=1, ...)
norm.length <- function(r) r / sqrt(sum(r^2))
norm.length.matrix <- function(m) t(apply(m, 1, norm.length))
norm.relative <- function(r) (r - min(r)) / (max(r) - min(r))
norm.relative.matrix <- function(m) t(apply(m, 1, norm.relative))
similarity.cosine <- function(x, y) crossprod(x, y) / sqrt(crossprod(x) * crossprod(y))
space.metric <- function(x, y) 1 - similarity.cosine(x, y)

clusterStats <- function(dataset) {
    data.norm <- cbind(dataset[, 1], norm.length.matrix(dataset[, 2:ncol(dataset)]))

    clusters.counts <- max(data.norm[, 1])
    clusters.means <- matrix(ncol=(ncol(data.norm) - 1), nrow=clusters.counts)

    for (i in 1:clusters.counts) {
        subset <- data.norm[data.norm[, 1] == i, 2:ncol(data.norm)]
        means <- colMeans(subset)
        clusters.means[i, ] <- means / sqrt(sum(means^2))
    }

    dist.center <- lapply(1:clusters.counts, function(i) {
        subset <- data.norm[data.norm[, 1] == i, 2:ncol(data.norm)]
        apply(subset, 1, function(r) {
            space.metric(r, clusters.means[i, ])
        })
    })

    dist.center.means <- sapply(dist.center, mean)
    dist.center.max <- sapply(dist.center, max)


    dist.pairwise <- lapply(1:clusters.counts, function(i) {
        subset <- as.matrix(data.norm[data.norm[, 1] == i, 2:ncol(data.norm)])
        sapply(1:10000, function(i) {
            pair <- sample(1:nrow(subset), 2)
            space.metric(subset[pair[1], ], subset[pair[2], ])
        })
    })

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


    dist.pairwise.means <- sapply(dist.pairwise, mean)
    dist.pairwise.max <- sapply(dist.pairwise, max)
    cos.pairwise.means <- sapply(cos.pairwise, mean)
    cos.pairwise.max <- sapply(cos.pairwise, max)

    cor.pairwise.means <- sapply(cor.pairwise, mean)

    result <- data.frame(
        center_mean=dist.center.means,
        center_max=dist.center.max,
        pairwise_mean=dist.pairwise.means,
        pairwise_max=dist.pairwise.max,
        pairwise_cor=cor.pairwise.means,
        pairwise_cos_mean=cos.pairwise.means,
        pairwise_cos_max=cos.pairwise.max
    )
    rownames(result) <- paste0("cluster", 1:clusters.counts)
    result <- result[order(-result[, 5]), ]
    return(result)
}

linearizeDataset <- function(ge) {
    if (is_logscale(ge)) return(2^ge + 1)
    return(ge)
}

logDataset <- function(ge) {
    if (is_logscale(ge)) return(ge)
    return(log2(ge + 1))
}

createReport <- function() {
    # TODO
}
