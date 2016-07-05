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
        return(2^ge + 1)
    return(ge)
}

logDataset <- function(ge) {
    if (is_logscale(ge)) 
        return(ge)
    return(log2(ge + 1))
}

createReport <- function(dataset, accuracy, where = ".") {
    clustersDirName <- paste0(where, "/clusters")
    dir.create(clustersDirName)
    clusterNumber <- max(dataset[, 1])
    lapply(1:clusterNumber, function(i) {
        subset <- dataset[dataset[, 1] == i, ]
        write.table(rownames(subset), paste0(clustersDirName, "/cluster", i, ".txt"), 
            sep = "\n", col.names = FALSE, row.names = FALSE, quote = FALSE)
    })
    stats <- clusterStats(dataset)
    write.table.mine(stats, "stats.tsv")
    results <- chooseBest(dataset, accuracy)
    plotProportions(results)
    ggsave("results.png", height = 5, width = 5, units = "in")
    write.table.mine(basis(results), "basis.tsv")
    write.table.mine(coef(results), "coef.tsv")
}

is_logscale <- function(x) {
    qx <- quantile(as.numeric(x))
    if (qx[5] - qx[1] > 100 || qx[5] > 100) {
        return(FALSE)
    } else {
        return(TRUE)
    }
}
