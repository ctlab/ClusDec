# #
# library(clusdec)
# data("datasetLiverBrainLung")
# prep <- preprocessDataset(datasetLiverBrainLung, topGenes=100, samples=10:42)
#
# start <- proc.time()
# res <- pairwiseDemingRegression(prep)
# print(proc.time() - start)
# rownames(res$ev) <- rownames(prep)
# colnames(res$ev) <- rownames(prep)
# #
# ev <- (res$ev + t(res$ev)) / 2
# edges <- ev
# edges <- (edges - 0.5) * 2
# edges[edges < 0] <- 0
#
#
# pcaStyle <- function(gene1, gene2, slope) {
#     gg <- cbind(gene1, gene2)
#     par(mfrow=c(1 ,2))
#     mmax <- max(c(gene1, gene2))
#     plot(gene1, gene2, xlim=c(0, mmax), ylim=c(0, mmax))
#     abline(0, slope, col="green")
#     p1 <- sqrt(1 / (slope^2 + 1))
#     p2 <- slope * p1
#     rotation <- matrix(c(p1, p2, -p2, p1), nrow=2, ncol=2)
#     tt <- gg %*% rotation
#     plot(tt)
#     abline(h=0)
#     abline(v=0)
#     tc <- cov(tt)
#     print(tc)
#     tc[1, 1] / sum(abs(tc[c(1, 2, 4)]))
# }
#
# pcaStyleNorm <- function(gene1, gene2, slope) {
#     gg <- cbind(gene1, gene2)
#     par(mfrow=c(1 ,2))
#     plot(gene1, gene2, xlim=c(0, max(gene1)), ylim=c(0, max(gene2)))
#     abline(0, slope, col="green")
#     p1 <- sqrt(1 / (slope^2 + 1))
#     p2 <- slope * p1
#     rotation <- matrix(c(p1, p2, -p2, p1), nrow=2, ncol=2)
#     tt <- gg %*% rotation
#     tt <- apply(tt, 2, function(cc){
#         if (max(cc) - min(cc) > 0) {
#             (cc - min(cc)) / (max(cc) - min(cc)) * 2 - 1
#         } else {
#             rep(0, length(cc))
#         }
#     })
#
#     plot(tt)
#     abline(h=0)
#     abline(v=0)
#
#     tc <- cov(tt)
#     print(tc)
#     tc[1, 1] / sum(abs(tc[c(1, 2, 4)]))
# }
#
# pcaStyle(gene1, gene2, 3.62)
# pcaStyle(gene2, gene1, 1/3.62)
# pcaStyle(mix[1, ], mix[2, ])
# norm_gene1 <- rnorm(2000) * 1000 + 10000
# norm_gene2 <- rnorm(2000) * 1000 + 10000
# pcaStyleNorm(norm_gene1, norm_gene2, 1)
# pcaStyleNorm(norm_gene2, norm_gene1, 1)
#
# liver_genes <- readLines("/media/askmebefore/1EC8FA3AC8FA102F/ClusDec/wgcna_tests/GSE19830_check/liver_genee_500")
# lg1 <- prep[liver_genes[2], ]
# lg2 <- prep[liver_genes[3], ]
# fit <- lm(lg2 ~ lg1 + 0)
# slope <- fit$coefficients
# plot(lg1, lg2, xlim=c(0, max(lg1)), ylim=c(0, max(lg2)))
# abline(fit)
#
# pcaStyle(lg1, lg2, slope)
# pcaStyleNorm(lg1, lg2, slope)
#
#
# genes <- c("COX1", "Fads1")
# lg1 <- prep[genes[1], ]
# lg2 <- prep[genes[2], ]
# fit <- lm(lg2 ~ lg1 + 0)
# slope <- fit$coefficients
# mmax <- max(c(lg1, lg2))
# plot(lg1, lg2, xlim=c(0, max(lg1)), ylim=c(0, max(lg2)))
# abline(fit)
#
# plot(lg1, lg2, xlim=c(0, max(lg1)), ylim=c(0, max(lg2)))
# abline(fit)
#
# pcaStyleNorm(lg1, lg2, slope)
#
# #
# # pheatmap(edges, show_rownames = F, show_colnames = F)
# #
# # clusters <- runSpici(edges)
# # dataToCheck <- provideClusterInfo(clusters, prep)
# # clusdec:::write.table.mine(dataToCheck, "~/toCheck.tsv")
# #
# # cppFunction('double p1(double slope) {
# #     return std::sqrt(1.0 / (slope * slope + 1.0));
# # }')
