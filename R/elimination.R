# # full cluster analysis
# source("~/backup/ClusDec/clusdec/R/comb_new.R")
#
#
# fca <- function(clustered) {
#     clusters <- sort(unique(clustered[, 1]))
#     clusters <- clusters[clusters > 0]
#
#     results_all <- data.frame(Source=character(), Dest=character(), Score=numeric(), stringsAsFactors = F)
#
#     for (i in clusters) {
#         dest <- i
#         sources <- clusters[clusters != i]
#         x <- as.matrix(do.call(expand.grid, rep(list(c(FALSE, TRUE)), length(sources))))
#         x <- x[rowSums(x) > 0, ]
#         results <- apply(x, 1, function(j) {
#             source <- sources[j]
#             c(paste(source, collapse=" "), dest, checkCollinearFull(clustered, source, dest))
#         })
#         results_all <- rbind(results_all, t(results))
#     }
#     colnames(results_all) <- c("Source", "Dest", "Score")
#     print(as.numeric(results_all[, 3]))
#     results_all[, 3] <- as.numeric(as.character(results_all[, 3]))
#     return(results_all)
# }
#
#
# fca_max <- function(clustered, clusters) {
#     results_all <- lapply(clusters, function(i) {
#         dest <- i
#         sources <- clusters[clusters != i]
#         results <- checkCollinearFull(clustered, sources, dest)
#         list(source=paste(sources, collapse=" "), dest=dest, score=results)
#     })
#     # print(results_all)
#     results_all <- rbindlist(results_all)
#     results_all <- as.data.frame(results_all)
#     print(results_all)
#     colnames(results_all) <- c("Source", "Dest", "Score")
#     return(results_all)
# }
#
# elimination_process <- function(clustered, threshold=1, clusters=NULL) {
#     if (is.null(clusters)) {
#         clusters <- sort(unique(clustered[, 1]))
#         clusters <- clusters[clusters > 0]
#     }
#     while (T) {
#         results <- fca_max(clustered, clusters)
#         i <- which.max(results$Score)
#         elimination_candidate <- results[i, ]
#
#         if (elimination_candidate[3] >= threshold) {
#             print(results[i, ,drop=F], row.names=F)
#             cat("----")
#             clusters <- setdiff(clusters, elimination_candidate[2])
#         } else {
#             break
#         }
#     }
#     return(clusters)
# }
#
# # source("~/backup/ClusDec/stimulation/3d/3d_stimulation.R")
# # subset <- clustered[31:nrow(clustered), ]
# # fca(subset)
# #
# # clustered <- clusdec:::read.table.mine("~/reports/GSE27563_training/GSE27563_clustered.tsv")
# # clustered <- as.matrix(clustered)
# # all_combs <- fca_max(clustered)
# # all_combs <- all_combs[order(all_combs$Score, decreasing = T), ]
# # all_combs
# #
# #
# # clustered <- clusdec:::read.table.mine("~/backup/ClusDec/datasets/GSE52245_mcv4/mcv4_clustered.tsv")
# # clustered <- as.matrix(clustered)
# # basis <- elimination_process(clustered, 0.7)
# # cat("basis is:")
# # cat(basis)
