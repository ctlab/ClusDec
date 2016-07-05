#' #' Draw a heatmap of clustered dataset
#' #'
#' #' Draws a heatmap of dataset clustered by preprocessDataset
#' #' If pheatmap is presented, heatmapClusters uses it
#' #' Othwerise stats::heatmap is used
#' #'
#' #'@param dataset matrix or data.frame of clustered dataset or path to file
#' #'@examples
#' #'\dontrun{
#' #'clusteredDataset <- preprocessDataset('GSE19830')
#' #'heatmapClusters(clusteredDataset)
#' #'heatmapClusters('/path/To/Dataset.tsv')
#' #'}
#' #'@import RColorBrewer
#' #'@export
#' setGeneric('heatmapClusters', function(dataset) {
#'     standardGeneric('heatmapClusters')
#' })
#' setMethod('heatmapClusters', c('data.frame'), function(dataset) heatmapClusters(as.matrix(dataset)))
#' setMethod('heatmapClusters', c('character'), function(dataset) heatmapClusters(read.table.mine(dataset)))
#' setMethod('heatmapClusters', c('matrix'), function(dataset) {
#'     if (requireNamespace('pheatmap')) {
#'         dataset_ordered <- dataset[order(dataset[, 1]), ]
#'         toDraw <- norm.length.matrix(dataset_ordered[, 2:ncol(dataset_ordered)])
#'         toDraw <- norm.relative.matrix(toDraw)
#'
#'         annotationRow <- as.data.frame(dataset[, 1, drop=FALSE])
#'         class(annotationRow[, 1]) <- 'character'
#'
#'         pheatmap(toDraw, annotation_row = annotationRow,
#'                  cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE,
#'                  color = colorRampPalette(c('blue', 'white', 'red'))(50))
#'     } else {
#'         dataset_ordered <- dataset[order(dataset[, 1]), ]
#'         toDraw <- norm.length.matrix(dataset_ordered[, 2:ncol(dataset_ordered)])
#'         toDraw <- norm.relative.matrix(toDraw)
#'         clusters <- max(dataset[, 1])
#'         colors <- sample(dscale(factor(1:clusters), hue_pal(l = 75)))[dataset_ordered[, 1]]
#'         heatmap(toDraw, Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(50), labRow='',
#'                 RowSideColors = colors)
#'     }
#' })
