#' Draw a plot of estimated proportions
#'
#' Draws a plot of estimated proprotions
#' If ggplot2 and reshape2 are installed will use them and return ggplot object
#' Otherwise will use standart R functions
#'
#'@param proportions matricies, data frames or NMF objects of estimated proportions or paths to file
#'@examples
#'\dontrun{
#'# lets assume we already have clustered dataset and estimated accuracy
#'deconvolutionResults <- chooseBest(dataset, accuracy)
#'plotProportions(deconvolutionResults)
#'plotProportions(coef(deconvolutionResults)) # same as previous call
#'plotProportions("/path/to/deconvolution_results.tsv")
#'}
#'@export
plotProportions <- function(..., pnames=NULL) {
    requireNamespace(ggplot2)
    requireNamespace(reshape2)


    proportions <- list(...)
    proportions <- lapply(proportions, toMatrix)

    newRowNames <- do.call(function(...) {
        mapply(function(...) {
            dots <- list(...)
            rn <- unlist(dots)
            paste0(rn, collapse = "\n")
        }, ...)
    }, lapply(proportions, rownames))

    proportions <- lapply(proportions, function(p) {
        rownames(p) <- newRowNames
        p
    })

    names(proportions) <- pnames


    cellTypes <- nrow(proportions[[1]])
    results.m <- melt(proportions)
    results.m[, 4] <- as.factor(results.m[, 4])

    gplot <-ggplot(results.m, aes(x=as.numeric(Var2), y=value, fill=Var1, color=L1)) +
            geom_line() +
            geom_point() +
            scale_x_discrete(labels=colnames(proportions[[1]])) +
            facet_grid(Var1~.) + ylab("proportions") + ylim(0, 1.1) + theme_bw() +
            theme(axis.title.x = element_blank(),
                  axis.text.x  = element_text(angle=45, hjust=1)) + guides(fill=FALSE)
    if (length(proportions) > 1) {
        gplot <- gplot + theme(legend.title = element_blank(), legend.position="top")

    } else {
        gplot <- gplot + theme(legend.position="none")
    }
    gplot
}



toMatrix <- function(x) {
    if (class(x) == "data.frame") return(as.matrix(x))
    if (inherits(x, "NMF")) return(toMatrix(coef(x)))
    if (class(x) == "character") return(toMatrix(read.table.mine(x)))
    if (class(x) == "matrix") return(x)
    stop("invalid type for plotting")
}

