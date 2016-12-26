

runSpici <- function(r2Table, k=10, spiciPath="spici", d=0.5, g=0.5) {
    requireNamespace("reshape2")
    infile <- tempfile()
    outfile <- tempfile()

    r2TableReshaped <- reshape2::melt(r2Table)
    filtered <-r2TableReshaped[r2TableReshaped[, 3] > 0, ]

    write.table(filtered, infile, sep="\t", col.names=F, row.names=F, quote=F)

    command <- paste0(spiciPath, " -s ", k, " -d ", d, " -g ", g, " -m 0 -i ", infile, " -o ", outfile)
    message(command)
    system(command)

    clusters <- readLines(outfile)
    clusters <- lapply(clusters, function(x) strsplit(x, "\t")[[1]])

    file.remove(infile)
    file.remove(outfile)

    return(clusters)
}
