
data <- read.table.mine("tests/data/tmps.tsv")
samples <- colnames(data)
samples <- samples[! samples %in% c("GSM1841072", "GSM1841074")]

dataCl <- preprocessDataset(data, samples=samples)
png("tests/pheatmap1.png", height = 10, width=5, units="in", res=300)
heatmapClusters(dataCl)
dev.off()

dataCl <- preprocessDataset("tests/data/tmps.tsv", samples=samples)
png("tests/pheatmap2.png", height = 10, width=5, units="in", res=300)
heatmapClusters(dataCl)
dev.off()

gseLuBrLi <- preprocessDataset("GSE19830")
png("tests/pheatmap3.png", height = 10, width=5, units="in", res=300)
heatmapClusters(gseLuBrLi)
dev.off()

samples <- as.character(read.table("../samples/GSE19830_samples.txt", sep="\n")$V1)
check <- preprocessDataset("GSE19830", samples=samples)
accuracyTable <- clusdec_combinations_accuracy(check, 2, cores=2)
accuracyTable <- accuracyTable[order(accuracyTable[, 5]), ]
deconvolution <- chooseBest(check, accuracyTable)
write.table.mine(coef(deconvolution), "tests/results.tsv")
plotProportions("tests/results.tsv")
