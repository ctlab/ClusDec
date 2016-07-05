## ---- message=FALSE, warning=FALSE---------------------------------------
library(clusdec)
set.seed(31)

## ---- message=FALSE, warning=FALSE---------------------------------------
data("datasetLiverBrainLung")
data("proportionsLiverBrainLung")

head(datasetLiverBrainLung[, 10:15])
head(proportionsLiverBrainLung[, 10:15])

## ---- message=FALSE, warning=FALSE---------------------------------------
mixedGed <- datasetLiverBrainLung[, 10:42]
mixedProportions <- proportionsLiverBrainLung[, 10:42]

## ---- message=FALSE, warning=FALSE---------------------------------------
clusteredGed <- preprocessDataset(mixedGed, k=5)
head(clusteredGed[, 1:6]) # every gene associated with cluster

## ---- message=FALSE, warning=FALSE---------------------------------------
accuracy <- clusdecAccuracy(clusteredGed, 3) 
head(accuracy)

## ---- message=FALSE, warning=FALSE---------------------------------------
results <- chooseBest(clusteredGed, accuracy)
plotProportions(results$H, mixedProportions[c(1, 3, 2), ],
                pnames=c("Estimated", "Actual"))

## ---- message=FALSE, warning=FALSE---------------------------------------
library(clusdec)
data("datasetLiverBrainLung")

mixedGed <- datasetLiverBrainLung[, 10:42]
head(mixedGed[, 1:5])
nrow(mixedGed)

set.seed(31)
clusteredGed <- preprocessDataset(mixedGed, k=4)
head(clusteredGed[, 1:5], 10)
nrow(clusteredGed)


## ------------------------------------------------------------------------
accuracy <- clusdecAccuracy(clusteredGed, 3, cores=2)
accuracy

