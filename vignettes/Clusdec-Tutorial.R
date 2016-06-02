## ---- message=FALSE, warning=FALSE---------------------------------------
library(clusdec)
data("datasetLiverBrainLung")
mixedGed <- datasetLiverBrainLung[, 10:42]
clusteredGed <- preprocessDataset(mixedGed, k=5)
accuracy <- clusdecAccuracy(clusteredGed, 3, cores=1)

## ------------------------------------------------------------------------
results <- chooseBest(clusteredGed, accuracy)
plotProportions(results)

## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

