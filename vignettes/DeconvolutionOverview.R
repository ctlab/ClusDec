## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
H = matrix(c(1.0, 0.0, 0.5, 0.5, 0.0, 1.0), nrow=2, ncol=3)
colnames(H) <- c("Sample 1", "Sample 2", "Sample 3")
rownames(H) <- c("Brain", "Liver")
H

## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
library(clusdec)
data("datasetLiverBrainLung")
W = 2^datasetLiverBrainLung[1:6, c(4, 1)]
colnames(W) <- c("Brain", "Liver")
W

## ---- message=FALSE, warning=FALSE---------------------------------------
X <- W %*% H
X

