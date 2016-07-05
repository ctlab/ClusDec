![Travis-CI Build Status](https://api.travis-ci.org/ctlab/ClusDec.svg?branch=master)](https://travis-ci.org/ctlab/ClusDec)
[![codecov](https://codecov.io/gh/ctlab/ClusDec/branch/master/graph/badge.svg)](https://codecov.io/gh/ctlab/ClusDec)


# ClusDec
An R-package implementing de novo deconvolution method that allows identification of individual cell types in the mixture without knowing either cell types proportions or their corresponding cell-specific markers.

The basic idea of ClusDec is to find clusters of genes with linear expression profiles that then will be used as putative signatures for Digital Sortiing Algorithm (DSA) described in [(Zhong et al. BMC Bioinformatics 2013, 14:89)][http://dx.doi.org/10.1186/1471-2105-14-89].

The basic workflow consists several steps: preprocessing (including clustering), evaluating accuracies of combination of clusters and then using best combination as putative signatures for DSA algorithm (link).


## Installation

```{r}
library(devtools)
install_github("ctlab/ClusDec")
```

## Quick run

Loading library

```{r}
library(clusdec)
```

Loading example data:
```{r}
data("datasetLiverBrainLung")
data("proportionsLiverBrainLung")
```

We don't take first 9 samples: these are pure samples, and we only leave mixed samples.
```{r, message=FALSE, warning=FALSE}
mixedGed <- datasetLiverBrainLung[, 10:42]
mixedProportions <- proportionsLiverBrainLung[, 10:42]

head(mixedGed[, 1:6])
head(mixedProportions[, 1:6])
```

Data looks like
```
##         GSM495218 GSM495219 GSM495220 GSM495221 GSM495222 GSM495223
## A1bg    10.267250 10.374572 10.322145 13.064351 13.016001 13.018780
## A1cf     4.059308  4.458606  4.359954  7.676283  7.696134  7.653000
## A2bp1    9.490740  9.566670  9.422953  7.792056  7.831022  7.892111
## A2ld1    7.655956  7.912126  7.948867  8.187795  8.454980  8.261322
## A2m      6.127628  6.046218  6.080013  6.539769  6.443252  6.373277
## A3galt2  5.126264  5.325947  5.224698  4.968571  4.663900  4.890523

##       GSM495218 GSM495219 GSM495220 GSM495221 GSM495222 GSM495223
## Liver      0.05      0.05      0.05      0.70      0.70      0.70
## Brain      0.25      0.25      0.25      0.05      0.05      0.05
## Lung       0.70      0.70      0.70      0.25      0.25      0.25
```


Now lets run ClusDec and perform deconvolution:
```{r}
set.seed(31)
lusteredGed <- preprocessDataset(mixedGed, k=5)
accuracy <- clusdecAccuracy(clusteredGed, 3)
results <- chooseBest(clusteredGed, accuracy)
head(results$H[, 1:6])
```

```
          GSM495218 GSM495219 GSM495220 GSM495221 GSM495222 GSM495223
Cluster 1 0.1407633 0.1417054 0.1426791 0.5317571 0.5254577 0.5288620
Cluster 2 0.4741810 0.4675254 0.4708257 0.2650922 0.2688183 0.2654242
Cluster 3 0.3811201 0.3759724 0.3755612 0.1922302 0.1916823 0.1883313
```
ClusDec chooses first three clusters as putative signatures 
for performing deconvolution with DSA. Now we can compare these 
estimated results with acutal proporions.

```{r}
plotProportions(results$H, mixedProportions[c(1, 3, 2), ],
                pnames=c("Estimated", "Actual"))
```

![proportions.png](https://dl.dropboxusercontent.com/u/38245921/ClusDec/proportions.png)

