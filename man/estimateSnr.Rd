% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vca.R
\name{estimateSnr}
\alias{estimateSnr}
\title{Signal To Noise estimation}
\usage{
estimateSnr(R, rM, x)
}
\arguments{
\item{R}{matrix containing points in high dimensional space}

\item{rM}{vector of feature means}

\item{x}{projection of R (shifted to zero) to lower dimensional space produced by SVD}
}
\value{
numeric, signal to noise ratio
}
\description{
Signal To Noise estimation
}
