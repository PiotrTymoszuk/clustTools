% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cross_distance_utils.R
\name{cross_single_homolog}
\alias{cross_single_homolog}
\alias{cross_single_heterolog}
\alias{cross_multi_homolog}
\alias{cross_multi_heterolog}
\title{Cross-distances calculation helpers.}
\usage{
cross_single_homolog(x, method)

cross_single_heterolog(x, y, method)

cross_multi_homolog(x)

cross_multi_heterolog(x, y)
}
\arguments{
\item{x}{a \code{clust_analysis} or \code{combi_analysis} object.}

\item{method}{name of the distance metric. If NULL, it will be extracted
from the \code{x} object.}

\item{y}{a \code{clust_analysis} or \code{combi_analysis} object.}
}
\value{
an instance of the \code{cross_dist} class.
}
\description{
The functions compute homologous and heterologous distance between the
clusters of a single or two clustering analysis objects, respectively.
}
\details{
Designed solely for internal use. \code{cross_single_homolog()} and
\code{cross_single_heterolog()} compute cross-distances for single layer data.
\code{cross_multi_homolog()} and \code{cross_multi_heterolog()} tackle multi-layer SOM.
In multi-layer analyses, the cross-distance matrix is a weighted sum of
cross-distances between the corresponding data layers. The weights are
extracted from the SOM analysis object.
}
