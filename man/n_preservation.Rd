% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neighborhood_utils.R
\name{n_preservation}
\alias{n_preservation}
\alias{np_reduction}
\title{Neighborhood preservation.}
\usage{
n_preservation(
  data_dist,
  clust_assignment,
  clust_dist = NULL,
  kNN_data = 5,
  kNN_cluster = 1
)

np_reduction(data_dist, layout_dist, kNN_data = 5)
}
\arguments{
\item{data_dist}{a matrix of distances between data points.}

\item{clust_assignment}{a data frame with the cluster assignment scheme.}

\item{clust_dist}{an optional matrix of distances between the clusters.
Ignored if \code{kNN_cluster = 1}.}

\item{kNN_data}{number of k-nearest neighbors for data points.}

\item{kNN_cluster}{number of k-nearest neighbors for clusters}

\item{layout_dist}{a matrix of distances between data points following
dimensionality reduction.}
}
\value{
An object of the \code{\link{knb}} class with
\code{\link{summary.knb}} and
\code{\link{plot.knb}} methods.
}
\description{
Computes neighborhood preservation stats.
}
\details{
For internal use. For computation details, see: \code{\link{np}}.
}
