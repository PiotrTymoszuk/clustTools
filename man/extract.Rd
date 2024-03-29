% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extraction.R
\name{extract}
\alias{extract}
\alias{extract.clust_analysis}
\alias{model.frame.clust_analysis}
\alias{dist.clust_analysis}
\alias{dist.min_analysis}
\alias{extract.combi_analysis}
\alias{model.frame.combi_analysis}
\alias{dist.combi_analysis}
\alias{extract.red_analysis}
\alias{model.frame.red_analysis}
\alias{dist.red_analysis}
\title{Extract features of a clust_analysis object.}
\usage{
extract(x, ...)

\method{extract}{clust_analysis}(
  x,
  type = c("distance", "assignment", "clust_object", "data", "object", "umatrix"),
  ...
)

\method{model.frame}{clust_analysis}(formula, ...)

\method{dist}{clust_analysis}(x, type = c("distance", "umatrix"), ...)

\method{dist}{min_analysis}(x, ...)

\method{extract}{combi_analysis}(
  x,
  type = c("distance", "assignment", "clust_object", "data", "object", "umatrix"),
  ...
)

\method{model.frame}{combi_analysis}(formula, ...)

\method{dist}{combi_analysis}(x, type = c("distance", "umatrix"), ...)

\method{extract}{red_analysis}(
  x,
  type = c("component_tbl", "scores", "loadings", "data", "sdev", "object", "distance",
    "layout_distance"),
  ...
)

\method{model.frame}{red_analysis}(formula, ...)

\method{dist}{red_analysis}(x, type = c("distance", "layout"), ...)
}
\arguments{
\item{x}{an object.}

\item{...}{extra arguments, currently none.}

\item{type}{the feature name:
\itemize{
\item \code{distance} extracts the matrix with distances between the observations,
\item \code{data} the data set used for the analysis,
\item \code{assignment} assignment of the observations to the clusters,
\item \code{clust_object} or \code{object} returns the wrapped clustering object.
\item \code{clust_object} or \code{scores} return the component pr score tables for the
observations,
\item \code{loadings} retrieves the table of variable loadings,
\item \code{sdev} returns standard deviations, associated with the
components.
\item \code{umatrix} computes the U-matrix, i.e. weighted distance between
the self-organizing map (SOM) nodes. Available only for clustering analyses
done with SOM.
The U matrix is computed with \code{\link[kohonen]{object.distances}}.
\item \code{layout} and \code{layout_distance} computes distances between observations in
the dimensionality reduction layout.
}}

\item{formula}{an object.}
}
\value{
the requested feature/property.
}
\description{
A general extractor method for accessing properties and features
of a \code{clust_analysis}, \code{combi_analysis} and \code{red_analysis} object,
and specific methods for accessing the modeling data frame
and distance matrix.
}
\details{
\code{extract()} is a S3 generic function.
}
