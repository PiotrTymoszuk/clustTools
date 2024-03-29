% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_heat_maps.R
\name{plot_clust_hm}
\alias{plot_clust_hm}
\alias{plot_clust_hm.clust_analysis}
\alias{plot_clust_hm.min_analysis}
\alias{plot_clust_hm.combi_analysis}
\alias{plot_clust_hm.umatrix_analysis}
\title{Plot levels of clustering features in a heat map.}
\usage{
plot_clust_hm(x_object, ...)

\method{plot_clust_hm}{clust_analysis}(
  x_object,
  y_object = NULL,
  line_color = NA,
  plot_title = NULL,
  plot_subtitle = NULL,
  x_lab = "Sample",
  fill_lab = "Feature level",
  cust_theme = ggplot2::theme_classic(),
  discrete_fill = FALSE,
  ...
)

\method{plot_clust_hm}{min_analysis}(x_object, ...)

\method{plot_clust_hm}{combi_analysis}(
  x_object,
  y_object = NULL,
  line_color = NA,
  plot_title = NULL,
  plot_subtitle = NULL,
  x_lab = "Sample",
  fill_lab = "Feature level",
  cust_theme = ggplot2::theme_classic(),
  discrete_fill = FALSE,
  ...
)

\method{plot_clust_hm}{umatrix_analysis}(
  x_object,
  line_color = NA,
  plot_title = NULL,
  plot_subtitle = NULL,
  x_lab = "Sample",
  fill_lab = "Feature level",
  cust_theme = ggplot2::theme_classic(),
  discrete_fill = FALSE,
  ...
)
}
\arguments{
\item{x_object}{a \code{clust_analysis} or \code{combi_analysis} object, specifies
clustering of the observations.}

\item{...}{extra arguments passed to methods.}

\item{y_object}{a \code{clust_analysis} or \code{combi_analysis} object, specifies
clustering of the features, an optional parameter. Ignored in case of
multi-layer SOM.}

\item{line_color}{color of the line around heat map tiles.}

\item{plot_title}{plot title. If \code{NULL}, the plots generated for multi-layer
SOM analyses will be named after the data layers.}

\item{plot_subtitle}{plot subtitle.}

\item{x_lab}{x axis title.}

\item{fill_lab}{fill scale title.}

\item{cust_theme}{a ggplot theme.}

\item{discrete_fill}{logical, force a discrete fill scale?}
}
\value{
a \code{ggplot} object (single-layer analysis) or a list of
\code{ggplot} objects (multi-layer cases).
}
\description{
Generates a heat map of the clustering features, cluster
assignment is indicated by the plot faceting.
}
\details{
\code{plot_clust_hm()} is a S3 generic function.
Note that it is not possible to visualize clustering variable levels for
\code{clust_analysis} objects generated with user-provided dissimilarity matrices.
In such cases, \code{NULL} is returned with a warning.
}
