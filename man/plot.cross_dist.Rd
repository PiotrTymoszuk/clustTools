% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cross_dist_oop.R
\name{plot.cross_dist}
\alias{plot.cross_dist}
\title{Plots of cross-distances.}
\usage{
\method{plot}{cross_dist}(
  x,
  type = c("heat_map", "mean", "histogram"),
  reorder = FALSE,
  upper = TRUE,
  signif_digits = 2,
  line_color = "black",
  show_txt = TRUE,
  txt_size = 2.75,
  labeller = NULL,
  cust_theme = ggplot2::theme_classic(),
  ...
)
}
\arguments{
\item{x}{a \code{cross_distance} class object.}

\item{type}{type of the plot:
\code{heat_map} (default) generates a heat map of
homologous or heterologous cross-distances for observation pairs with
mean distances and 95\% ranges of distances,
\code{mean} plots mean distances with 95\% ranges as a heat map, and
\code{histogram} generates a faceted panel of cross-distance histograms
(for heterologous distances, x object clusters are represented by horizontal
facets, y object clusters are represented by vertical facets) .}

\item{reorder}{logical: should distances in the heat maps be ordered
by mean distance? Defaults to FALSE.}

\item{upper}{should the upper half of the distance heat map be plotted?
Defaults to TRUE.}

\item{signif_digits}{significant digits for mean distances and distance
ranges presented in the heat map.}

\item{line_color}{color of the line around the tile, used only if \code{type}
is set to 'mean'.}

\item{show_txt}{logical, should the mean distance be presented in the plot?}

\item{txt_size}{of the mean distance text.}

\item{labeller}{a \code{\link[ggplot2]{labeller}} object to provide
customized labels of the facets of the histogram panel.}

\item{cust_theme}{a custom ggplot theme.}

\item{...}{extra arguments, such as color or number of bins, passed to
\code{\link[ggplot2]{geom_histogram}}.}
}
\value{
a \code{ggplot} graphic, whose elements like themes or fill scales can be easily
modified by the user.
}
\description{
Visualizes pairwise cross-distances as heat maps for observation pairs,
heat maps of average cross-distances between the clusters or histograms.
}
