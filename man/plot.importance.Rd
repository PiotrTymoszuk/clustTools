% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importance_oop.R
\name{plot.importance}
\alias{plot.importance}
\title{Plot feature importance as a scatter or bar plot.}
\usage{
\method{plot}{importance}(
  x,
  type = c("scatter", "bar"),
  fill_color = "cornsilk3",
  point_color = fill_color,
  point_size = 2,
  point_alpha = 0.5,
  point_wjitter = 0,
  point_hjitter = 0.1,
  box_alpha = 0.25,
  label = TRUE,
  txt_size = 2.75,
  signif_digits = 2,
  plot_title = NULL,
  plot_subtitle = NULL,
  plot_tag = NULL,
  cust_theme = ggplot2::theme_classic(),
  ...
)
}
\arguments{
\item{x}{an \code{importance} object.}

\item{type}{type of the plot: scatter or bar. Defaults to scatter.
This parameter is silently ignored, if evaluation of the importance was done
in multiple iterations
(e.g. \code{n_iter} set to > 1 in \code{\link{impact.clust_analysis}}).
In such cases, a box plot of importance metrics obtained in algorithm
iteration is generated.}

\item{fill_color}{fill color for the bars or boxes.}

\item{point_color}{size of the points, refers only to scatter and box plots.}

\item{point_size}{size of the points, refers only to scatter and box plots.}

\item{point_alpha}{alpha of the points, refers only to box plots.}

\item{point_wjitter}{width of the data point jittering, refers only to
box plots.}

\item{point_hjitter}{height of the data point jittering, refers only to
box plots.}

\item{box_alpha}{alpha of the boxes, refers only to box plots.}

\item{label}{logical, should the points be labeled with the importance
stat value?}

\item{txt_size}{label text size.}

\item{signif_digits}{significant digits for rounding of the statistic value.}

\item{plot_title}{plot title.}

\item{plot_subtitle}{plot subtitle.}

\item{plot_tag}{plot tag.}

\item{cust_theme}{a ggplot theme.}

\item{...}{extra arguments, currently none.}
}
\value{
a ggplot bar or scatter plot.
}
\description{
Generates a bar, scatter or box plot with the importance
statistic for the clustering variables.
The importance statistic is the difference in the
clustering variance fraction between the original clustering structure and
the clustering objects with the given variable reshuffled randomly.
}
