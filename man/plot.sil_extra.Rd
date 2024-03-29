% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sil_oop.R
\name{plot.sil_extra}
\alias{plot.sil_extra}
\title{Plots of silhouette statistics.}
\usage{
\method{plot}{sil_extra}(
  x,
  fill_by = c("cluster", "neighbor", "value", "sign"),
  show_stats = TRUE,
  signif_digits = 2,
  cust_theme = ggplot2::theme_classic(),
  bar_color = "black",
  ...
)
}
\arguments{
\item{x}{an object of the \code{\link{sil_extra}} class.}

\item{fill_by}{defines the color coding of the bar fill color.
For \code{cluster}, the bars are colored after cluster assignment of the
observations (default).
For \code{neighbor}, the bar color codes for the nearest neighbor cluster.
For \code{value}, the bar color codes for the silhouette width.
For \code{sign}, the bar color represents the sign of the silhouette width.}

\item{show_stats}{logical, should the number of observations in the cluster,
percentage of negative silhouette widths and average silhouette statistic
be shown in the plot? Defaults to TRUE.}

\item{signif_digits}{significant digits used for rounding of the statistics
presented in the plot.}

\item{cust_theme}{custom ggplot theme.}

\item{bar_color}{color of the bar line.}

\item{...}{extra arguments passed to \code{\link[ggplot2]{geom_bar}}.}
}
\value{
a \code{ggplot} class graphic.
}
\description{
Generates a classical bar plot of silhouette width distribution in clusters.
}
\references{
Rousseeuw PJ. Silhouettes: A graphical aid to the interpretation and
validation of cluster analysis. J Comput Appl Math (1987) 20:53–65.
doi:10.1016/0377-0427(87)90125-7
}
