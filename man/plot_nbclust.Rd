% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_utils.R
\name{plot_nbclust}
\alias{plot_nbclust}
\title{Plot WSS curve and silhouette statistic values as a function of cluster
number.}
\usage{
plot_nbclust(
  data,
  k,
  FUNcluster = NULL,
  method = c("wss", "silhouette", "gap_stat"),
  plot_title = NULL,
  plot_subtitle = NULL,
  plot_tag = NULL,
  cust_theme = ggplot2::theme_classic(),
  ...
)
}
\arguments{
\item{data}{a numeric matrix with the distances or a data frame.}

\item{k}{an integer, the cluster number.}

\item{FUNcluster}{a clustering function. See:
\code{\link[factoextra]{fviz_nbclust}} for details.}

\item{method}{a statistic to be plotted. See:
\code{\link[factoextra]{fviz_nbclust}} for details.}

\item{plot_title}{plot title.}

\item{plot_subtitle}{plot subtitle.}

\item{plot_tag}{plot tag.}

\item{cust_theme}{custom plot theme, a ggplot2 theme object.}

\item{...}{extra arguments passed to \code{\link[factoextra]{fviz_nbclust}}.}
}
\value{
a ggplot object.
}
\description{
Plots the values of the total within-cluster sum-of-squares and
silhouette statistic as a function of the cluster number.
}
\details{
Takes a distance matrix (e.g. the \code{\link{get_kernel_info}}
output) and a clustering function, for the details, see:
\code{\link[factoextra]{fviz_nbclust}}.
}
\references{
Kassambara A, Mundt F. factoextra: Extract and Visualize the Results
of Multivariate Data Analyses. (2020) Available at:
https://cran.r-project.org/web/packages/factoextra/index.html

Rousseeuw PJ. Silhouettes: A graphical aid to the interpretation and
validation of cluster analysis. J Comput Appl Math (1987) 20:53–65.
doi:10.1016/0377-0427(87)90125-7
}
