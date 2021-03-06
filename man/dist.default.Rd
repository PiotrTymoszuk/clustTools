% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{dist.default}
\alias{dist.default}
\title{Distance Matrix Computation.}
\usage{
\method{dist}{default}(x, ...)
}
\arguments{
\item{x}{a numeric matrix, data frame or 'dist' object.}

\item{...}{extra arguments passed to \code{\link[stats]{dist}}.}
}
\description{
Default distance calculation method, a wrapper around
\code{\link[stats]{dist}}.
}
