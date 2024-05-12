# Provides S3 generics and the default methods for functions 'kidnapped'
# from base R and other packages, and `as_` functions used for type conversion.

# Variance --------

#' Object's variance
#'
#' @description Computes variance statistic specific for the given object.
#' @details The default `var()` method is a wrapper around
#' \code{\link[stats]{var}}.
#' @param x an object. For the default method a numeric vector, matrix
#' or a data frame.
#' @param ... extra arguments passed to methods, e.g. \code{\link[stats]{var}}.
#' @export

  var <- function(x, ...) UseMethod('var')

#' @rdname var
#' @export

  var.default <- function(x, ...) stats::var(x, ...)

# Distance --------

#' Distance between observations.
#'
#' @description
#' Computes the distance between observations in a matrix, data frame
#' or other compatible objects.
#' @details The default `dist()` method is a wrapper around
#' \code{\link[stats]{dist}}.
#' @param x an object. For the default method a numeric matrix, data frame
#' or `dist` object.
#' @param ... arguments for methods, e.g. passed to \code{\link[stats]{dist}}.
#' @export

  dist <- function(x, ...) UseMethod('dist')

#' @rdname dist
#' @export

  dist.default <- function(x, ...) stats::dist(x, ...)

# END ------
