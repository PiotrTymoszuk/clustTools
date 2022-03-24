# Provides S3 generics and the default methods.

# Distance ------

#' Distance Matrix Computation.
#'
#' @description Calculates or extracts a distance matrix.
#' @param x	a numeric matrix, data frame or 'dist' object.
#' @param ... extra arguments passed to methods.
#' @export

  dist <- function(x, ...) {

    UseMethod('dist', x)

  }

#' Distance Matrix Computation.
#'
#' @description Default distance calculation method, a wrapper around
#' \code{\link[stats]{dist}}.
#' @param x	a numeric matrix, data frame or 'dist' object.
#' @param ... extra arguments passed to \code{\link[stats]{dist}}.
#' @export dist.default
#' @export

  dist.default <- function(x, ...) {

    stats::dist(x, ...)

  }

# Object features ------

#' Extract object features.
#'
#' @description Extracts requested features from an object.
#' @param x an object.
#' @param ... extra arguments passed to methods.
#' @export

  extract <- function(x, ...) {

    UseMethod('extract', x)

  }

# Dimensionality reduction of an object ------

#' Reduce object data dimensions.
#'
#' @description Performs dimensionality reduction analysis with the given
#' object's data.
#' @param x an object.
#' @param ... extra arguments passed to methods.
#' @export

  components <- function(x, ...) {

    UseMethod('components', x)

  }

# Object variance ------

#' Calculate object's variance.
#'
#' @description Calculates variance of an object.
#' @param x an object.
#' @param ... extra arguments passed to methods.
#' @export

  var <- function(x, ...) {

    UseMethod('var', x)

  }

#' Calculate object's variance.
#'
#' @description Default variance calculation method, a wrapper around
#' \code{\link[stats]{var}}.
#' @param x an object.
#' @param ... extra arguments passed to \code{\link[stats]{var}}.
#' @export var.default
#' @export

  var.default <- function(x, ...) {

    stats::var(x, ...)

  }

# Observation and group N numbers -----

#' Get observation number.
#'
#' @description Returns the number of complete observations.
#' @param object an object.
#' @param ... extra arguments passed to methods.
#' @export

  nobs <- function(object, ...) {

    UseMethod('nobs', object)

  }

#' Get observation number.
#'
#' @description Returns the number of complete observations, a wrapper around
#' \code{\link[stats]{nobs}}.
#' @param object an object.
#' @param ... extra arguments passed to \code{\link[stats]{nobs}}.
#' @export nobs.default
#' @export

  nobs.default <- function(object, ...) {

    stats::nobs(object, ...)

  }

#' Get group assignment numbers.
#'
#' @description Returns the numbers of observations assigned to analysis groups
#' or clusters.
#' @param object an object.
#' @param ... extra arguments passed to methods.
#' @export

  ngroups <- function(x, ...) {

    UseMethod('ngroups', x)

  }

# Cross validation and variable importance -------

#' Cross-validate an object.
#'
#' @description Performs cross-validation of an object.
#' @param x an object.
#' @param ... extra arguments passed to methods.
#' @export

  cv <- function(x, ...) {

    UseMethod('cv', x)

  }

#' Get variable importance statistics.
#'
#' @description Computes variable importance statistics.
#' @param x an object.
#' @param ... extra arguments passed to methods.
#' @export

  impact <- function(x, ...) {

    UseMethod('impact', x)

  }

# END ------
