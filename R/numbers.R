# Numbers of observations and observations in the clusters.

# Number of observations and variables ------

#' Number of observations and variables for dimensionality reduction and clustering.
#'
#' @description
#' Computes numbers of observations and variables used in the analyses.
#'
#' @param object an object.
#' @param ... extra arguments, currently none.
#'
#' @return a list with the numbers of observations and variables.
#'
#' @export nobs.clust_analysis
#' @export

  nobs.clust_analysis <- function(object, ...) {

    stopifnot(is_clust_analysis(object))

    get_data_dim(eval_tidy(object$data))

  }

#' @rdname nobs.clust_analysis
#' @export nobs.red_analysis
#' @export

  nobs.red_analysis <- function(object, ...) {

    stopifnot(is_red_analysis(object))

    get_data_dim(eval_tidy(object$data))

  }

#' @rdname nobs.clust_analysis
#' @export nobs.combi_analysis
#' @export

  nobs.combi_analysis <- function(object, ...) {

    stopifnot(is_combi_analysis(object))

    map(object$clust_analyses, nobs)

  }

# Number of clusters ------

#' Numbers of observations in the clusters.
#'
#' @description
#' Compute numbers of observations in the clusters or, for `combi_analysis`
#' objects, numbers of observations in the SOM nodes and clusters.
#'
#' @details
#' `ngroups()` is a S3 generic function.
#'
#'
#' @param x an object.
#' @param ... extra arguments passed to methods, currently none.
#'
#' @export

  ngroups <- function(x, ...) {

    UseMethod('ngroups')

  }

#' @rdname ngroups
#' @export ngroups.clust_analysis
#' @export

  ngroups.clust_analysis <- function(x, ...) {

    stopifnot(is_clust_analysis(x))

    clust_id <- NULL

    count(x$clust_assignment, clust_id)

  }

#' @rdname ngroups
#' @export ngroups.combi_analysis
#' @export

  ngroups.combi_analysis <- function(x, ...) {

    stopifnot(is_combi_analysis(x))

    clust_id <- NULL

    c(map(x$clust_analyses, ngroups),
      list(final = count(x$clust_assignment, clust_id)))

  }

# END -------
