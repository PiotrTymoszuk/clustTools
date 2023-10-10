# Clustering variance

#' Calculate clustering variance.
#'
#' @description
#' Calculates the clustering sum of squares (total, within
#' clusters, total within clusters and between clusters) as well as the
#' fraction of 'explained' clustering variance. The later is the ratio of
#' the total between-cluster sum of squares to the total sum of squares.
#'
#' @details
#' `var()` is a S3 generic function.
#' `var()` overwrites the `var()` function provided by the stats package,
#' but provides a handy default method, so that `var()` is expected to behave
#' the same way as in base R.
#' For `combi_analysis` objects, variance is calculated based on distances
#' between the observations and the final cluster assignment -
#' nodes are ignored.
#'
#' @param x an object.
#' @param ... extra arguments, currently none.
#'
#' @return a list with the total, within-cluster, between-cluster sum of
#' squares and explained clustering variance (named `frac_var`).
#'
#' @export var.clust_analysis
#' @export

  var.clust_analysis <- function(x, ...) {

    stopifnot(is_clust_analysis(x))

    get_sum_sq(dist_mtx = x$dist_mtx,
               assignment = x$clust_assignment)

  }

#' @rdname var.clust_analysis
#' @export var.combi_analysis
#' @export

  var.combi_analysis <- function(x, ...) {

    stopifnot(is_combi_analysis(x))

    get_sum_sq(dist_mtx = as.matrix(dist(x, 'distance')[[1]]),
               assignment = x$clust_assignment)

  }

# END --------
