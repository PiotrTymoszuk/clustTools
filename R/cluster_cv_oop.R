# S3 methods for the cluster_cv class

# Summary ------

#' Cross-validation predictions and result summary.
#'
#' @description
#' The functions `extract.cluster_cv()` and `summary.cluster_cv()` retrieve
#' out-of-fold predictions, fold-means and global means with 95% confidence
#' intervals from results of cross-validation.
#'
#' @param x results of clustering cross-validation, e.g. generated with
#' \code{\link{cv}} for `clust_analysis` or `combi_analysis` objects.
#' @param object results of clustering cross-validation, e.g. generated with
#' \code{\link{cv}} for `clust_analysis` or `combi_analysis` objects.
#' @param type the element to be extracted: 'predictions' or 'fold_stats'.
#' @param ... extra arguments, currently none.
#'
#' @export extract.cluster_cv
#' @export

  extract.cluster_cv <- function(x,
                                 type = c('predictions', 'fold_stats'), ...) {

    stopifnot(is_cluster_cv(x))

    type <- match.arg(type[1], c('predictions', 'fold_stats'))

    x[[type]]

  }

#' @rdname extract.cluster_cv
#' @export

  summary.cluster_cv <- function(object, ...) {

    stopifnot(is_cluster_cv(object))

    object$summary

  }

# END ----
