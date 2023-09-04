# Class-inheritance testing functions

#' Test class inheritance.
#'
#' @description
#' Tests if the object is an instance of the `red_analysis`, `clust_analysis`,
#' `combi_analysis`, `importance`, `cross_dist` or `sil_extra` class.
#'
#' @return a logical value.
#'
#' @param x an object.
#'
#' @export

  is_clust_analysis <- function(x) inherits(x, 'clust_analysis')

#' @rdname is_clust_analysis
#' @export

  is_combi_analysis <- function(x) inherits(x, 'combi_analysis')

#' @rdname is_clust_analysis
#' @export

  is_red_analysis <- function(x) inherits(x, 'red_analysis')

#' @rdname is_clust_analysis
#' @export

  is_importance <- function(x) inherits(x, 'importance')

#' @rdname is_clust_analysis
#' @export

  is_cross_dist <- function(x) inherits(x, 'cross_dist')

#' @rdname is_clust_analysis
#' @export

  is_sil_extra <- function(x) inherits(x, 'sil_extra')

# END ------
