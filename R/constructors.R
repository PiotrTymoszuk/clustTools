# Contains the class constructor functions.

# Dimensionality reduction analysis class ------

#' Construct a red_analysis class object.
#'
#' @description Constructs a red_analysis class object given a list storing
#' results of a dimensionality reduction analysis.
#' @details A named list with the following elements ir required:
#' red_obj (analysis output),
#' red_fun (name of the reduction function),
#' component_tbl (component values for the observations),
#' loadings (variable loadings, relevant only for PCA),
#' data (a quosure calling the original data set).
#' @param x a named list, see details.
#' @return a `red_analysis` object.
#' @export

  red_analysis <- function(x) {

    ## entry control

    if(!is.list(x)) {

      stop('Only named lists can be converted to a clust_analysis objects',
           call. = FALSE)

    }

    if(is.null(names(x))) {

      stop('Only named lists can be converted to a clust_analysis objects',
           call. = FALSE)

    }

    if(!all(c('red_obj',
              'red_fun',
              'component_tbl',
              'loadings',
              'data') %in% names(x))) {

      stop(paste('The input list needs to provide red_obj,',
                 'red_fun, component_tbl, loadings and data elements'),
           call. = FALSE)

    }

    stopifnot(rlang::is_quosure(x$data))

    ## output

    structure(x, class = 'red_analysis')

  }

# Simple clustering class -----

#' Construct a clust_analysis class object.
#'
#' @description Constructs a `clust_analysis` class object given a list with
#' results of a clustering analysis.
#' @details A named list with the following elements is required:
#' data (a quosure calling the original data set),
#' dist_mtx (a matrix with the distances between the observations),
#' dist_method (name of the distance statistic),
#' clust_obj (the output object of the clustering analysis),
#' clust_fun (the name of the clustering function),
#' clust_assignment (a tibble with the cluster assignment of the
#' observations with the observation and clust_id variables),
#' dots (additional arguments passed to the clustering function).
#' @param x a named list, see details.
#' @return a `clust_analysis` object.
#' @export

  clust_analysis <- function(x) {

    ## entry control

    if(!is.list(x)) {

      stop('Only named lists can be converted to a clust_analysis objects',
           call. = FALSE)

    }

    if(is.null(names(x))) {

      stop('Only named lists can be converted to a clust_analysis objects',
           call. = FALSE)

    }

    if(!all(c('data',
              'dist_mtx',
              'dist_method',
              'clust_obj',
              'clust_fun',
              'clust_assignment',
              'dots') %in% names(x))) {

      stop(paste('The input list needs to provide data, dist_mtx,',
                 'dist_method, clust_fun, clust_obj, clust_assignment',
                 'and dots elements'),
           call. = FALSE)

    }

    stopifnot(rlang::is_quosure(x$data))

    stopifnot(is.matrix(x$dist_mtx))

    stopifnot(x$dist.method %in% get_kernel_info())

    stopifnot(x$clust_fun %in% c('hclust',
                                 'kmeans',
                                 'dbscan',
                                 'som',
                                 'pam',
                                 'prediction'))

    stopifnot(is_tibble(x$clust_assignment))

    stopifnot(all(c('observation',
                    'clust_id') %in% names(x$clust_assignment)))

    ## output

    structure(x, class = 'clust_analysis')

  }

# Combined reduction analysis clustering class -----

#' Construct a combi_analysis class object.
#'
#' @description Constructs a `combi_anlysis()` class object given a list with
#' results of reduction analysis or self-organizing map and clustering analysis.
#' @details A named list with the following elements is required:
#' clust_analyses (`red_analysis` or `clust_analysis` objects),
#' clust_assignment (a tibble with the cluster assignment of the
  #' observations with the observation and clust_id variables).
#' @param x a named list, see details.
#' @return a `combi_analysis` object.
#' @export

  combi_analysis <- function(x) {

    ## entry control

    if(!is.list(x)) stop('Only named lists can be converted to a clust_analysis objects',
                         call. = FALSE)

    if(is.null(names(x))) stop('Only named lists can be converted to a clust_analysis objects',
                               call. = FALSE)

    if(!all(c('clust_analyses', 'clust_assignment') %in% names(x))) {

      stop('The input list needs to provide clust_analyses and clust_assignment elements',
           call. = FALSE)

    }

    stopifnot(is_tibble(x$clust_assignment))

    stopifnot(all(c('observation',
                    'clust_id') %in% names(x$clust_assignment)))

    ## output

    structure(x, class = 'combi_analysis')

  }

# Importance table class -----

#' Construct an importance object.
#'
#' @description Constructs an object of class `importance` on the top of
#' a tibble with the importance testing results.
#' @param x a data frame with the columns 'total_wss', 'total_ss', 'between_ss',
#' 'frac_var', 'variable' and 'frac_diff'.
#' @return a tibble of the `importance` class.
#' @export

  importance <- function(x) {

    if(!is.data.frame(x)) {

      stop('Please provide a data frame as an input.',
           call. = FALSE)

    }

    if(!all(c('total_wss',
              'total_ss',
              'between_ss',
              'frac_var',
              'variable',
              'frac_diff') %in% names(x))) {

      stop(paste('The input data frame needs to have total_wss, total_ss,',
                 'between_ss, frac_var, variable and frac_diff columns.'),
           call. = FALSE)

    }

    ## output

    x <- as_tibble(x)

    structure(x, class = c('importance', class(x)))

  }

# Cross-distance class -------

#' Construct a cross_dist object.
#'
#' @description
#' The `cross_dist` class objects a created on the top of a named list
#' of matrices of cross-distances between clusters.
#' The list elements have to be named after the compared clusters following
#' the 'clust1 vs clust2' scheme.
#'
#' @param x a named list of cross-distance matrices.
#' @param type type of the cross-distances:
#' `homologous` for comparison of within the same clustering structure or
#' `heterologous` for comparison of two clustering structures.
#' @param method name of the distance metric.
#' @param x_levels order of the clusters of the
#' first cluster/combi analysis object.
#' @param y_levels order of the clusters of the
#' second cluster/combi analysis object.
#' @param ... extra arguments, currently none defined.
#'
#' @return an object of class `cross-distance` being a list of
#' cross-distance matrices.
#' Information on the comparison type and distance
#' metric are stored as the `type` and `dist_method` attributes.
#' Information on the cluster order is stored as the
#' `x_levels` and `y_levels` attributes.
#'
#' @export

  cross_dist <- function(x,
                         type = c('homologous', 'heterologous'),
                         method = 'euclidean',
                         x_levels = NULL,
                         y_levels = NULL, ...) {

    ## entry control -------

    err_txt <- "The 'x' argument has to be a named list of numeric matrices."

    if(!is.list(x)) stop(err_txt, call. = FALSE)

    mtx_check <- map_lgl(x, is.matrix)

    if(!all(mtx_check)) stop(err_txt, call. = FALSE)

    num_check <- map_lgl(x, is.numeric)

    if(!all(num_check)) stop(err_txt, call. = FALSE)

    type <- match.arg(type[1], c('homologous', 'heterologous'))

    if(!method %in% get_kernel_info()) {

      stop(paste("Unsupported distance measure.",
                 "Please consult 'get_kernel_info()."),
           call. = FALSE)

    }

    ## object -------

    x <- structure(x, class = c('cross_dist', class(x)))

    attr(x, 'type') <- type
    attr(x, 'dist_method') <- method
    attr(x, 'x_levels') <- x_levels
    attr(x, 'y_levels') <- y_levels

    x

  }

# Silhouette ------

#' Generate a sil_extra object.
#'
#' @description
#' Extends the \code{\link[cluster]{silhouette}} object by
#' cluster order and names.
#'
#' @details
#' The `sil_extra` class has \code{\link{summary.sil_extra}}
#' and \code{\link{plot.sil_extra}}
#' methods compatible with the tidyverse environment.
#'
#' @param x an object of the \code{\link[cluster]{silhouette}} class.
#' @param assignment an data frame with the `clust_id` and `observation`
#' columns defining the cluster assignment, e.g. obtained by the
#' \code{\link{extract}} function applied to a `clust_analysis` or
#' `combi_analysis` object.
#'
#' @return an object of the `sil_extra` class. Technically, a tibble with
#' the observation ID, cluster name, neighbor cluster name and silhouette width.
#'
#' @export

  sil_extra <- function(x, assignment) {

    ## entry control -----

    if(!inherits(x, 'silhouette')) {

      stop("'x' must be a 'silhouette' class object.",
           call. = FALSE)

    }

    err_txt <- paste("'assignment' needs to be a data frame with the
                     columns 'observation' and 'clust_id'.")

    if(!is.data.frame(assignment)) stop(err_txt, call. = FALSE)

    if(!all(c('observation', 'clust_id') %in% names(assignment))) {

      stop(err_txt, call. = FALSE)

    }

    if(nrow(x) != nrow(assignment)) {

      stop("The number of rows in 'x' and 'assignment' must be equal.",
           call. = FALSE)

    }

    ## the output tibble ------

    clust_levs <- levels(assignment$clust_id)

    observation <- NULL
    clust_id <- NULL
    neighbor_id <- NULL
    sil_width <- NULL

    sil_extra <-
      tibble(observation = assignment$observation,
             clust_id = clust_levs[x[, 'cluster']],
             neighbor_id = clust_levs[x[, 'neighbor']],
             sil_width = x[, 'sil_width'])

    sil_extra <- mutate(sil_extra,
                        clust_id = factor(clust_id, clust_levs),
                        neighbor_id = factor(neighbor_id, clust_levs))

    structure(sil_extra,
              class = c('sil_extra', class(sil_extra)))

  }

# END ------
