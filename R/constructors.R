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

      stop('The input list needs to provide red_obj, red_fun, component_tbl, loadings and data elements',
           call. = FALSE)

    }

    stopifnot(rlang::is_quosure(x$data))

    ## output

    structure(x, class = 'red_analysis')

  }

# Simple clustering class -----

#' Construct a clust_analysis class object.
#'
#' @description Constructs a clust_analysis class object given a list storing
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

      stop('The input list needs to provide data, dist_mtx, dist_method, clust_fun, clust_obj, clust_assignment and dots elements', call. = FALSE)

    }

    stopifnot(rlang::is_quosure(x$data))

    stopifnot(is.matrix(x$dist_mtx))

    stopifnot(x$dist.method %in% clustTools::get_kernel_info())

    stopifnot(x$clust_fun %in% c('hclust',
                                 'kmeans',
                                 'dbscan',
                                 'som',
                                 'pam',
                                 'prediction'))

    stopifnot(tibble::is_tibble(x$clust_assignment))

    stopifnot(all(c('observation',
                    'clust_id') %in% names(x$clust_assignment)))

    ## output

    structure(x, class = 'clust_analysis')

  }

# Combined reduction analysis clustering class -----

#' Construct a combi_analysis class object.
#'
#' @description Constructs a clust_analysis class object given a list storing
#' results of reduction analysis or self-organizing map and clustering analysis.
#' @details A named list with the following elements is required:
#' clust_analyses (sores red_analysis or clust_analysis objects),
#' clust_assignment (a tibble with the cluster assignment of the
  #' observations with the observation and clust_id variables).
#' @param x a named list, see details.
#' @export

  combi_analysis <- function(x) {

    ## entry control

    if(!is.list(x)) stop('Only named lists can be converted to a clust_analysis objects', call. = FALSE)

    if(is.null(names(x))) stop('Only named lists can be converted to a clust_analysis objects', call. = FALSE)

    if(!all(c('clust_analyses', 'clust_assignment') %in% names(x))) {

      stop('The input list needs to provide clust_analyses and clust_assignment elements')

    }

    stopifnot(tibble::is_tibble(x$clust_assignment))

    stopifnot(all(c('observation',
                    'clust_id') %in% names(x$clust_assignment)))

    ## output

    structure(x, class = 'combi_analysis')

  }

# Importance table class -----

#' Construct an importance object.
#'
#' @description Constructs an object of class 'importance' on the top of
#' a tibble with the importance testing results.
#' @param x a data frame with the columns 'total_wss', 'total_ss', 'between_ss',
#' 'frac_var', 'variable' and 'frac_diff'.
#' @export

  importance <- function(x) {

    if(!'data.frame' %in% class(x)) {

      stop('Please provide a data frame as an input.', call. = FALSE)

    }

    if(!all(c('total_wss',
              'total_ss',
              'between_ss',
              'frac_var',
              'variable',
              'frac_diff') %in% names(x))) {

      stop('The input data frame needs to have total_wss, total_ss, between_ss, frac_var, variable and frac_diff columns.',
           call. = FALSE)

    }

    ## output

    x <- tibble::as_tibble(x)

    structure(x, class = c('importance', class(x)))

  }


# END ------
