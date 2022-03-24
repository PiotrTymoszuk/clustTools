# Provides functions for dimensionality reduction anlysis.

# Utils ------

#' Perform multi-dimensional scaling.
#'
#' @description Performs multi-dimensional scaling of the data. A wrapper around
#' \code{\link[stats]{cmdscale}}.
#' @param data a data frame.
#' @param distance_method name of the distance metric, see:
#' \code{\link{get_kernel_info}}
#' @param kdim dimension number.
#' @param ... extra arguments passed to \code{\link[stats]{cmdscale}}.

  mds <- function(data,
                  distance_method = 'euclidean',
                  kdim = 2, ...) {

    ## the entry control is managed by the upstream analysis function.

    ## analysis

    model_frame <- rlang::enexpr(data)

    dist_mtx <- clustTools::calculate_dist(data = data,
                                           method = distance_method)

    dist_mtx <- as.dist(dist_mtx)

    component_tbl <- stats::cmdscale(d = dist_mtx,
                                     k = kdim, ...)

    component_tbl <- as.data.frame(component_tbl)

    component_tbl <- dplyr::mutate(component_tbl,
                                   observation = rownames(data))

    component_tbl <- tibble::as_tibble(component_tbl)

    component_tbl <- rlang::set_names(component_tbl,
                                      c(paste0('comp_', 1:kdim)),
                                      'observation')

    ## output

    clustTools::red_analysis(list(data = rlang::quo(!!model_frame),
                                  red_obj = NULL,
                                  red_fun = 'mds',
                                  component_tbl = component_tbl,
                                  loadings = NULL))


  }

#' Principal component analysis.
#'
#' @description Performs principal component analysis of the data. A wrapper
#' around \code{\link[pcaPP]{PCAproj}}.
#' @param data a data frame.
#' @param kdim dimension number.
#' @param ... extra arguments passed to \code{\link[pcaPP]{PCAproj}}.

  pca <- function(data, kdim = 2, ...) {

    ## the entry control is accomplished by the upstream function

    ## analysis

    model_frame <- rlang::enexpr(data)

    red_obj <- pcaPP::PCAproj(x = data,
                              k = kdim, ...)

    component_tbl <- red_obj$scores

    rownames(component_tbl) <- rownames(data)

    component_tbl <- as.data.frame(component_tbl)

    component_tbl <- dplyr::mutate(component_tbl,
                                   observation = rownames(data))

    component_tbl <- tibble::as_tibble(component_tbl)

    component_tbl <- rlang::set_names(component_tbl,
                                      c(paste0('comp_', 1:kdim)),
                                      'observation')

    loadings <- as.data.frame(unclass(red_obj$loadings))

    loadings <- tibble::rownames_to_column(loadings, 'variable')

    loadings <- tibble::as_tibble(loadings)

    loadings <- rlang::set_names(loadings,
                                 c('variable',
                                   paste0('comp_', 1:kdim)))

    ## output

    clustTools::red_analysis(list(data = rlang::quo(!!model_frame),
                                  red_obj = red_obj,
                                  red_fun = 'pca',
                                  component_tbl = component_tbl,
                                  loadings = loadings))


  }

#' Uniform Manifold Approximation and Projection (UMAP).
#'
#' @description Performs UMAP of the data. A wrapper
#' around \code{\link[umap]{umap}}.
#' @details UMAP parameters such as dimension number or distance are provided
#' as a \code{\link[umap]{umap.defaults}} object.
#' @param data a data frame.
#' @param distance_method name of the distance metric specified by a
#' \code{\link[umap]{umap.defaults}} object.
#' @param kdim dimension number.
#' @param ... extra arguments passed to \code{\link[umap]{umap}}.

  umap <- function(data, distance_method, kdim, ...) {

    ## the entry control is accomplished by the upstream function

    ## analysis

    model_frame <- rlang::enexpr(data)

    red_obj <- umap::umap(d = as.matrix(data),
                          n_components = kdim,
                          metric = distance_method, ...)

    component_tbl <- as.data.frame(red_obj$layout)

    component_tbl <- dplyr::mutate(component_tbl,
                                   observation = rownames(data))

    component_tbl <- tibble::as_tibble(component_tbl)

    component_tbl <- rlang::set_names(component_tbl,
                                      c(paste0('comp_', 1:(ncol(component_tbl) - 1)),
                                        'observation'))

    ## output

    clustTools::red_analysis(list(data = rlang::quo(!!model_frame),
                                  red_obj = red_obj,
                                  red_fun = 'umap',
                                  component_tbl = component_tbl,
                                  loadings = NULL))

  }

# Upstream analysis function ------

#' Dimensionality reduction of a data set.
#'
#' @description Performs dimensionality reduction of a data frame with principal
#' component analysis (PCA), multi-dimensional scaling (MDS) or Uniform Manifold
#' Approximation and Projection (UMAP).
#' @details A wrapper around \code{\link[pcaPP]{PCAproj}} (PCA),
#' \code{\link[stats]{cmdscale}} (MDS) and \code{\link[umap]{umap}}. Note:
#' the distances and other UMAP parameters are specified by a
#' \code{\link[umap]{umap.defaults}} object. Hence, not all distance measures
#' returned by \code{\link{get_kernel_info}} are available for UMAP computation.
#' @return a 'red_analysis' object.
#' @param data a numeric data frame or a matrix.
#' @param distance_method name of the distance metric, see:
#' \code{\link{get_kernel_info}}. Valid only for MDS and UMAP.
#' For UMAP, the distance is specified by a
#' \code{\link[umap]{umap.defaults}} object.
#' @param kdim dimension number.
#' @param red_fun name of the dimensionality reduction function.
#' @param ... extra arguments passed to \code{\link[pcaPP]{PCAproj}} (PCA),
#' \code{\link[stats]{cmdscale}} (MDS) and \code{\link[umap]{umap}} (UMAP),
#' like the \code{\link[umap]{umap.defaults}} object for UMAP.
#' @export

  reduce_data <- function(data,
                          distance_method = 'euclidean',
                          kdim = 2,
                          red_fun = c('pca', 'mds', 'umap'), ...) {

    ## entry control

    clustTools:::check_numeric(data)

    red_fun <- match.arg(red_fun[1], c('pca', 'mds', 'umap'))

    if(!distance_method %in% clustTools::get_kernel_info()) {

      stop('Invalid distance method.', call. = FALSE)

    }

    kdim <- as.integer(kdim)

    ## analysis

    if(red_fun == 'mds') {

      return(clustTools:::mds(data,
                              distance_method,
                              kdim, ...))

    } else if(red_fun == 'pca') {

      return(clustTools:::pca(data,
                              kdim, ...))

    } else {

      return(clustTools:::umap(data,
                               distance_method,
                               kdim, ...))

    }

  }

# END -----
