# Provides functions for dimensionality reduction analysis.

# Utils ------

#' Perform dimensionality reduction.
#'
#' @description Performs
#' multi-dimensional scaling (`mds()` via ' \code{\link[stats]{cmdscale}}),
#' principal component analysis (`pca()` via \code{\link[pcaPP]{PCAproj}}),
#' UMAP (`umap()` via \code{\link[umap]{umap}}),
#' or factor analysis (`fa()` via \code{\link[stats]{factanal}}).
#'
#' @references
#' McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and
#' Projection for Dimension Reduction. (2018)
#' Available at: https://arxiv.org/abs/1802.03426v3
#' @references
#' Croux C, Filzmoser P, Oliveira MR. Algorithms for Projection-Pursuit robust
#' principal component analysis. Chemom Intell Lab Syst (2007) 87:218–225.
#' doi:10.1016/j.chemolab.2007.01.004
#' @references
#' BARTLETT MS. THE STATISTICAL CONCEPTION OF MENTAL FACTORS. Br J Psychol
#' Gen Sect (1937) 28:97–104. doi:10.1111/j.2044-8295.1937.tb00863.x
#'
#' @details
#' UMAP parameters such as dimension number or distance are provided
#' as a \code{\link[umap]{umap.defaults}} object.
#'
#' @param data a data frame or a distance object(class `dist`).
#' @param distance_method name of the distance metric, see:
#' \code{\link{get_kernel_info}}. Ignored if `data` is a distance object.
#' @param kdim dimension number.
#' @param ... extra arguments passed to
#' \code{\link[stats]{cmdscale}}, \code{\link[pcaPP]{PCAproj}},
#' \code{\link[umap]{umap}} or \code{\link[stats]{factanal}}.
#'
#' @return an object of the class \code{\link{red_analysis}}.

  mds <- function(data,
                  distance_method = 'euclidean',
                  kdim = 2, ...) {

    ## the entry control is managed by the upstream analysis function.

    ## analysis

    model_frame <- enexpr(data)

    if(!inherits(data, 'dist')) {

      dist_mtx <- calculate_dist(data = data,
                                 method = distance_method)

      dist_mtx <- as.dist(dist_mtx)

    } else {

      dist_mtx <- data

      distance_method <- 'custom'

    }

    component_tbl <- stats::cmdscale(d = dist_mtx,
                                     k = kdim, ...)

    component_tbl <- as.data.frame(component_tbl)

    component_tbl <- mutate(component_tbl,
                            observation = rownames(as.data.frame(as.matrix(data))))

    component_tbl <- as_tibble(component_tbl)

    component_tbl <- set_names(component_tbl,
                               c(paste0('comp_', 1:kdim)),
                               'observation')

    ## output

    red_analysis(list(data = quo(!!model_frame),
                      red_obj = NULL,
                      red_fun = 'mds',
                      dist_method = distance_method,
                      component_tbl = component_tbl,
                      loadings = NULL))

  }

#' @rdname mds

  pca <- function(data, kdim = 2, ...) {

    ## the entry control is accomplished by the upstream function

    ## analysis

    model_frame <- enexpr(data)

    if(inherits(data, 'dist')) {

      model_data <- as.matrix(data)

    } else {

      model_data <- data

    }

    red_obj <- pcaPP::PCAproj(x = model_data,
                              k = kdim, ...)

    component_tbl <- red_obj$scores

    rownames(component_tbl) <- rownames(model_data)

    component_tbl <- as.data.frame(component_tbl)

    component_tbl <- mutate(component_tbl,
                            observation = rownames(model_data))

    component_tbl <- as_tibble(component_tbl)

    component_tbl <- set_names(component_tbl,
                               c(paste0('comp_', 1:kdim)),
                               'observation')

    loadings <- as.data.frame(unclass(red_obj$loadings))

    loadings <- rownames_to_column(loadings, 'variable')

    loadings <- as_tibble(loadings)

    loadings <- set_names(loadings,
                          c('variable',
                            paste0('comp_', 1:kdim)))

    ## output

    red_analysis(list(data = quo(!!model_frame),
                      red_obj = red_obj,
                      red_fun = 'pca',
                      dist_method = 'custom',
                      component_tbl = component_tbl,
                      loadings = loadings))

  }

#' @rdname mds

  umap <- function(data, distance_method, kdim, ...) {

    ## the entry control is accomplished by the upstream function

    ## analysis

    model_frame <- enexpr(data)

    if(inherits(data, 'dist')) {

      distance_method <- 'custom'

      ## metric name provided only for the umap()
      ## function to work - it won't be used anyway, because
      ## the data argument is a distance object

      metric <- 'euclidean'

      input <- 'dist'

    } else if(!distance_method %in% c('euclidean', 'manhattan', 'cosine')){

      data <- stats::as.dist(calculate_dist(data, distance_method))

      ## as above

      metric <- 'euclidean'

      input <- 'dist'

    } else {

      input <- 'data'

      metric <- distance_method

    }

    red_obj <- umap::umap(d = as.matrix(data),
                          n_components = kdim,
                          metric = metric,
                          input = input, ...)

    component_tbl <- as.data.frame(red_obj$layout)

    component_tbl <- mutate(component_tbl,
                            observation = rownames(as.data.frame(as.matrix(data))))

    component_tbl <- as_tibble(component_tbl)

    component_tbl <- set_names(component_tbl,
                               c(paste0('comp_', 1:(ncol(component_tbl) - 1)),
                                 'observation'))

    ## output

    red_analysis(list(data = quo(!!model_frame),
                      red_obj = red_obj,
                      red_fun = 'umap',
                      dist_method = distance_method,
                      component_tbl = component_tbl,
                      loadings = NULL))

  }

#' @rdname mds

  fa <- function(data, kdim = 2, ...) {

    ## the entry control is accomplished by the upstream function

    ## analysis

    model_frame <- enexpr(data)

    if(inherits(data, 'dist')) {

      warning('Factor analysis with distances in not recommended.',
              call. = FALSE)

      model_data <- as.matrix(data)

    } else {

      model_data <- data

    }

    red_obj <- stats::factanal(x = model_data,
                               factors = kdim, ...)

    if(!is.null(red_obj$scores)) {

      component_tbl <- red_obj$scores

      rownames(component_tbl) <- rownames(model_data)

      component_tbl <- as.data.frame(component_tbl)

      component_tbl <- mutate(component_tbl,
                              observation = rownames(model_data))

      component_tbl <- as_tibble(component_tbl)

      component_tbl <- set_names(component_tbl,
                                 c(paste0('comp_', 1:kdim)),
                                 'observation')

    } else {

      component_tbl <- NULL

    }

    loadings <- as.data.frame(unclass(red_obj$loadings))

    loadings <- rownames_to_column(loadings, 'variable')

    loadings <- as_tibble(loadings)

    loadings <- set_names(loadings,
                          c('variable',
                            paste0('comp_', 1:kdim)))

    ## output

    red_analysis(list(data = quo(!!model_frame),
                      red_obj = red_obj,
                      red_fun = 'fa',
                      dist_method = 'custom',
                      component_tbl = component_tbl,
                      loadings = loadings))

  }

# Upstream analysis function ------

#' Dimensionality reduction of a data set.
#'
#' @description
#' Performs dimensionality reduction of a data frame with principal
#' component analysis (PCA), multi-dimensional scaling (MDS), Uniform Manifold
#' Approximation and Projection (UMAP) or factor analysis (FA).
#'
#' @details
#' A wrapper around \code{\link[pcaPP]{PCAproj}} (PCA),
#' \code{\link[stats]{cmdscale}} (MDS), \code{\link[umap]{umap}} (UMAP)
#' and \code{\link[stats]{factanal}} (FA). Note:
#' the distances and other UMAP parameters are specified by a
#' \code{\link[umap]{umap.defaults}} object. Hence, not all distance measures
#' returned by \code{\link{get_kernel_info}} are available for UMAP computation.
#'
#' @references
#' McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and
#' Projection for Dimension Reduction. (2018)
#' Available at: https://arxiv.org/abs/1802.03426v3
#' @references
#' Croux C, Filzmoser P, Oliveira MR. Algorithms for Projection-Pursuit robust
#' principal component analysis. Chemom Intell Lab Syst (2007) 87:218–225.
#' doi:10.1016/j.chemolab.2007.01.004
#' @references
#' BARTLETT MS. THE STATISTICAL CONCEPTION OF MENTAL FACTORS. Br J Psychol
#' Gen Sect (1937) 28:97–104. doi:10.1111/j.2044-8295.1937.tb00863.x
#'
#' @return a \code{\link{red_analysis}} object with, among others,
#' \code{\link{plot.red_analysis}}, \code{\link{summary.red_analysis}} and
#' \code{\link{predict.red_analysis}} methods.
#'
#' @param data a numeric data frame, a matrix or a distance object
#' (class `dist`).
#' @param distance_method name of the distance metric, see:
#' \code{\link{get_kernel_info}}. Explicitly used only by MDS and UMAP.
#' @param kdim dimension number.
#' @param red_fun name of the dimensionality reduction function.
#' @param ... extra arguments passed to \code{\link[pcaPP]{PCAproj}} (PCA),
#' \code{\link[stats]{cmdscale}} (MDS) and \code{\link[umap]{umap}} (UMAP),
#' like the \code{\link[umap]{umap.defaults}} object for UMAP,
#' \code{\link[stats]{factanal}} (FA).
#'
#' @export

  reduce_data <- function(data,
                          distance_method = 'euclidean',
                          kdim = 2,
                          red_fun = c('pca', 'mds', 'umap', 'fa'), ...) {

    ## entry control --------

    if(!inherits(data, 'dist')) {

      check_numeric(data)

      if(!distance_method %in% get_kernel_info()) {

        stop('Invalid distance method.', call. = FALSE)

      }

    }

    red_fun <- match.arg(red_fun[1], c('pca', 'mds', 'umap', 'fa'))

    kdim <- as.integer(kdim)

    ## analysis -------

    switch(red_fun,
           mds = mds(data,
                     distance_method,
                     kdim, ...),
           pca = pca(data,
                     kdim, ...),
           umap = umap(data,
                       distance_method,
                       kdim, ...),
           fa = fa(data,
                   kdim, ...))

  }

# SOM as dimensionality reduction method ------

#' Dimensionality reduction with self-organizing maps.
#'
#' @description
#' The `som_reduce()` function applies self-organizing maps (SOM) as a classical
#' dimensionality reduction method, i.e. to reduce the number of variables.
#'
#' @details
#' In contrast to \code{\link{som_cluster}}, which reduces the observation
#' number (observation -> nodes), `som_reduce()` diminishes the variable number
#' by representing associated variables as SOM nodes or 'meta-variables'.
#' The data frame of such meta-variables (to be strict, codebook vectors of
#' node positions) is stored as the score/component table of the `red_analysis`
#' object returned by the function.
#'
#' @return a \code{\link{red_analysis}} object with, among others,
#' \code{\link{plot.red_analysis}}, \code{\link{summary.red_analysis}} and
#' \code{\link{predict.red_analysis}} methods.
#'
#' @param data a numeric data frame or matrix.
#' @param distance_method name of the distance metric, see:
#' \code{\link{get_kernel_info}}.
#' @inheritParams som_cluster
#'
#' @export

  som_reduce <- function(data,
                         distance_method = 'euclidean',
                         xdim = 5,
                         ydim = 4,
                         topo = c('hexagonal', 'rectangular'),
                         neighbourhood.fct = c('gaussian', 'bubble'),
                         toroidal = FALSE,
                         seed = 1234, ...) {

    ## input control for the data and distance -------

    check_numeric(data)

    distance_method <- distance_method[1]

    if(!distance_method %in% get_kernel_info()) {

      stop('Invalid distance method.', call. = FALSE)

    }

    model_frame <- enexpr(data)

    ## entry control for the grid parameters -----

    stopifnot(is.numeric(xdim))
    stopifnot(is.numeric(ydim))

    xdim <- as.integer(xdim)
    ydim <- as.integer(ydim)

    topo <- match.arg(topo[1],
                      c('hexagonal', 'rectangular'))

    neighbourhood.fct <- match.arg(neighbourhood.fct[1],
                                   c('gaussian', 'bubble'))

    stopifnot(is.logical(toroidal))

    ## grid and RNG ----------

    som_grid <- kohonen::somgrid(xdim = xdim,
                                 ydim = ydim,
                                 topo = topo,
                                 neighbourhood.fct = neighbourhood.fct,
                                 toroidal = toroidal)

    set.seed(seed = seed)

    ## fitting the SOM ----------

    kohonen_obj <- kohonen::som(t(as.matrix(data)),
                                grid = som_grid,
                                dist.fcts = distance_method, ...)

    score_tbl <- as.data.frame(t(kohonen_obj$codes[[1]]))

    score_tbl <- set_names(score_tbl, paste0('comp_', 1:ncol(score_tbl)))

    score_tbl <- rownames_to_column(score_tbl, 'observation')

    ## the output red_analysis class object ------

    red_analysis(list(data = quo(!!model_frame),
                      red_obj = kohonen_obj,
                      red_fun = 'som',
                      dist_method = distance_method,
                      component_tbl = as_tibble(score_tbl),
                      loadings = NULL))

  }

# END -----
