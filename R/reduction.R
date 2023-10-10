# Dimensionality reduction for results of clustering analyses

#' Dimensionality reduction analysis of the analysis data or distance matrix.
#'
#' @description
#' Performs principal component analysis (PCA), multi-dimensional
#' scaling (MDS) or uniform manifold approximation and projection (UMAP) of the
#' analysis data set used for clustering or distance matrix.
#'
#' @details
#' See \code{\link{reduce_data}} for the implementation details.
#' The distance method, relevant for MDS and UMAP. is taken over from the
#' `clust_object`. Hence, some distances may crash the analysis with UMAP, see:
#' \code{\link[umap]{umap.defaults}} for the compatible distances.
#' For `combi_analysis` objects, the analysis is done for the global clustering,
#' i.e. assignment of observations to the clusters and not to the SOM nodes.
#' In cases, when the clustering analysis was done with an user-provided
#' dissimilarity object (subclass `min_analysis` of `clust_analysis`) it is not
#' possible to perform a dimensionality reduction analysis with the genuine data
#' set - the analysis can be performed only for the distance matrix.
#' In case of multi-layer SOM analyses or combined multi-layer SOM - clustering
#' analyses and `with` set to 'data', the dimensionality reduction analysis
#' is done separately for each data layer.
#'
#' @param object an object.
#' @param kdim number of dimensions. If NULL, kdim is set to the number of
#' clusters.
#' @param red_fun reduction analysis function: 'pca' (PCA), 'mds' (MDS) or
#' 'umap' (UMAP).
#' @param with type of the input data for the reduction analysis:
#' the clustering data ('data'), the matrix of distances between observations
#' ('distance') or U matrix between SOM nodes ('umatrix').
#' @param ... extra arguments passed to \code{\link{reduce_data}}.
#'
#' @return a `red_analysis` object with the component/score table containing
#' the cluster assignment information (`clust_id` variable). In case of
#' multi-layer SOM analyses with the clustering data set, a list of
#' `red_analysis` objects is returned.
#'
#' @export components.clust_analysis
#' @export

  components.clust_analysis <- function(object,
                                        kdim = NULL,
                                        red_fun = c('pca', 'mds', 'umap'),
                                        with = c('distance', 'data', 'umatrix'), ...) {

    ## entry control -------

    stopifnot(is_clust_analysis(object))

    with <- match.arg(with[1], c('distance', 'data', 'umatrix'))

    red_fun <- match.arg(red_fun[1], c('pca', 'mds', 'umap'))

    if(is.null(kdim)) {

      kdim <- length(unique(object$clust_assignment$clust_id))

    }

    if(!is.numeric(kdim)) stop("'kdim' must be numeric.", call. = FALSE)

    check_supersom <-
      object$clust_fun %in% c('supersom', 'supersom_prediction')

    if(with == 'umatrix' & !object$clust_fun %in% c('som', 'supersom')) {

      warning("The 'umatrix' option is only available for SOM analyses.",
              call. = FALSE)

      return(NULL)

    }

    observation <- NULL
    node <- NULL
    clust_id <- NULL

    ## reduction: most cases --------
    ## concerns clustering without multi-layer SOM and distances

    if(!check_supersom | with %in% c('distance', 'umatrix')) {

      red_obj <- reduce_data(extract.clust_analysis(object, type = with),
                             distance_method = object$dist_method,
                             kdim = kdim,
                             red_fun = red_fun, ...)

      if(with == 'umatrix') {

        red_obj$component_tbl <-
          mutate(red_obj$component_tbl,
                 node = factor(observation),
                 clust_id = node)

      } else {

        red_obj$component_tbl <- left_join(red_obj$component_tbl,
                                           object$clust_assignment,
                                           by = 'observation')

      }

      return(red_obj)

    }

    ## reduction: multi-layer SOM and its predictions ---------

    if(object$clust_fun %in% c('supersom')) {

      layer_dist <- object$clust_obj$dist.fcts

    } else {

      layer_dist <- object$dots$dist.fcts

    }

    red_objects <-
      pmap(list(data = model.frame(object),
                distance_method = layer_dist),
           reduce_data,
           kdim = kdim,
           red_fun = red_fun, ...)

    for(i in seq_along(red_objects)) {

      red_objects[[i]]$component_tbl <-
        left_join(red_objects[[i]]$component_tbl,
                  object$clust_assignment,
                  by = 'observation')

    }

    red_objects

  }

#' @rdname components.clust_analysis
#' @export

  components.min_analysis <- function(object,
                                      kdim = NULL,
                                      red_fun = c('pca', 'mds', 'umap'), ...) {
    NextMethod()

  }
#' @rdname components.clust_analysis
#' @export components.combi_analysis
#' @export

  components.combi_analysis <- function(object,
                                        kdim = NULL,
                                        red_fun = c('pca', 'mds', 'umap'),
                                        with = c('distance', 'data', 'umatrix'), ...) {

    ## entry control -------

    stopifnot(is_combi_analysis(object))

    with <- match.arg(with[1], c('distance', 'data', 'umatrix'))

    red_fun <- match.arg(red_fun[1], c('pca', 'mds', 'umap'))

    ## reduction analysis -------

    red_analysis <- components(object$clust_analyses$observation,
                               kdim = kdim,
                               red_fun = red_fun,
                               with = with, ...)

    observation <- NULL

    red_analysis$component_tbl <- select(red_analysis$component_tbl,
                                         observation,
                                         dplyr::starts_with('comp'))

    if(with == 'umatrix') {

      red_analysis$component_tbl <-
        mutate(red_analysis$component_tbl,
               node = observation,
               clust_id = factor(observation))

    } else {

      red_analysis$component_tbl <- left_join(red_analysis$component_tbl,
                                              object$clust_assignment,
                                              by = 'observation')

    }

    red_analysis

  }

#' @rdname components.clust_analysis
#' @export

  components.umatrix_analysis <- function(object,
                                          kdim = NULL,
                                          red_fun = c('pca', 'mds', 'umap'),
                                          with = c('distance', 'data', 'umatrix'), ...) {

    ## entry control ------

    stopifnot(is_umatrix_analysis(object))

    red_fun <- match.arg(red_fun[1],
                         c('pca', 'mds', 'umap'))

    with <- match.arg(with[1],
                      c('distance', 'data', 'umatrix'))

    if(with != 'data') return(NextMethod())

    data_comps <- components(object$clust_analyses$observation,
                             kdim = kdim,
                             red_fun = red_fun,
                             with = 'data')

    observation <- NULL

    for(i in seq_along(data_comps)) {

      data_comps[[i]]$component_tbl <-
        select(data_comps[[i]]$component_tbl,
               observation,
               dplyr::starts_with('comp'))

      data_comps[[i]]$component_tbl <-
        left_join(data_comps[[i]]$component_tbl,
                  object$clust_assignment,
                  by = 'observation')

    }

    data_comps

  }

# END ------
