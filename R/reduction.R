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
#' clust_object. Hence, some distances may crash the analysis with UMAP, see:
#' \code{\link[umap]{umap.defaults}} for the compatible distances.
#' For `combi_analysis` objects, the analysis is done for the global clustering,
#' i.e. assignment of observations to the clusters and not to the SOM nodes.
#'
#' @param object an object.
#' @param kdim number of dimensions. If NULL, kdim is set to the number of
#' clusters.
#' @param red_fun reduction analysis function: 'pca' (PCA), 'mds' (MDS) or
#' 'umap' (UMAP).
#' @param with type of the input data for the reduction analysis:
#' the clustering data ('data') or the matrix of distances ('distance').
#' @param ... extra arguments passed to \code{\link{reduce_data}}.
#'
#' @return a `red_analysis` object with the component/score table containing
#' the cluster assignment information ('clust_id' variable).
#'
#' @export components.clust_analysis
#' @export

  components.clust_analysis <- function(object,
                                        kdim = NULL,
                                        red_fun = c('pca', 'mds', 'umap'),
                                        with = c('distance', 'data'), ...) {

    ## entry control -------

    stopifnot(is_clust_analysis(object))

    with <- match.arg(with[1], c('distance', 'data', 'umap'))

    red_fun <- match.arg(red_fun[1], c('pca', 'mds', 'umap'))

    if(is.null(kdim)) {

      kdim <- length(unique(object$clust_assignment$clust_id))

    }

    if(!is.numeric(kdim)) stop("'kdim' must be numeric.", call. = FALSE)

    ## reduction --------

    red_obj <- reduce_data(extract.clust_analysis(object, type = with),
                           distance_method = object$dist_method,
                           kdim = kdim,
                           red_fun = red_fun, ...)

    red_obj$component_tbl <- left_join(red_obj$component_tbl,
                                       object$clust_assignment,
                                       by = 'observation')

    red_obj

  }

#' @rdname components.clust_analysis
#' @export components.combi_analysis
#' @export

  components.combi_analysis <- function(object,
                                        kdim = NULL,
                                        red_fun = c('pca', 'mds', 'umap'),
                                        with = c('distance', 'data'), ...) {

    ## entry control -------

    stopifnot(is_combi_analysis(object))

    with <- match.arg(with[1], c('distance', 'data'))

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


    red_analysis$component_tbl <- left_join(red_analysis$component_tbl,
                                            object$clust_assignment,
                                            by = 'observation')

    red_analysis

  }

# END ------
