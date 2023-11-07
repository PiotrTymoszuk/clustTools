# Tools for assessment of the projection of neighborhood into the clusters

# Neighborhood preservation -------

#' Neighborhood preservation.
#'
#' @description
#'
#' Checks the fraction of the nearest neighbors of data points, that are
#' located in the same cluster or in the node of a self-organizing map (SOM).
#'
#' @details
#' The function computes fractions of nearest neighbors of each data point that
#' are preserved after dimensionality reduction or clustering analysis.
#' In case of reduction analyses, the neighborhoods in the genuine data set and
#' the data set following dimensionality reduction ('layout') are compared.
#' In case of simple clustering analyses, the function calculates fractions of
#' the k-nearest neighbors located in k-nearest clusters of the given data
#' point. By default, the number of nearest clusters is set to one
#' (`kNN_cluster`), which means that `np()` simply checks which fraction of the
#' nearest data point neighbors is played in the same cluster. In such case,
#' the neighborhood preservation fraction averaged for the entire clustering
#' object and particular clusters gives a measure of cluster separation, similar
#' to \code{\link[cluster]{silhouette}}.
#'
#' @references
#' Venna J, Kaski S. Neighborhood preservation in nonlinear projection methods:
#' An experimental study. Lect Notes Comput Sci (including Subser Lect Notes
#' Artif Intell Lect Notes Bioinformatics) (2001) 2130:485–491.
#' doi:10.1007/3-540-44668-0_68#'
#' @references
#' Breard GT. Evaluating Self-Organizing Map Quality Measures as Convergence
#' Criteria Criteria. Open Access Master’s Theses. Paper 1033.
#' Available at: https://digitalcommons.uri.edu/theses/1033
#'
#' @param x an object.
#' @param kNN_data number of k-nearest neighbors for data points.
#' @param kNN_cluster number of k-nearest neighbors for clusters. If `NULL`,
#' the values will be determined automatically. In this case `kNN_cluster = 1`
#' for non-SOM cluster analyses and `kNN_cluster = kNN_data` for SOM analyses.
#' @param type type of data used for calculation of the neighborhood
#' preservation. For `type = 'data'`, the comparison of neighborhoods is done
#' between the data points and their SOM node assignment. For `type = 'node'`,
#' neighborhoods of the nodes and final clusters are compared.
#' For `type = 'final'`, the analysis is done at the top level, i.e.
#' neighborhoods of the data points and the final cluster assignment
#' are evaluated.
#' @param ... extra arguments passed to methods.
#'
#' @return An object of the \code{\link{knb}} class with
#' \code{\link{summary.knb}} and
#' \code{\link{plot.knb}} methods.
#'
#' @export

  np <- function(x, ...) UseMethod('np')

#' @rdname np
#' @export

  np.clust_analysis <- function(x,
                                kNN_data = 5,
                                kNN_cluster = NULL, ...) {

    ## input control -------

    stopifnot(is_clust_analysis(x))
    stopifnot(is.numeric(kNN_data))

    kNN_data <- as.integer(kNN_data)

    if(is.null(kNN_cluster)) {

      if(x$clust_fun %in% c('som', 'supersom')) {

        kNN_cluster <- kNN_data

      } else {

        kNN_cluster <- 1

      }

    }

    stopifnot(is.numeric(kNN_cluster))

    kNN_cluster <- as.integer(kNN_cluster)

    ## cluster assignment and distances -------

    clust_assignment <- extract(x, 'assignment')

    data_dist <- dist(x, type = 'distance')

    if(x$clust_fun %in% c('som', 'supersom')) {

      clust_dist <- dist(x, type = 'umatrix')

    } else if(kNN_cluster > 1) {

      clust_dist <- cross2dist(cross_distance(x), zero_diag = TRUE)

    } else {

      clust_dist <- NULL

    }

    ## neighborhood comparison ------

    n_preservation(data_dist = data_dist,
                   clust_assignment = clust_assignment,
                   clust_dist = clust_dist,
                   kNN_data = kNN_data,
                   kNN_cluster = kNN_cluster)

  }

#' @rdname np
#' @export

  np.combi_analysis <- function(x,
                                kNN_data = 5,
                                kNN_cluster = NULL,
                                type = c('data', 'node', 'final'), ...) {

    ## entry check -------

    stopifnot(is_combi_analysis(x))
    stopifnot(is.numeric(kNN_data))

    type <- match.arg(type[1], c('data', 'node', 'final'))

    ## observation and node neighborhood ------

    if(type != 'final') {

      return(switch(type,
                    data = np(x$clust_analyses$observation,
                              kNN_data = kNN_data,
                              kNN_cluster = kNN_cluster),
                    node = np(x$clust_analyses$node,
                              kNN_data = kNN_data,
                              kNN_cluster = kNN_cluster)))

    }

    ## final neighborhood -------

    if(is.null(kNN_cluster)) kNN_cluster <- 1

    kNN_cluster <- as.integer(kNN_cluster)

    stopifnot(is.numeric(kNN_cluster))

    data_dist <- dist(x, type = 'distance')[['observation']]

    clust_dist <- cross2dist(cross_distance(x), zero_diag = TRUE)

    n_preservation(data_dist = data_dist,
                   clust_assignment = extract(x, 'assignment'),
                   clust_dist = clust_dist,
                   kNN_data = kNN_data,
                   kNN_cluster = kNN_cluster)

  }

#' @rdname np
#' @export

  np.red_analysis <- function(x,
                              kNN_data = 5, ...) {

    stopifnot(is_red_analysis(x))

    np_reduction(data_dist = dist(x, type = 'distance'),
                 layout_dist = dist(x, type = 'layout'),
                 kNN_data = kNN_data)

  }

# Topology error --------

#' Topology error of self-organizing maps.
#'
#' @description
#' Calculates the topology error for data points in a self-organizing map (SOM).
#'
#' @details
#' The procedure of topology error computation is as follows: for each
#' observation, two nearest self-organizing map (SOM) nodes are identified.
#' If such nodes are neighbors in the initial layout of SOM prior to data
#' fitting, correct topology (coded with 0) is returned and an error otherwise
#' (coded as 1).
#' The function returns NULL with a warning when called for a non-SOM analysis
#' object.
#' `te()` is a S3 generic function.
#'
#' @return An object of the \code{\link{knb}} class with
#' \code{\link{summary.knb}} and \code{\link{plot.knb}} methods.
#'
#' @param x a `clust_analysis` or `combi_analysis` object.
#' @param type type reference clusters. For `type = 'node'`, topology error
#' within SOM nodes is computed, for `type = 'final'`, topology error in the
#' final clusters (clusters of SOM nodes) is calculated.
#' @param ... extra arguments passed to methods.
#'
#' @return An object of the \code{\link{knb}} class with
#' \code{\link{summary.knb}} and \code{\link{plot.knb}} methods.
#'
#' @export

  te <- function(x, ...) UseMethod('te')

#' @rdname te
#' @export

  te.clust_analysis <- function(x, ...) {

    ## input control ------

    stopifnot(is_clust_analysis(x))

    if(!x$clust_fun %in% c('som', 'supersom')) {

      warning('Topology error are returned only for SOM analysis objects.',
              call. = FALSE)

      return(NULL)

    }

    return(check_topo(x))

  }

#' @rdname te
#' @export

  te.combi_analysis <- function(x, type = c('node', 'final'), ...) {

    stopifnot(is_combi_analysis(x))

    type <- match.arg(type[1], c('node', 'final'))

    topo_res <- check_topo(x$clust_analyses$observation)

    if(type == 'node') return(topo_res)

    clust_assignment <- extract(x, 'assignment')[c('node', 'clust_id')]

    node <- NULL
    clust_id <- NULL

    clust_assignment <- filter(clust_assignment, !duplicated(node))

    topo_res <- mutate(topo_res,
                       node = as.character(clust_id))

    topo_res <- left_join(select(topo_res, -clust_id),
                          clust_assignment,
                          by = 'node')

    select(topo_res, -node)

  }

# END -----
