# Utilities for neighborhood analyses

# Cross distance to distance matrix ------

#' Convert cross-distance objects to a distance matrix.
#'
#' @description
#' Converts a `cross_dist` class object to a matrix of mean
#' cross-distances between the clusters.
#'
#' @param x a `cross_dist` object.
#' @param zero_diag logical, should diagonals be filled with zeros?
#'
#' @return a `dist` object.
#'
#' @export

  cross2dist <- function(x, zero_diag = TRUE) {

    ## input control ------

    stopifnot(is_cross_dist(x))
    stopifnot(is.logical(zero_diag))

    ## conversion ------

    cluster1 <- NULL
    cluster2 <- NULL

    dist_data <-
      mutate(summary(x),
             cluster1 = as.character(cluster1),
             cluster2 = as.character(cluster2))

    upper_data <- filter(dist_data,
                         cluster1 != cluster2)

    prov_cluster <- NULL

    upper_data <-
      mutate(upper_data,
             prov_cluster = cluster1,
             cluster1 = cluster2,
             cluster2 = prov_cluster)

    upper_data <- select(upper_data, - prov_cluster)

    dist_data <- rbind(dist_data, upper_data)[c('cluster1', 'cluster2', 'mean')]

    clust_dist <- tidyr::pivot_wider(dist_data,
                                     names_from = 'cluster2',
                                     values_from = 'mean')

    clust_labels <- clust_dist$cluster1

    clust_dist <- as.matrix(column_to_rownames(clust_dist, 'cluster1'))

    clust_dist <- clust_dist[clust_labels, clust_labels]

    if(zero_diag) diag(clust_dist) <- 0

    as.dist(clust_dist)

  }

# Neighborhood preservation --------

#' Neighborhood preservation.
#'
#' @description
#' Computes neighborhood preservation stats.
#'
#' @details
#' For internal use. For computation details, see: \code{\link{np}}.
#'
#' @param data_dist a matrix of distances between data points.
#' @param layout_dist a matrix of distances between data points following
#' dimensionality reduction.
#' @param clust_assignment a data frame with the cluster assignment scheme.
#' @param clust_dist an optional matrix of distances between the clusters.
#' Ignored if `kNN_cluster = 1`.
#' @param kNN_data number of k-nearest neighbors for data points.
#' @param kNN_cluster number of k-nearest neighbors for clusters
#'
#' @return An object of the \code{\link{knb}} class with
#' \code{\link{summary.knb}} and
#' \code{\link{plot.knb}} methods.

  n_preservation <- function(data_dist,
                             clust_assignment,
                             clust_dist = NULL,
                             kNN_data = 5,
                             kNN_cluster = 1) {

    ## input control ------

    ## done by upstream functions

    stopifnot(inherits(data_dist, 'dist'))

    if(!is.null(clust_dist)) {

      stopifnot(inherits(clust_dist, 'dist'))

    }

    ## nearest data points and their node assignment ------

    data_labels <- rownames(as.matrix(data_dist))

    clust_vct <-
      set_names(as.character(clust_assignment$clust_id),
                clust_assignment$observation)

    data_kNN <- dbscan::kNN(data_dist, k = kNN_data)$id

    data_neighbors <-
      map(rownames(data_kNN),
          ~data_labels[data_kNN[.x, ]])

    data_neighbors <- set_names(data_neighbors, rownames(data_kNN))

    data_neighbor_assignment <- map(data_neighbors, ~clust_vct[.x])

    ## nearest cluster neighbors -------

    clust_levs <- levels(clust_assignment$clust_id)

    clust_id <- NULL

    clust_assignment <- mutate(clust_assignment,
                               clust_id = as.character(clust_id))

    data_cluster_assignment <-
      split(clust_assignment,
            factor(clust_assignment$observation,
                   clust_assignment$observation))

    if(kNN_cluster > 1) {

      clust_labels <- rownames(as.matrix(clust_dist))

      cluster_neighbors <- dbscan::kNN(clust_dist, k = kNN_cluster)$id

      rownames(cluster_neighbors) <- clust_labels

      cluster_neighbors <-
        map(rownames(cluster_neighbors),
            ~clust_labels[cluster_neighbors[.x, ]])

      cluster_neighbors <- set_names(cluster_neighbors,
                                     clust_labels)

      data_cluster_assignment <-
        map(data_cluster_assignment,
            ~c(.x$clust_id[[1]], cluster_neighbors[[.x$clust_id[[1]]]]))

    } else {

      data_cluster_assignment <-
        map(data_cluster_assignment, ~.x$clust_id[[1]])

    }

    ## comparison of the neighborhoods -----

    ## mean fraction of k-nearest data points belonging to the k-nearest
    ## neighbor clusters

    kNN_fractions <-
      map2(data_neighbor_assignment,
           data_cluster_assignment,
           function(x, y) x %in% y)

    kNN_fractions <-
      tibble(observation = names(kNN_fractions),
             kNN_data = kNN_data,
             kNN_cluster = kNN_cluster,
             frac_np = map_dbl(kNN_fractions, mean))

    kNN_fractions <-
      left_join(kNN_fractions,
                clust_assignment[c('observation', 'clust_id')],
                by = 'observation')

    kNN_fractions <- mutate(kNN_fractions,
                            clust_id = factor(clust_id, clust_levs))

    knb(kNN_fractions[c('observation',
                        'clust_id',
                        'kNN_data',
                        'kNN_cluster',
                        'frac_np')])

  }

#' @rdname n_preservation

  np_reduction <- function(data_dist,
                           layout_dist,
                           kNN_data = 5) {

    ## input control ------

    stopifnot(inherits(data_dist, 'dist'))
    stopifnot(inherits(layout_dist, 'dist'))

    stopifnot(is.numeric(kNN_data))

    kNN_data <- as.integer(kNN_data)

    ## knearest neighbors in the data set and reduction layout -------

    dist_lst <- list(data = data_dist,
                     layout = layout_dist)

    point_labels <- map(dist_lst, ~rownames(as.matrix(.x)))

    neighbor_lst <- map(dist_lst, dbscan::kNN, k = kNN_data)

    for(i in names(neighbor_lst)) {

      data_points <- rownames(neighbor_lst[[i]]$id)

      neighbor_lst[[i]] <- map(data_points,
                               ~point_labels[[i]][neighbor_lst[[i]]$id[.x, ]])

      neighbor_lst[[i]] <- set_names(neighbor_lst[[i]],
                                     data_points)

    }

    ## preservation of the neighborhood --------

    ## fraction of nearest neighbors shared between the data set and the layout

    shared_neighbors <-
      map2(neighbor_lst[['data']], neighbor_lst[['layout']],
           intersect)

    knb(tibble(observation = names(shared_neighbors),
               kNN_data = kNN_data,
               kNN_cluster = NA,
               clust_id = NA,
               frac_np = map_dbl(shared_neighbors, length)/kNN_data))

  }

# Distance between observations and SOM nodes -------

#' Distance between observations and nodes of a self-organizing map.
#'
#' @description
#' Computes distances between observations and nodes of a self-organizing map
#' (SOM) with distance metrics extracted from a SOM clustering object.
#'
#' @return a numeric matrix with observations in rows and nodes in columns.
#' If `kNN` is specified a list of nearest neighbor nodes for observations.
#'
#' @param x a `clust_analysis` object.
#' @param kNN optional, number of the nearest neighbor nodes. If `NULL`, a numeric
#' matrix is returned. If provided, a list of nearest nodes named after
#' observations is returned.

  data_node_distance <- function(x, kNN = NULL) {

    stopifnot(is_clust_analysis(x))

    if(!x$clust_fun %in% c('som', 'supersom')) {

      warning("Data - node distances are implemented only for SOM.",
              call. = FALSE)

      return(NULL)

    }

    ## single layer case --------

    if(x$clust_fun == 'som') {

      dist_method <- x$dist_method

      data_mtx <- as.matrix(model.frame(x))

      node_mtx <- x$clust_obj$codes[[1]]

      rownames(node_mtx) <-
        stringi::stri_replace(rownames(node_mtx),
                              regex = '^V',
                              replacement = '')

      mix_mtx <- rbind(data_mtx, node_mtx)

      dists <- calculate_dist(mix_mtx, method = dist_method)

      dists <- dists[rownames(data_mtx), rownames(node_mtx)]


    } else {

      dist_methods <- x$clust_obj$dist.fcts

      wt <-
        x$clust_obj$user.weights *
        x$clust_obj$distance.weights

      data_layers <- map(model.frame(x), as.matrix)

      node_layers <- x$clust_obj$codes

      for(i in seq_along(node_layers)) {

        rownames(node_layers[[1]]) <-
          stringi::stri_replace(rownames(node_layers[[1]]),
                                regex = '^V',
                                replacement = '')

      }

      mix_mtx <- map2(data_layers, node_layers, rbind)

      dists <-
        map2(mix_mtx, dist_methods,
             ~calculate_dist(data = .x, method = .y))

      dists <- map2(dists, wt, `*`)

      dists <- reduce(dists, `+`)

      dists <-
        dists[rownames(data_layers[[1]]), rownames(node_layers[[1]])]

    }

    ## output -------

    if(is.null(kNN)) return(dists)

    neighbor_nodes <-
      map(rownames(dists), ~sort(dists[.x, ]))

    neighbor_nodes <- map(neighbor_nodes, ~names(.x[1:kNN]))

    set_names(neighbor_nodes, rownames(dists))

  }

# SOM grid neighbors ---------

#' Find neighbor nodes in the initial layout of a self-organizing map.
#'
#' @description
#' Identifies nearest nodes in the initial layout of a self-organizing map,
#' i.e. prior to fitting to the data points.
#'
#' @details
#' The grid architecture information is extracted from the `kohonen` object.
#' The function returns NULL and a warning is applied to a non-SOM clustering
#' analysis object.
#'
#' @param x a `clust_analysis` object.
#'
#' @return a list named after node identifiers with vectors of identifiers of
#' the nearest nodes.

  node_neighbors <- function(x) {

    ## entry control ------

    stopifnot(is_clust_analysis(x))

    if(!x$clust_fun %in% c('som', 'supersom')) {

      warning('Nearest nodes can be find only for SOM.',
              call. = FALSE)

      return(NULL)
    }

    ## node architecture information -------

    grid <- x$clust_obj$grid

    node_arch <- grid$pts

    toroid <- grid$toroidal

    xdim <- grid$xdim
    ydim <- grid$xdim

    topo <- grid$topo

    node_dists <- kohonen::unit.distances(grid)

    node_labels <- as.character(1:nrow(node_dists))

    rownames(node_dists) <- node_labels
    colnames(node_dists) <- node_labels

    node_dists <- map(node_labels, ~node_dists[.x, ])

    ## resorting to round, since the distances between neighbors
    ## never exactly one for hexagonal or toroidal grids

    node_dists <- map(node_dists, ~.x[round(.x) == 1])

    node_dists <- map(node_dists, names)

    set_names(node_dists, node_labels)

  }

# Topology error -------

#' Topology error of self-organizing maps.
#'
#' @description
#' Checks for topology errors for consecutive data points.
#'
#' @details
#' The procedure of topology error computation is as follows: for each
#' observation, two nearest self-organizing map (SOM) nodes are identified.
#' If such nodes are neighbors in the initial layout of SOM prior to data
#' fitting, correct topology (coded with 0) is returned and an error otherwise
#' (coded as 1).
#'
#' @return An object of the \code{\link{knb}} class with
#' \code{\link{summary.knb}} and \code{\link{plot.knb}} methods.
#'
#' @param x a `clust_analysis` object.

  check_topo <- function(x) {

    nearest_nodes <- data_node_distance(x, kNN = 2)

    neighbors <- node_neighbors(x)

    topo_res <-
      map_lgl(nearest_nodes,
              function(node_pair) !node_pair[2] %in% neighbors[[node_pair[1]]])

    topo_res <-
      tibble(observation = names(topo_res),
             topo_error = unname(as.numeric(topo_res)),
             kNN_data = NA,
             kNN_cluster = NA)

    topo_res <-
      left_join(topo_res,
                extract(x, 'assignment')[c('observation', 'clust_id')],
                by = 'observation')

    knb(topo_res[c('observation',
                   'clust_id',
                   'kNN_data',
                   'kNN_cluster',
                   'topo_error')])

  }

# END ------
