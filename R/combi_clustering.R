# Combined SOM + simple clustering functions.

#' Cluster self-organizing map nodes.
#'
#' @description
#' Performs clustering of the self-organizing map (SOM) with
#' one of the clustering functions provided by the clustTools package.
#'
#' @details
#' The clustering procedure involves construction of SOM with the user-provided
#' data followed by unsupervised clustering of the inter-node distance matrix.
#' For `combi_cluster()` tackling single-layer SOM, the user is allowed to
#' specify distances both for the input data and the nodes.
#' In case when both distance methods are the same, the inter-node
#' distance matrix corresponds to a classical U-matrix as computed by
#' \code{\link[kohonen]{object.distances}} - this is also the
#' recommended default option.
#' `multi_cluster()`, which takes a multi-layer data set as the `data` argument,
#' the U matrix constructed by \code{\link[kohonen]{object.distances}} is
#' always used as the inter-node distance matrix subjected to unsupervised
#' clustering.
#'
#' @references
#' Vesanto J, Alhoniemi E. Clustering of the self-organizing map.
#' IEEE Trans Neural Networks (2000) 11:586–600. doi:10.1109/72.846731
#'
#' @param data for `combi_clust()`, a numeric data frame, matrix or
#' a `red_analysis` object. If a `red_analysis` object is provided, its
#' component/score table will be clustered. For `multi_clust()` a list
#' of such objects.
#' @param distance_som metric of distance between the observations, used for SOM
#' development. See: \code{\link{get_kernel_info}}.
#' @param distance_method a vector of distance names, that matches elemnts of
#' the `data` list.
#' @param xdim x dimension of the SOM grid,
#' see: \code{\link[kohonen]{somgrid}} for details.
#' @param ydim y dimension of the SOM grid,
#' #' see: \code{\link[kohonen]{somgrid}} for details.
#' @param topo SOM grid topology, see: \code{\link[kohonen]{somgrid}}
#' for details. 'hexagonal' for default.
#' @param neighbourhood.fct neighborhood function, 'gaussian' for default.
#' @param toroidal logical, should toroidal grid be used?
#' @param rlen number of the SOM algorithm iterations.
#' @param som_args a list of extra arguments passed to
#' \code{\link{som_cluster}}, \code{\link[kohonen]{som}} or
#' \code{\link[kohonen]{supersom}}. They may include weights for data layers
#' or the learning rate.
#' @param node_clust_fun a function provided by the clustTools package used to
#' cluster the SOM nodes. An alternative for `combi_cluster()`: a user-provided
#' function that takes a numeric data frame or matrix and returns a
#' `clust_analysis` object. An alternative for `multi_cluster()`: a
#' user-provided function that takes a dissimilarity object (R's `dist` class)
#' and returns a `clust_analysis` object.
#' @param distance_nodes metric of distance between the nodes, used for SOM
#' development. Defaults to `distance_som`. See: \code{\link{get_kernel_info}}.
#' @param seed initial setting of the random number generator.
#' @param ... extra arguments. For `combi_clust()`, they are passed to
#' `node_clust_fun` and may include e.g. `k` number of clusters.
#'
#' @return an object of the class \code{\link{combi_analysis}}.
#'
#' @export

  combi_cluster <- function(data,
                            distance_som = 'euclidean',
                            xdim = 5,
                            ydim = 4,
                            topo = 'hexagonal',
                            neighbourhood.fct = 'gaussian',
                            toroidal = FALSE,
                            rlen = 500,
                            som_args = NULL,
                            node_clust_fun = hcluster,
                            distance_nodes = distance_som,
                            seed = 1234, ...) {

    ## entry control -------

    if(is.matrix(data) | is.data.frame(data)) {

      if(is_red_analysis(data)) {

        data <- column_to_rownames(data$component_tbl,
                                   'observation')

      }

      check_numeric(data)

    } else if(is.list(data)) {

      stop(paste("Your 'data' seems to be a list of data sets. Consider",
                 "using 'multi_cluster' which tackles multi-layer data."),
           call. = FALSE)

    } else {

      stop("Unsupported 'data' format.", call. = FALSE)

    }

    if(!distance_som %in% get_kernel_info()) {

      stop('Invalid observation distance method.', call. = FALSE)

    }

    rlen <- as.integer(rlen)
    stopifnot(is.function(node_clust_fun))

    ## the remaining input checks are done by 'som_cluster()'

    ## SOM clustering of the observations ------

    som_call <-
      rlang::call2(.fn = 'som_cluster',
                   data = data,
                   distance_method = distance_som,
                   xdim = xdim,
                   ydim = ydim,
                   topo = topo,
                   neighbourhood.fct = neighbourhood.fct,
                   toroidal = toroidal,
                   rlen = rlen,
                   seed = seed,
                   !!!som_args)

    som_clust <- eval(som_call)

    ## clustering of the nodes -------

    node_clust <- node_clust_fun(data = som_clust$clust_obj$codes[[1]],
                                 distance_method = distance_nodes,
                                 seed = seed, ...)

    ## assignment to the nodes and node clusters -------

    som_ass <- som_clust$clust_assignment[c('observation', 'node')]

    node_ass <- set_names(node_clust$clust_assignment[c('observation',
                                                        'clust_id')],
                          c('node', 'clust_id'))

    node <- NULL

    node_ass <- mutate(node_ass,
                       node = stringi::stri_replace(node,
                                                    fixed = 'V',
                                                    replacement = ''))

    combi_ass <- left_join(som_ass,
                           node_ass,
                           by = 'node')

    ## output -------

    dots <- list2(...)

    ## special handling of lambda vectors for regularized
    ## cluster functions of the nodes

    dots$lambdas <- NULL

    if(!is.null(som_args)) {

      dots <- c(dots,
                list(som_args = som_args))

    }

    dots <- compact(dots)

    combi_analysis(list(clust_analyses = list(observation = som_clust,
                                              node = node_clust),
                        clust_assignment = combi_ass,
                        dots = dots))

  }

#' @rdname combi_cluster
#' @export

  multi_cluster <- function(data,
                            distance_method = 'euclidean',
                            xdim = 5,
                            ydim = 4,
                            topo = 'hexagonal',
                            neighbourhood.fct = 'gaussian',
                            toroidal = FALSE,
                            rlen = 500,
                            som_args = NULL,
                            node_clust_fun = hcluster,
                            seed = 1234, ...) {

    ## entry control for the data ------

    if(is.matrix(data) | is.data.frame(data)) {

      stop(paste("The function does not handle single data frames or matrices.",
                 "Consider using 'combi_cluster()' instead."),
           call. = FALSE)

    } else if(is.list(data)) {

      data <-
        map(data,
            function(x) if(is_red_analysis(x)) column_to_rownames(data$component_tbl,
                                                                  'observation') else x)

      purrr::walk(data, check_numeric)

    } else {

      stop('Unsupported data format.', call. = FALSE)

    }

    if(is.null(names(data))) {

      data <- set_names(data, paste0('layer_', 1:length(data)))

    }

    data <- map(data, as.matrix)

    ## distance - data consistency -------

    num_layers <- length(data)

    if(length(distance_method) != 1) {

      if(length(distance_method) != num_layers) {

        stop(paste("'distance_method' has to be a single name of a",
                   "distance or a vector of distance names."),
             call. = FALSE)

      }

    }

    av_distances <- get_kernel_info()

    dist_check <- map_lgl(distance_method, ~.x %in% av_distances)

    if(any(!dist_check)) {

      stop('At least one of distances in no supported.', call. = FALSE)

    }

    rlen <- as.integer(rlen)
    stopifnot(is.function(node_clust_fun))

    ## the remaining input checks are done by 'som_cluster()'

    ## SOM clustering of the observations ------

    som_call <-
      rlang::call2(.fn = 'som_cluster',
                   data = data,
                   distance_method = distance_method,
                   xdim = xdim,
                   ydim = ydim,
                   topo = topo,
                   neighbourhood.fct = neighbourhood.fct,
                   toroidal = toroidal,
                   rlen = rlen,
                   seed = seed,
                   !!!som_args)

    som_clust <- eval(som_call)

    ## clustering of the U matrix -------

    node_clust <- node_clust_fun(data = dist(som_clust, 'umatrix'),
                                 seed = seed, ...)

    ## node and cluster assignment --------

    som_ass <- som_clust$clust_assignment[c('observation', 'node')]

    node_ass <- set_names(node_clust$clust_assignment[c('observation',
                                                        'clust_id')],
                          c('node', 'clust_id'))

    node <- NULL

    node_ass <- mutate(node_ass,
                       node = stringi::stri_replace(node,
                                                    fixed = 'V',
                                                    replacement = ''))

    combi_ass <- left_join(som_ass,
                           node_ass,
                           by = 'node')

    ## output -------

    dots <- list2(...)

    if(!is.null(som_args)) {

      dots <- c(dots,
                list(som_args = som_args))

    }

    combi_analysis(list(clust_analyses = list(observation = som_clust,
                                              node = node_clust),
                        clust_assignment = combi_ass,
                        dots = dots))

  }

# END -----
