# Combined SOM + simple clustering functions.

#' Cluster self-organizing map nodes.
#'
#' @description Performs clustering of the self-orgnizing map (SOM) nodes with
#' one of the clustering functions provided by the clustTools package.
#' @param data a numeric data frame, matrix or a red_analysis object.
#' If a red_analysis object is provided, its component/score table will be
#' clustered.
#' @param distance_som metric of distance between the observations, used for SOM
#' development. See: \code{\link{get_kernel_info}}.
#' @param xdim x dimension of the SOM grid,
#' see: \code{\link[kohonen]{somgrid}} for details.
#' @param ydim y dimension of the SOM grid,
#' #' see: \code{\link[kohonen]{somgrid}} for details.
#' @param topo SOM grid topology, see: \code{\link[kohonen]{somgrid}}
#' for details. 'hexagonal' for default.
#' @param neighbourhood.fct neighborhood function, 'gaussian' for default.
#' @param toroidal logical, should toroidal grid be used?
#' @param rlen number of the SOM algorithm iterations.
#' @param node_clust_fun a function provided by the clustTools package used to
#' cluster the SOM nodes. Alternatively, a custom function may be provided here,
#' which returns a clust_analysis class object.
#' @param distance_nodes metric of distance between the nodes, used for SOM
#' development. See: \code{\link{get_kernel_info}}.
#' @param seed initial setting of the random number generator.
#' @param ... extra arguments passed to node_clust_fun, such as k number of
#' clusters.
#' @export

  combi_cluster <- function(data,
                            distance_som = 'euclidean',
                            xdim = 5,
                            ydim = 4,
                            topo = 'hexagonal',
                            neighbourhood.fct = 'gaussian',
                            toroidal = FALSE,
                            rlen = 500,
                            node_clust_fun = hcluster,
                            distance_nodes = 'euclidean',
                            seed = 1234, ...) {

    ## entry control

    if(all(class(data) == 'red_analysis')) {

      data <- tibble::column_to_rownames(data$component_tbl,
                                         'observation')

    }

    clustTools:::check_numeric(data)

    if(!distance_som %in% clustTools::get_kernel_info()) {

      stop('Invalid observation distance method.', call. = FALSE)

    }

    if(!distance_nodes %in% clustTools::get_kernel_info()) {

      stop('Invalid node distance method.', call. = FALSE)

    }

    xdim <- as.integer(xdim)
    ydim <- as.integer(ydim)
    rlen <- as.integer(rlen)

    topo <- match.arg(topo[1],
                      c('hexagonal', 'rectangular'))

    neighbourhood.fct <- match.arg(neighbourhood.fct[1],
                                   c('gaussian', 'bubble'))

    stopifnot(is.logical(toroidal))

    stopifnot(is.function(node_clust_fun))

    ## SOM clustering of the observations

    som_clust <- clustTools::som_cluster(data = data,
                                         distance_method = distance_som,
                                         xdim = xdim,
                                         ydim = ydim,
                                         topo = topo,
                                         neighbourhood.fct = neighbourhood.fct,
                                         toroidal = toroidal,
                                         rlen = rlen,
                                         seed = seed)

    ## clustering of the nodes

    node_clust <- node_clust_fun(data = som_clust$clust_obj$codes[[1]],
                                 distance_method = distance_nodes,
                                 seed = seed, ...)

    ## assignment to the nodes and node clusters

    som_ass <- som_clust$clust_assignment[c('observation', 'node')]

    node_ass <- rlang::set_names(node_clust$clust_assignment[c('observation',
                                                               'clust_id')],
                                 c('node', 'clust_id'))

    node_ass <- dplyr::mutate(node_ass,
                              node = stringi::stri_replace(node,
                                                           fixed = 'V',
                                                           replacement = ''))

    combi_ass <- dplyr::left_join(som_ass,
                                  node_ass,
                                  by = 'node')

    ## output

    clustTools::combi_analysis(list(clust_analyses = list(observation = som_clust,
                                                          node = node_clust),
                                    clust_assignment = combi_ass,
                                    dots = rlang::list2(...)))

  }

# END -----
