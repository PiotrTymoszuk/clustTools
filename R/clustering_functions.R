# Functions of simple clustering of data frames or reduction analysis results.

# Investigation of the clustering tendency ------

#' Check clustering tendency of a data set.
#'
#' @description Check clustering tendency of a data set as compared with a
#' random data set using Hopkins statistic.
#' @details The p value for the Hopkins statistic is calculated based on the
#' beta distribution of its values. Technically, the function is an enriched
#' wrapper around \code{\link[factoextra]{get_clust_tendency}}.
#' @param data a data frame, tibble or a matrix. Numeric variables only.
#' @param n the number of points selected from sample space which is also the
#' number of points selected from the given sample (data).
#' @param seed initial setting of the random number generator.
#' @param ... extra arguments passed to
#' \code{\link[factoextra]{get_clust_tendency}}.
#' @return The values of the Hopkins statistic, p value and a heat map plot.
#' @export

  get_clust_tendency <- function(data, n, seed = 1234, ...) {

    ## entry control

    clustTools:::check_numeric(data)

    ## computation

    tend <- factoextra::get_clust_tendency(data = data,
                                           n = n,
                                           seed = seed, ...)

    tend$p_value <- 1 - stats::pbeta(tend$hopkins_stat,
                                     shape1 = n,
                                     shape2 = n)

    return(tend)

  }

# Hierarchical clustering ------

#' Hierarchical clustering.
#'
#' @description Performs hierarchical clustering analysis of a numeric data
#' frame, matrix or the results of a reduction analysis.
#' @details Technically, a wrapper around \code{\link[stats]{hclust}}. If a
#' red_analysis object is provided as the data argument, the observation
#' component/score table is subjected to clustering.
#' @param data a numeric data frame or matrix or a red_analysis object.
#' @param distance_method name of the distance metric, see:
#' \code{\link{get_kernel_info}}.
#' @param k number of clusters.
#' @param hc_method the hierarchical clustering algorithm, see:
#' \code{\link[stats]{hclust}} for details.
#' @param seed initial setting of the random number generator.
#' @param ... extra arguments passed to \code{\link[stats]{hclust}}.
#' @return a clust_analysis object.
#' @export

  hcluster <- function(data,
                       distance_method = 'euclidean',
                       k = 2,
                       hc_method = 'ward.D2',
                       seed = 1234, ...) {

    ## entry control

    if(all(class(data) == 'red_analysis')) {

      data <- tibble::column_to_rownames(data$component_tbl,
                                         'observation')

    }

    clustTools:::check_numeric(data)

    if(!distance_method %in% clustTools::get_kernel_info()) {

      stop('Invalid distance method.', call. = FALSE)

    }

    k <- as.integer(k)

    ## RNG

    set.seed(seed = seed)

    ## distance calculation

    dist_mtx <- clustTools::calculate_dist(data = data,
                                           method = distance_method)

    ## hierarchical clustering and cluster assignment table

    hclust_str <- stats::hclust(as.dist(dist_mtx),
                                method = hc_method, ...)

    hclust_ass <- stats::cutree(hclust_str, k = k)

    hclust_ass <- tibble::tibble(observation = names(hclust_ass),
                                 clust_id = factor(unname(hclust_ass)))

    ## output

    model_frame <- rlang::enexpr(data)

    clustTools::clust_analysis(list(data = rlang::quo(!!model_frame),
                                    dist_mtx = dist_mtx,
                                    dist_method = distance_method,
                                    clust_fun = 'hclust',
                                    clust_obj = hclust_str,
                                    clust_assignment = hclust_ass,
                                    hc_method = hc_method,
                                    dots = rlang::list2(...)))

  }

# K-means or medoid clustering -----

#' K-means or medoid clustering.
#'
#' @description Performs k-means PAM (partition around medoids) clustering
#' analysis of of a numeric data frame, matrix or the results of
#' a reduction analysis.
#' @details Technically, a wrapper around \code{\link[stats]{kmeans}} and
#' \code{\link[cluster]{pam}}. If a red_analysis object is provided as the
#' data argument, the observation component/score table is subjected to
#' clustering.
#' @param data a numeric data frame or matrix or a red_analysis object.
#' @param distance_method name of the distance metric, see:
#' \code{\link{get_kernel_info}}.
#' @param clust_fun the name of the clustering function, currently implemented
#' are 'kmeans' and 'pam'.
#' @param k number of clusters.
#' @param seed initial setting of the random number generator.
#' @param ... extra arguments passed to \code{\link[stats]{kmeans}} or
#' \code{\link[cluster]{pam}}.
#' @return a clust_analysis object.
#' @export

  kcluster <- function(data,
                       distance_method = 'euclidean',
                       clust_fun = c('kmeans', 'pam'),
                       k = 2,
                       seed = 1234, ...) {

    ## entry control

    if(all(class(data) == 'red_analysis')) {

      data <- tibble::column_to_rownames(data$component_tbl,
                                         'observation')

    }

    clustTools:::check_numeric(data)

    if(!distance_method %in% clustTools::get_kernel_info()) {

      stop('Invalid distance method.', call. = FALSE)

    }

    k <- as.integer(k)

    clust_fun <- match.arg(clust_fun[1], c('kmeans', 'pam'))

    fun <- switch(clust_fun,
                  kmeans = stats::kmeans,
                  pam= cluster::pam)

    ## RNG

    set.seed(seed = seed)

    ## distance calculation

    dist_mtx <- as.dist(clustTools::calculate_dist(data,
                                                   method = distance_method))

    ## kmeans/pam clustering and cluster assignment table

    kclust_str <- fun(dist_mtx, k, ...)

    kclust_ass <- tibble::tibble(observation = names(kclust_str$cluster),
                                 clust_id = factor(unname(kclust_str$cluster)))

    ## output

    model_frame <- rlang::enexpr(data)

    clustTools::clust_analysis(list(data = rlang::quo(!!model_frame),
                                    dist_mtx = as.matrix(dist_mtx),
                                    dist_method = distance_method,
                                    clust_fun = clust_fun,
                                    clust_obj = kclust_str,
                                    clust_assignment = kclust_ass,
                                    dots = rlang::list2(...)))

  }

# Density clustering -------

#' Density clustering with DBSCAN.
#'
#' @description Performs DBSCAN clustering analysis of a numeric data frame,
#' matrix or the results of a reduction analysis.
#' @details Technically, a wrapper around \code{\link[dbscan]{dbscan}}. If a
#' red_analysis object is provided as the data argument, the observation
#' component/score table is subjected to clustering.
#' @param data a numeric data frame or matrix or a red_analysis object.
#' @param distance_method name of the distance metric, see:
#' \code{\link{get_kernel_info}}.
#' @param eps size (radius) of the epsilon neighborhood.
#' @param minPts number of minimum points required in the eps neighborhood for
#' core points (including the point itself).
#' @param seed initial setting of the random number generator.
#' @param ... extra arguments passed to \code{\link[dbscan]{dbscan}}.
#' @return a clust_analysis object.
#' @export

  dbscan_cluster <- function(data,
                             distance_method = 'euclidean',
                             eps,
                             minPts = 5,
                             seed = 1234, ...) {

    ## entry control

    if(all(class(data) == 'red_analysis')) {

      data <- tibble::column_to_rownames(data$component_tbl,
                                         'observation')

    }

    clustTools:::check_numeric(data)

    if(!distance_method %in% clustTools::get_kernel_info()) {

      stop('Invalid distance method.', call. = FALSE)

    }

    minPts <- as.integer(minPts)


    ## RNG

    set.seed(seed)

    ## distance calculation

    dist_mtx <- clustTools::calculate_dist(data,
                                           method = distance_method)

    ## dbscan clustering and cluster assignment table

    dclust_str <- dbscan::dbscan(as.dist(dist_mtx),
                                 eps = eps,
                                 minPts = minPts, ...)

    dclust_ass <- tibble::tibble(observation = rownames(dist_mtx),
                                 clust_id = factor(dclust_str$cluster))

    ## output

    model_frame <- rlang::enexpr(data)

    clustTools::clust_analysis(list(data = rlang::quo(!!model_frame),
                                    dist_mtx = dist_mtx,
                                    dist_method = distance_method,
                                    clust_fun = 'dbscan',
                                    clust_obj = dclust_str,
                                    clust_assignment = dclust_ass,
                                    eps = eps,
                                    minPts = minPts,
                                    dots = rlang::list2(...)))

  }

# Self-organizing maps ------

#' Self-organizing maps.
#'
#' @description Performs self-organizing map (SOM) clustering of a numeric data
#' frame, matrix or the results of a reduction analysis.
#' @details Technically, a wrapper around \code{\link[kohonen]{som}}. If a
#' red_analysis object is provided as the data argument, the observation
#' component/score table is subjected to clustering. Note, in order to make use
#' of the full set of distance measures, the package 'somKernels' need to be
#' installed and loaded.
#' @param data a numeric data frame or matrix or a red_analysis object.
#' @param distance_method name of the distance metric, see:
#' \code{\link{get_kernel_info}}.
#' @param xdim x dimension of the SOM grid,
#' see: \code{\link[kohonen]{somgrid}} for details.
#' @param ydim y dimension of the SOM grid,
#' #' see: \code{\link[kohonen]{somgrid}} for details.
#' @param topo SOM grid topology, see: \code{\link[kohonen]{somgrid}}
#' for details. 'hexagonal' for default.
#' @param neighbourhood.fct neighborhood function, 'gaussian' for default.
#' @param toroidal logical, should toroidal grid be used?
#' @param seed initial setting of the random number generator.
#' @param ... extra arguments passed to \code{\link[kohonen]{som}}.
#' @export

  som_cluster <- function(data,
                          distance_method = 'euclidean',
                          xdim = 5,
                          ydim = 4,
                          topo = c('hexagonal', 'rectangular'),
                          neighbourhood.fct = c('gaussian', 'bubble'),
                          toroidal = FALSE,
                          seed = 1234, ...) {

    ## entry control

    if(all(class(data) == 'red_analysis')) {

      data <- tibble::column_to_rownames(data$component_tbl,
                                         'observation')

    }

    clustTools:::check_numeric(data)

    if(!distance_method %in% clustTools::get_kernel_info()) {

      stop('Invalid distance method.', call. = FALSE)

    }

    xdim <- as.integer(xdim)
    ydim <- as.integer(ydim)

    topo <- match.arg(topo[1],
                      c('hexagonal', 'rectangular'))

    neighbourhood.fct <- match.arg(neighbourhood.fct[1],
                                   c('gaussian', 'bubble'))

    stopifnot(is.logical(toroidal))

    ## grid

    som_grid <- kohonen::somgrid(xdim = xdim,
                                 ydim = ydim,
                                 topo = topo,
                                 neighbourhood.fct = neighbourhood.fct,
                                 toroidal = toroidal)

    ## RNG

    set.seed(seed = seed)

    ## fitting

    kohonen_obj <- kohonen::som(as.matrix(data),
                                grid = som_grid,
                                dist.fcts = distance_method, ...)

    node_ass <- tibble::tibble(observation = rownames(kohonen_obj$data[[1]]),
                               clust_id = factor(kohonen_obj$unit.classif),
                               node = factor(kohonen_obj$unit.classif),
                               neuro_dist = kohonen_obj$distances)

    ## output

    model_frame <- rlang::enexpr(data)

    clustTools::clust_analysis(list(data = rlang::quo(!!model_frame),
                                    dist_mtx = calculate_dist(data = data,
                                                              method = distance_method),
                                    dist_method = distance_method,
                                    clust_fun = 'som',
                                    clust_obj = kohonen_obj,
                                    clust_assignment = node_ass,
                                    grid = som_grid,
                                    dots = rlang::list2(...)))

  }

# END -----
