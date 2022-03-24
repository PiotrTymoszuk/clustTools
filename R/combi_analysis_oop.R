# S3 OOP methods for the combin_analysis class.

# General extractor function --------

#' Extract features of a combi_analysis object.
#'
#' @description A general extractor method for accessing properties and features
#' of a combi_analysis object.
#' @param x a combi_analysis object.
#' @inheritParams extract.clust_analysis
#' @return a list with the requested feature for the observation and node clustering
#' @export extract.combi_analysis
#' @export

  extract.combi_analysis <- function(x,
                                     type = c('distance',
                                              'assignment',
                                              'clust_object',
                                              'data',
                                              'object')) {

    stopifnot(all(class(x) == 'combi_analysis'))

    if(type == 'assignment') {

      x$clust_assignment

    } else {

      purrr::map(x$clust_analyses,
                 extract,
                 type = type)

    }

  }

# Appearance and class testing ------

#' Test for the combi_analysis class.
#'
#' @description Tests if an object belongs to the combi_analysis class.
#' @param x an object.
#' @return a logical value.
#' @export

  is_combi_analysis <- function(x) {

    all(class(x) == 'combi_analysis')

  }

#' Printing of a clust_analysis object.
#'
#' @description Prints a clust_analysis object.
#' @param x a combi_analysis object.
#' @param ... extra arguments, currently none.
#' @export

  print.combi_analysis <- function(x, ...) {

    stopifnot(all(class(x) == 'combi_analysis'))

    purrr::walk(x$clust_analyses, print)

  }

# Data, observation number and distance access ----

#' Access the clust_analysis data frame.
#'
#' @description Retrieves the data frame used for clustering analysis.
#' @param formula a combi_analysis object.
#' @param ... extra arguments, currently none.
#' @return a list of data frames.
#' @export

  model.frame.combi_analysis <- function(formula, ...) {

    stopifnot(all(class(formula) == 'combi_analysis'))

    purrr::map(formula$clust_analyses, model.frame)

  }

#' Access the distances of a combi_analysis object.
#'
#' @description Retrieves the matrices with distances between the observations
#' and SOM nodes.
#' @param x a combi_analysis object.
#' @param ... extra arguments, currently none.
#' @return a list of numeric matrices.
#' @export

  dist.combi_analysis <- function(x, ...) {

    stopifnot(all(class(x) == 'combi_analysis'))

    purrr::map(x$clust_analyses, dist)

  }

#' Number of observations of a combi_analysis object.
#'
#' @description Retrieves the number of observations and SOM nodes used to generate
#' a combi_analysis object.
#' @param object a combi_analysis object.
#' @param ... extra arguments, currently none.
#' @return a list with the numbers of observations and clustering features.
#' @export

  nobs.combi_analysis <- function(object, ...) {

    stopifnot(all(class(object) == 'combi_analysis'))

    purrr::map(object$clust_analyses, nobs)

  }

#' Numbers of observations assigned to the clusters.
#'
#' @description Retrieves the numbers of observations assigned to the clusters.
#' @param x a combi_analysis object.
#' @param ... extra arguments, currently none.
#' @return a list of data frames with the numbers of observations in
#' the clusters
#' @export

  ngroups.combi_analysis <- function(x, ...) {

    stopifnot(all(class(x) == 'combi_analysis'))

    purrr::map(x$clust_analyses, ngroups)

  }

# Dimensionality reduction ---------

#' Dimensionality reduction analysis of the combi_analysis object data
#' or distance matrix.
#'
#' @description Performs principal component analysis (PCA), multi-dimensional
#' scaling (MDS) or uniform manifold approximation and projection (UMAP) of the
#' combi_analysis object data or distance matrix.
#' The analysis is done for the global clustering.
#' @details See \code{\link{reduce_data}} for the implementation details.
#' The distance methods, relevant for MDS and UMAP, are taken over from the
#' combi_object. Hence, some distances may crash the analysis with UMAP, see:
#' \code{\link[umap]{umap.defaults}} for the compatible distances.
#' @param x a combi_analysis object.
#' @param kdim number of dimensions. If NULL, kdim is set to the number of
#' clusters.
#' @param red_fun reduction analysis function: 'pca' (PCA), 'mds' (MDS) or
#' 'umap' (UMAP).
#' @param with type of the input data for the reduction analysis:
#' the clustering data ('data') or the matrix of distances ('distance').
#' @param ... extra arguments passed to \code{\link{reduce_data}}.
#' @return a red_analysis object with the component/score table
#' containing the cluster assignment information ('clust_id' variable).
#' @export components.combi_analysis
#' @export

  components.combi_analysis <- function(x,
                                        kdim = NULL,
                                        red_fun = c('pca', 'mds', 'umap'),
                                        with = c('distance', 'data'), ...) {

    ## entry control

    stopifnot(all(class(x) == 'combi_analysis'))

    with <- match.arg(with[1], c('distance', 'data'))

    red_fun <- match.arg(red_fun[1], c('pca', 'mds', 'umap'))

    ## reduction analysis

    red_analysis <- components(x$clust_analyses$observation,
                               kdim = kdim,
                               red_fun = red_fun,
                               with = with, ...)

    red_analysis$component_tbl <- dplyr::select(red_analysis$component_tbl,
                                                observation,
                                                dplyr::starts_with('comp'))


    red_analysis$component_tbl <- dplyr::left_join(red_analysis$component_tbl,
                                                   x$clust_assignment,
                                                   by = 'observation')

    red_analysis

  }

# Clustering variance -------

#' Calculate clustering variance.
#'
#' @description Calculates the clustering sum of squares (total, within
#' clusters, total within clusters and between clusters) as well as the
#' fraction of 'explained' clustering variance. The later is the ratio of
#' between-cluster sum of squares to the total sum of squares statistic.
#' @param x a combi_analysis object.
#' @param ... extra arguments, currently none.
#' @export

  var.combi_analysis <- function(x, ...) {

    stopifnot(all(class(x) == 'combi_analysis'))

    get_sum_sq(dist_mtx = dist.combi_analysis(x, 'distance')[[1]],
               assignment = x$clust_assignment)

  }

# Plotting --------

#' Plot selected features of a combi_analysis object.
#'
#' @description The plotting method for the clust_analysis class. Enables
#' plotting of the standard diagnostic plots used for the optimal cluster number
#' determination for the node clustering (dendrogram, WSS- and
#' silhouette curve), results of the reduction analysis,
#' heat map of the distances between the observations and SOM nodes as
#' well as the self-organizing map training process. It is also possible to plot
#' the first two variables of the clustering data frame, an option which is
#' attractive, if the clustering of reduction analysis was performed.
#' @param x a combi_analysis object.
#' @param type the type of plots:
#' 'diagnostic' returns a series of diagnostic plots for the SOM construction
#' and clustering of the nodes;
#' 'components' plots the results of reduction analysis done with the clustering
#' data or the distance matrix (see: \code{\link{components.clust_analysis}});
#' 'heat_map' plots the distances between observations and nodes as a heat maps,
#' 'training' plots the mean distance to the SOM winning unit as a function
#' of the iteration number; 'data' works only if reduction analysis results
#' were used for clustering and plots the first two comonents/dimensions.
#' @inheritParams plot.clust_analysis
#' @return a ggplot object or a list of ggplot objects, as specified by the
#' 'type' argument
#' @export plot.combi_analysis
#' @export

  plot.combi_analysis <- function(x,
                                  type = c('diagnostic',
                                           'components',
                                           'heat_map',
                                           'training',
                                           'data'),
                                  cust_theme = ggplot2::theme_classic(),
                                  jitter_width = 0,
                                  jitter_height = 0,
                                  point_alpha = 1, ...) {

    ## entry control

    stopifnot(all(class(x) == 'combi_analysis'))

    type <- match.arg(type[1],
                      c('diagnostic',
                        'components',
                        'heat_map',
                        'training',
                        'data'))

    stopifnot(any(class(cust_theme) == 'theme'))

    ## plotting

    plot_list <- purrr::map(x$clust_analyses, plot.clust_analysis,
                            type = type,
                            cust_theme = cust_theme,
                            jitter_height = jitter_height,
                            jitter_width = jitter_width, ...)

    ## summary component plots for the final clustering results

    if(type == 'components') {

      obs_red <- components(x$clust_analyses$observation, ...)

      score_tbl <- dplyr::select(obs_red$component_tbl,
                                 observation,
                                 dplyr::starts_with('comp'))

      score_tbl <- dplyr::left_join(score_tbl,
                                    x$clust_assignment,
                                    by = 'observation')

      sdevs <- var(obs_red)

      if(obs_red$red_fun == 'pca') {

        ax_labs <- purrr::map2(c('PC1', 'PC2'),
                               signif(sdevs$perc_var[1:2], 3),
                               ~paste0(.x, ', ', .y, '%'))

      } else {

        ax_labs <- purrr::map2(c('Dim 1', 'Dim 2'),
                               signif(sdevs$perc_var[1:2], 3),
                               ~paste0(.x, ', ', .y, '%'))

      }

      plot_list$final <- clustTools::plot_point(data = score_tbl,
                                                x_var = 'comp_1',
                                                y_var = 'comp_2',
                                                fill_var = 'clust_id',
                                                plot_title = switch(obs_red$red_fun ,
                                                                    pca = 'PCA',
                                                                    mds = 'MDS',
                                                                    umap = 'UMAP'),
                                                plot_tag = plot_list$observation$labels$tag,
                                                x_lab = ax_labs[[1]],
                                                y_lab = ax_labs[[2]],
                                                cust_theme = cust_theme,
                                                jitter_height = jitter_height,
                                                jitter_width = jitter_width,
                                                fill_lab = 'Cluster ID',
                                                point_alpha = point_alpha)

    }

    plot_list

  }

# Semi-supervised clustering ------

#' Semi-supervised clustering.
#'
#' @description Projects the cluster assignment onto new data using simple
#' observation matching or a k-nearest neighbor (kNN) label propagation
#' algorithm.
#' @details For the implementation details, see: \code{\link{propagate}}.
#' The default distance metric is extracted from the combi_analysis object.
#' The cluster projection is done on the top level, i.e. takes into account the
#' final assignment of the observations to the clusters and ignoring
#' the SOM nodes.
#' @param object a combi_analysis object.
#' @inheritParams predict.clust_analysis
#' @return a clust_analysis object.
#' @export predict.combi_analysis
#' @export

  predict.combi_analysis <- function(object,
                                     newdata = NULL,
                                     type = c('class', 'propagation'), ...) {

    ## entry control

    stopifnot(all(class(object) == 'combi_analysis'))

    type <- match.arg(type[1],
                      c('class', 'propagation'))

    if(is.null(newdata)) {

      return(object$clust_assignment)

    }

    if(all(class(newdata) == 'red_analysis')) {

      newdata <- tibble::column_to_rownames(newdata$component_tbl,
                                            'observation')

    }

    clustTools:::check_numeric(newdata)

    ## prediction

    if(type == 'class') {

      train_assignment <- tibble::column_to_rownames(object$clust_assignment,
                                                     'observation')

      if(nrow(newdata) != nrow(train_assignment)) {

        stop('The numbers of rows in new data and the table used for cluster development must be equal',
             call. = FALSE)

      }

      newdata <- as.data.frame(newdata)

      if(!is.null(rownames(newdata))) {

        test_assignment <- tibble::tibble(observation = rownames(newdata),
                                          clust_id = train_assignment[rownames(newdata),
                                                                      'clust_id'])

      } else {

        warning('Unnamed observations in new data')

        test_assignment <- tibble::rownames_to_column(train_assignment,
                                                      'observation')

        train_assignment <- tibble::as_tibble(train_assignment)

      }

      ## output

      model_frame <- rlang::enexpr(newdata)

      clustTools::clust_analysis(list(data = rlang::quo(!!model_frame),
                                      dist_mtx = clustTools::calculate_dist(newdata,
                                                                            method = object$clust_analyses$observation$dist_method),
                                      dist_method = object$clust_analyses$observation$dist_method,
                                      clust_fun = 'prediction',
                                      clust_obj = NULL,
                                      clust_assignment = test_assignment))

    } else {

      clustTools::propagate(object = object,
                            newdata = newdata, ...)

    }

  }

# Cross-validation ------

#' Cross-validate the clustering object.
#'
#' @description Checks the stability of a clustering solution by
#' cross-validation (CV) and the classification error as a measure of the cluster
#' stability.
#' @details By principle similar to cross-validation of any machine learning
#' multi-level classifier. The training portion of a CV split is used to develop
#' of a cluster structure and the projection on the test portion is accomplished
#' by k-nearest neighbor (kNN) label propagation algorithm. For its
#' implementation details, see: \code{\link{propagate}}.
#' The fold are generated with \code{\link[caret]{createFolds}}.
#' @param x a combi_analysis object.
#' @inheritParams cv.clust_analysis
#' @return a list containing the global clust_analysis object, projection
#' (prediction) results and prediction summary for each fold and a prediction
#' summary for the whole CV.
#' @export cv.combi_analysis
#' @export

  cv.combi_analysis <- function(x,
                                nfolds = 5,
                                kNN = 5,
                                simple_vote = TRUE,
                                resolve_ties = FALSE,
                                kernel_fun = function(x) 1/x,
                                seed = 1234,
                                .parallel = FALSE) {

    ## entry control

    stopifnot(all(class(x) == 'combi_analysis'))

    nfolds <- as.integer(nfolds)

    kNN <- as.integer(kNN)

    stopifnot(is.logical(simple_vote))
    stopifnot(is.logical(resolve_ties))
    stopifnot(is.logical(.parallel))

    stopifnot(is.function(kernel_fun))

    ## cross-validation

      args <- rlang::list2(data = rlang::eval_tidy(x$clust_analyses$observation$data),
                           nfolds = 10,
                           kNN = 5,
                           clustering_fun = combi_cluster,
                           distance_som = x$clust_analyses$observation$dist_method,
                           xdim = x$clust_analyses$observation$grid$xdim,
                           ydim = x$clust_analyses$observation$grid$ydim,
                           topo = x$clust_analyses$observation$grid$topo,
                           neighbourhood.fct = as.character(x$clust_analyses$observation$grid$neighbourhood.fct),
                           toroidal = x$clust_analyses$observation$grid$toroidal,
                           rlen = nrow(x$clust_analyses$observation$clust_obj$changes),
                           node_clust_fun = switch(x$clust_analyses$node$clust_fun,
                                                   hclust = hcluster,
                                                   kmeans = kcluster,
                                                   pam = kcluster,
                                                   dbscan = dbscan_cluster,
                                                   som = som_cluster),
                           kernel_fun = kernel_fun,
                           simple_vote = simple_vote,
                           resolve_ties = resolve_ties,
                           seed = seed,
                           .parallel = .parallel)

    args <- c(args, x$dots)

    cv_call <- rlang::call2('cv_cluster',
                            !!!args)

    eval(cv_call)

  }

# Variable importance -------

#' Determine clustering feature importance.
#'
#' @description Determines importance of specific clustering variables by
#' comparing the fraction of 'explained' clustering variance of the input
#' clustering object and the object generated with the variable
#' re-shuffled randomly.
#' @param x a combi_analysis object.
#' @inheritParams impact.clust_analysis
#' @return a data frame of class 'importance'
#' @export impact.combi_analysis
#' @export

  impact.combi_analysis <- function(x,
                                    seed = 1234,
                                    .parallel = FALSE) {

    ## entry control

    stopifnot(all(class(x) == 'combi_analysis'))

    ## common parameters

    args <- rlang::list2(data = rlang::eval_tidy(x$clust_analyses$observation$data),
                         clustering_fun = combi_cluster,
                         distance_som = x$clust_analyses$observation$dist_method,
                         xdim = x$clust_analyses$observation$grid$xdim,
                         ydim = x$clust_analyses$observation$grid$ydim,
                         topo = x$clust_analyses$observation$grid$topo,
                         neighbourhood.fct = as.character(x$clust_analyses$observation$grid$neighbourhood.fct),
                         toroidal = x$clust_analyses$observation$grid$toroidal,
                         rlen = nrow(x$clust_analyses$observation$clust_obj$changes),
                         node_clust_fun = switch(x$clust_analyses$node$clust_fun,
                                                 hclust = hcluster,
                                                 kmeans = kcluster,
                                                 pam = kcluster,
                                                 dbscan = dbscan_cluster,
                                                 som = som_cluster),
                         seed = seed,
                         .parallel = .parallel)

    args <- c(args, x$dots)

    imp_call <- rlang::call2('importance_cluster',
                             !!!args)

    eval(imp_call)

  }

# END ------
