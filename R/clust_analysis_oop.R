# S3 OOP for the clust_analysis class.

# General extractor function -----

#' Extract features of a clust_analysis object.
#'
#' @description A general extractor method for accessing properties and features
#' of a clust_analysis object.
#' @param x a clust_analysis object.
#' @param type the feature name:
#' 'distance' extracts the matrix with distances between the observations,
#' 'data' the data set used for the clust_object generation,
#' 'assignment' assignment of the observations to the clusters,
#' 'clust_object' or 'object' returns the clustering object.
#' @param ... extra arguments, currently none.
#' @return the requested feature/property.
#' @export extract.clust_analysis
#' @export

  extract.clust_analysis <- function(x,
                                     type = c('distance',
                                              'assignment',
                                              'clust_object',
                                              'data',
                                              'object'), ...) {

    ## entry control

    stopifnot(all(class(x) == 'clust_analysis'))

    type <- match.arg(type[1],
                      c('distance',
                        'assignment',
                        'clust_object',
                        'data',
                        'object'))

    ## output

    switch(type,
           distance = x$dist_mtx,
           assignment = x$clust_assignment,
           clust_object = x$clust_obj,
           data = rlang::eval_tidy(x$data),
           object = x$clust_obj)

  }

# Appearance and class testing -----

#' Test for the clust_analysis class.
#'
#' @description Tests if an object belongs to the clust_analysis class.
#' @param x an object.
#' @return a logical value.
#' @export

  is_clust_analysis <- function(x) {

    all(class(x) == 'clust_analysis')

  }

#' Printing of a clust_analysis object.
#'
#' @description Prints a clust_analysis object.
#' @param x a clust_analysis object.
#' @param ... extra arguments, currently none.
#' @export

  print.clust_analysis <- function(x, ...) {

    stopifnot(all(class(x) == 'clust_analysis'))

    cat(paste0('Clustering analysis with ',
               toupper(x$clust_fun), ' and ',
               x$dist_method, ' distance method.'))

    cat('\nCluster assignment:\n')

    print(tibble::as_tibble(x$clust_assignment))


  }

# Data, observation number and distance access ----

#' Access the combi_analysis data frames.
#'
#' @description Retrieves the data frame used for clustering analysis.
#' @param formula a combi_analysis object.
#' @param ... extra arguments, currently none.
#' @return a data frame.
#' @export

  model.frame.clust_analysis <- function(formula, ...) {

    stopifnot(all(class(formula) == 'clust_analysis'))

    rlang::eval_tidy(formula$data)

  }

#' Access the distance of a clust_analysis object.
#'
#' @description Retrieves the matrix with distances between the observations.
#' @param x a clust_analysis object.
#' @param ... extra arguments, currently none.
#' @return a numeric matrix.
#' @export

  dist.clust_analysis <- function(x, ...) {

    stopifnot(all(class(x) == 'clust_analysis'))

    x$dist_mtx

  }

#' Number of observations of a clust_analysis object.
#'
#' @description Retrieves the number of observations used to generate
#' a clust_analysis object.
#' @param object a clust_analysis object.
#' @param ... extra arguments, currently none.
#' @return a list with the numbers of observations and clustering features.
#' @export

  nobs.clust_analysis <- function(object, ...) {

    stopifnot(all(class(object) == 'clust_analysis'))

    clustTools:::get_data_dim(rlang::eval_tidy(object$data))

  }

#' Numbers of observations assigned to the clusters.
#'
#' @description Retrieves the numbers of observations assigned to the clusters.
#' @param x a clust_analysis object.
#' @param ... extra arguments, currently none.
#' @return a data frame with the numbers of observations in the clusters
#' @export

  ngroups.clust_analysis <- function(x, ...) {

    stopifnot(all(class(x) == 'clust_analysis'))

    dplyr::count(x$clust_assignment, clust_id)

  }

# Dimensionality reduction for the clust_analysis class -----

#' Dimensionality reduction analysis of the clust_analysis object data
#' or distance matrix.
#'
#' @description Performs principal component analysis (PCA), multi-dimensional
#' scaling (MDS) or uniform manifold approximation and projection (UMAP) of the
#' clust_analysis object data or distance matrix.
#' @details See \code{\link{reduce_data}} for the implementation details.
#' The distance method, relevant for MDS and UMAP. is taken over from the
#' clust_object. Hence, some distances may crash the analysis with UMAP, see:
#' \code{\link[umap]{umap.defaults}} for the compatible distances.
#' @param x a clust_analysis object.
#' @param kdim number of dimensions. If NULL, kdim is set to the number of
#' clusters.
#' @param red_fun reduction analysis function: 'pca' (PCA), 'mds' (MDS) or
#' 'umap' (UMAP).
#' @param with type of the input data for the reduction analysis:
#' the clustering data ('data') or the matrix of distances ('distance').
#' @param ... extra arguments passed to \code{\link{reduce_data}}.
#' @return a red_analysis object with the component/score table containing the
#' cluster assignment information ('clust_id' variable).
#' @export components.clust_analysis
#' @export

  components.clust_analysis <- function(x,
                                        kdim = NULL,
                                        red_fun = c('pca', 'mds', 'umap'),
                                        with = c('distance', 'data'), ...) {

    ## entry control

    stopifnot(all(class(x) == 'clust_analysis'))

    with <- match.arg(with[1], c('distance', 'data', 'umap'))

    red_fun <- match.arg(red_fun[1], c('pca', 'mds', 'umap'))

    if(is.null(kdim)) {

      kdim <- length(unique(x$clust_assignment$clust_id))

    }

    ## reduction

    red_obj <- reduce_data(extract.clust_analysis(x, type = with),
                           distance_method = x$dist_method,
                           kdim = kdim,
                           red_fun = red_fun, ...)

    red_obj$component_tbl <- dplyr::left_join(red_obj$component_tbl,
                                              x$clust_assignment,
                                              by = 'observation')

    red_obj

  }

# Clustering variance -------

#' Calculate clustering variance.
#'
#' @description Calculates the clustering sum of squares (total, within
#' clusters, total within clusters and between clusters) as well as the
#' fraction of 'explained' clustering variance. The later is the ratio of
#' between-cluster sum of squares to the total sum of squares statistic.
#' @param x a clust_analysis object.
#' @param ... extra arguments, currently none.
#' @export

  var.clust_analysis <- function(x, ...) {

    ## calculates within and total sum of squares for any clustering object

    stopifnot(all(class(x) == 'clust_analysis'))

    clustTools:::get_sum_sq(dist_mtx = x$dist_mtx,
                            assignment = x$clust_assignment)

  }

# Plotting ------

#' Plot selected features of a clust_analysis object.
#'
#' @description The plotting method for the clust_analysis class. Enables
#' plotting of the standard diagnostic plots used for the optimal cluster number
#' determination (dendrogram, WSS- and silhouette curve), results of the
#' reduction analysis, heat map of the distances between the observations as
#' well as the self-organizing map training process. It is also possible to plot
#' the first two variables of the clustering data frame, an option which is
#' attractive, if the clustering of reduction analysis was performed.
#' @param x a clust_analysis object.
#' @param type the type of plots:
#' 'diagnostic' returns a series of diagnostic plots,
#' for non-SOM clustering those include a dendrogram (hierarchical clustering),
#' WSS and silhouette curve (see: \code{\link{plot_nbclust}}) or the complete
#' output of \code{\link[kohonen]{plot.kohonen}};
#' 'components' plots the results of reduction analysis done with the clustering
#' data or the distance matrix (see: \code{\link{components.clust_analysis}});
#' 'heat_map' plots the distances between observations as a heat map, 'training'
#' plots the mean distance to the SOM winning unit as a function of the
#' iteration number; 'data' works only if reduction analysis results were used
#' for clustering and plots the first two comonents/dimensions.
#' @param cust_theme a ggplot theme.
#' @param jitter_width horizontal jittering of the point in the plots.
#' @param jitter_height vertical jittering of the point in the plots.
#' @param point_alpha scatter plot's point alpha.
#' @param ... extra arguments passed to \code{\link{components.clust_analysis}}.
#' @return a ggplot object or a list of ggplot objects, as specified by the
#' 'type' argument
#' @export plot.clust_analysis
#' @export

  plot.clust_analysis <- function(x,
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

    stopifnot(all(class(x) == 'clust_analysis'))

    type <- match.arg(type[1],
                      c('diagnostic',
                        'components',
                        'heat_map',
                        'training',
                        'data'))

    stopifnot(any(class(cust_theme) == 'theme'))

    ## plot meta

    clust_algorithm <- switch(x$clust_fun,
                              hclust = paste(',', x$hc_method, 'method'),
                              kmeans = '',
                              pam = '',
                              som = '')

    clust_method <- switch(x$clust_fun,
                           hclust = 'Hierarchical clustering',
                           kmeans = 'Kmeans clustering',
                           pam = 'PAM clustering',
                           som = 'SOM clustering',
                           dbscan = 'DBSCAN clustering',
                           prediction = 'prediction')

    plot_subtitle <- paste0(clust_method,
                            ', ',
                            x$dist_method,
                            ' distance',
                            clust_algorithm)

    plot_n <- nobs(x)

    plot_tag <- paste0('\nObservations: n = ',
                       plot_n$observations,
                       '\nVariables: n = ',
                       plot_n$variables)

    ## plots

    if(type == 'diagnostic') {

      ## diagnostic plots

      plot_list <- list()

      if(x$clust_fun %in% c('hclust', 'kmeans', 'pam')) {

        plot_title <- 'Optimal cluster number'

        k <- length(unique(x$clust_assignment$clust_id))

        clust_args <-  switch(x$clust_fun,
                              hclust = list(FUNcluster = factoextra::hcut,
                                            hc_method = x$hc_method),
                              kmeans = list(FUNcluster = stats::kmeans),
                              pam = list(FUNcluster = cluster::pam))

        plot_list[c('wss', 'silhouette')] <- purrr::map(c('wss', 'silhouette'),
                                                        function(met) rlang::call2('plot_nbclust',
                                                                                   data = x$dist_mtx,
                                                                                   method = met,
                                                                                   k = k,
                                                                                   !!!clust_args,
                                                                                   plot_title = plot_title,
                                                                                   plot_subtitle = plot_subtitle,
                                                                                   plot_tag = plot_tag,
                                                                                   cust_theme = cust_theme))

        plot_list <- purrr::map(plot_list, purrr::safely(eval))

        plot_list <- purrr::map(plot_list, ~.x$result)

        plot_list <- purrr::compact(plot_list)

        if(x$clust_fun == 'hclust') {

          plot_list$dendrogram <- plot_dendro(clust_str = x$clust_obj,
                                              k = k,
                                              labels = TRUE,
                                              cust_theme = cust_theme,
                                              plot_tag = plot_tag)

        }

      } else if(x$clust_fun == 'dbscan') {

        plot_list$knn_dist <- clustTools::plot_knn_distance(as.dist(x$dist_mtx),
                                                            k = x$minPts - 1,
                                                            eps = x$eps,
                                                            plot_title = 'kNN distance plot',
                                                            plot_subtitle = paste('DBSCAN algorithm, eps =',
                                                                                  x$eps),
                                                            plot_tag = plot_tag,
                                                            cust_theme = cust_theme)

      } else if(x$clust_fun == 'som') {

        plot_list <- clustTools::plot_som(x$clust_obj)

      } else {

        warning('No training plots available for cluster predictions', call. = FALSE)

        return(NULL)

      }

      return(plot_list)

    } else if(type == 'components') {

      ## MDS, PCA or UMAP

      red_results <- clustTools::components(x, ...)

      ## providing a 1d PCA, MDS or UMAP should throw an error

      if(!all(c('comp_1', 'comp_2') %in% names(red_results$component_tbl))) {

        stop('Atempt to plot a 1-dimensional reduction analysis result. Adjust kdim?', call. = FALSE)

      }

      sdevs <- var(red_results)

      if(red_results$red_fun == 'pca') {

        ax_labs <- purrr::map2(c('PC1', 'PC2'),
                               signif(sdevs$perc_var[1:2], 3),
                               ~paste0(.x, ', ', .y, '%'))

      } else {

        ax_labs <- purrr::map2(c('Dim 1', 'Dim 2'),
                               signif(sdevs$perc_var[1:2], 3),
                               ~paste0(.x, ', ', .y, '%'))

      }

      plot_point(data = red_results$component_tbl,
                 x_var = 'comp_1',
                 y_var = 'comp_2',
                 fill_var = 'clust_id',
                 plot_title = switch(red_results$red_fun ,
                                     pca = 'PCA',
                                     mds = 'MDS',
                                     umap = 'UMAP'),
                 plot_subtitle = plot_subtitle,
                 plot_tag = plot_tag,
                 x_lab = ax_labs[[1]],
                 y_lab = ax_labs[[2]],
                 cust_theme = cust_theme,
                 jitter_height = jitter_height,
                 jitter_width = jitter_width,
                 fill_lab = 'Cluster ID',
                 point_alpha = point_alpha)

    } else if(type == 'heat_map') {

      ## heat map of the distances between the observations

      plotting_tbl <- as.data.frame(x$dist_mtx)

      clust_vars <- colnames(plotting_tbl)

      plotting_tbl <- tibble::rownames_to_column(plotting_tbl,
                                                 'observation')

      plotting_tbl <- dplyr::left_join(plotting_tbl,
                                       x$clust_assignment[c('observation',
                                                            'clust_id')],
                                       by = 'observation')

      plotting_tbl <- tidyr::gather(plotting_tbl,
                                    key = 'observation2',
                                    value = 'distance',
                                    all_of(clust_vars))

      plotting_tbl <- dplyr::left_join(plotting_tbl,
                                       rlang::set_names(x$clust_assignment[c('observation',
                                                                             'clust_id')],
                                                        c('observation2',
                                                          'clust_id2')),
                                       by = 'observation2')

      heat_map <- ggplot2::ggplot(plotting_tbl,
                                  ggplot2::aes(x = reorder(observation, distance),
                                               y = reorder(observation2, distance),
                                               fill = distance)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = 'firebrick',
                                      mid = 'white',
                                      high = 'steelblue',
                                      midpoint = mean(range(plotting_tbl$distance))) +
        cust_theme +
        ggplot2::theme(axis.title = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(angle = 90,
                                                           hjust = 1,
                                                           vjust = 0.5)) +
        ggplot2::labs(title = 'Distances between observations',
                      subtitle = plot_subtitle,
                      tag = plot_tag,
                      fill = 'Distance') +
        ggplot2::facet_grid(clust_id2 ~ clust_id,
                            scales = 'free',
                            space = 'free')

      return(heat_map)

    } else if(type == 'training') {

      if(x$clust_fun != 'som') {

        warning('The training plots available only for the SOM cluster analyses', call. = FALSE)

        return(NULL)

      }

      som_training <- clustTools::plot_train_som(kohonen_object = x$clust_obj,
                                                 plot_title = 'SOM clustering: training',
                                                 plot_subtitle = plot_subtitle,
                                                 cust_theme = cust_theme)

      return(som_training)

    } else {

      plotting_tbl <- model.frame(x)

      if(!all(c('comp_1', 'comp_2') %in% names(plotting_tbl))) {

        warning('The data plots are available only, if a red_analysis object was subjected to clustering.',
                call. = FALSE)

        return(NULL)

      }

      ## component variances

      var_list <- dplyr::select(plotting_tbl, dplyr::starts_with('comp'))

      var_list <- purrr::map_dbl(var_list, var)

      perc_var <- var_list/sum(var_list) * 100

      plotting_tbl <- tibble::rownames_to_column(plotting_tbl, 'observation')

      plotting_tbl <- dplyr::left_join(plotting_tbl,
                                       x$clust_assignment,
                                       by = 'observation')

      point_plot <- clustTools::plot_point(data = plotting_tbl,
                                           x_var = 'comp_1',
                                           y_var = 'comp_2',
                                           fill_var = 'clust_id',
                                           plot_title = 'Input reduction analysis',
                                           plot_subtitle = plot_subtitle,
                                           plot_tag = plot_tag,
                                           x_lab = paste0('Dim 1, ',
                                                          signif(perc_var[1], 2),
                                                          '%'),
                                           y_lab = paste0('Dim 2, ',
                                                          signif(perc_var[2], 2),
                                                          '%'),
                                           cust_theme = cust_theme,
                                           point_alpha = point_alpha,
                                           jitter_width = jitter_width,
                                           jitter_height = jitter_width)

      return(point_plot)

    }

  }

# Semi-supervised clustering ------

#' Semi-supervised clustering.
#'
#' @description Projects the cluster assignment onto new data using simple
#' observation matching or a k-nearest neighbor (kNN) label propagation
#' algorithm.
#' @details For the implementation details, see: \code{\link{propagate}}.
#' The default distance metric is extracted from the clust_analysis object.
#' @param object a clust_analysis object.
#' @param newdata a numeric data frame, matrix or a red_analysis object. If NULL
#' (default), the bare cluster assignment table is returned.
#' @param type type of the projection: simple observation matching
#' ('class', default) or kNN label propagation ('propagation').
#' @param ... extra arguments passed to \code{\link{propagate}}.
#' @return a clust_analysis object.
#' @export predict.clust_analysis
#' @export

  predict.clust_analysis <- function(object,
                                     newdata = NULL,
                                     type = c('class', 'propagation'), ...) {

    ## entry control

    stopifnot(all(class(object) == 'clust_analysis'))

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

    ## projections

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
                                          clust_id = train_assignment[rownames(newdata), 'clust_id'])

      } else {

        warning('Unnamed observations in new data')

        test_assignment <- tibble::rownames_to_column(train_assignment,
                                                      'observation')

        test_assignment <- tibble::as_tibble(test_assignment)

      }

      ## output

      model_frame <- rlang::enexpr(newdata)

      clustTools::clust_analysis(list(data = rlang::quo(!!model_frame),
                                      dist_mtx = calculate_dist(newdata, method = object$dist_method),
                                      dist_method = object$dist_method,
                                      clust_fun = 'prediction',
                                      clust_obj = NULL,
                                      clust_assignment = test_assignment,
                                      dots = rlang::list2()))

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
#' @param x a clust_analysis object.
#' @param nfolds number of CV folds.
#' @param kNN number of the nearest neighbors.
#' @param simple_vote logical, should classical unweighted k-NN classification
#' be applied? If FALSE, distance-weighted k-NN is used with the provided kernel
#' function.
#' @param resolve_ties logical, should the ties be resolved at random? Applies
#' only to the simple unweighted voting algorithm.
#' @param kernel_fun kernel function transforming the distance into weight.
#' @param seed initial setting of the random number generator.
#' @param .parallel logical, should the CV be run in parallel? Experimental.
#' @return a list containing the global clust_analysis object, projection
#' (prediction) results and prediction summary for each fold and a prediction
#' summary for the whole CV.
#' @export cv.clust_analysis
#' @export

  cv.clust_analysis <- function(x,
                                nfolds = 5,
                                kNN = 5,
                                simple_vote = TRUE,
                                resolve_ties = FALSE,
                                kernel_fun = function(x) 1/x,
                                seed = 1234,
                                .parallel = FALSE) {

    ## entry control

    stopifnot(all(class(x) == 'clust_analysis'))

    nfolds <- as.integer(nfolds)

    kNN <- as.integer(kNN)

    stopifnot(is.logical(simple_vote))
    stopifnot(is.logical(resolve_ties))
    stopifnot(is.logical(.parallel))

    stopifnot(is.function(kernel_fun))

    ## common parameters

    cmm_args <- rlang::list2(data = rlang::eval_tidy(x$data),
                             nfolds = nfolds,
                             kNN = kNN,
                             simple_vote = simple_vote,
                             resolve_ties = resolve_ties,
                             kernel_fun = kernel_fun,
                             distance_method = x$dist_method,
                             seed = seed,
                             .parallel = .parallel,
                             k = nrow(ngroups(x)),
                             hc_method = x$hc_method,
                             clust_fun = x$clust_fun,
                             xdim = x$grid$xdim,
                             ydim = x$grid$ydim,
                             topo = x$grid$topo,
                             neighbourhood.fct = as.character(x$grid$neighbourhood.fct),
                             toroidal = x$grid$toroidal,
                             eps = x$eps,
                             minPts = x$minPts)

    cmm_args <- c(cmm_args, x$dots)

    ## function calls

    if(x$clust_fun %in% c('hclust', 'som')) {

      cmm_args$clust_fun <- NULL

    } else if(x$clust_fun == 'dbscan') {

      cmm_args$k <- NULL

      cmm_args$clust_fun <- NULL

    }

    cv_call <- rlang::call2('cv_cluster',
                            !!!purrr::compact(cmm_args),
                            clustering_fun = switch(x$clust_fun,
                                                    hclust = hcluster,
                                                    kmeans = kcluster,
                                                    pam = kcluster,
                                                    dbscan = dbscan_cluster,
                                                    som = som_cluster))

    eval(cv_call)

  }

# Variable importance ------

#' Determine clustering feature importance.
#'
#' @description Determines importance of specific clustering variables by
#' comparing the fraction of 'explained' clustering variance of the input
#' clustering object and the object generated with the variable
#' re-shuffled randomly.
#' @param x a clust_analysis object.
#' @param seed initial setting of the random number generator.
#' @param .parallel logical, should the CV be run in parallel? Experimental.
#' @return a data frame of class 'importance'
#' @export impact.clust_analysis
#' @export

  impact.clust_analysis <- function(x,
                                    seed = 1234,
                                    .parallel = FALSE) {

    ## entry control

    stopifnot(all(class(x) == 'clust_analysis'))

    ## common parameters

    cmm_args <- rlang::list2(data = rlang::eval_tidy(x$data),
                             distance_method = x$dist_method,
                             seed = seed,
                             .parallel = .parallel,
                             k = nrow(ngroups(x)),
                             hc_method = x$hc_method,
                             clust_fun = x$clust_fun,
                             xdim = x$grid$xdim,
                             ydim = x$grid$ydim,
                             topo = x$grid$topo,
                             eps = x$eps,
                             minPts = x$minPts,
                             neighbourhood.fct = as.character(x$grid$neighbourhood.fct),
                             toroidal = x$grid$toroidal)

    cmm_args <- c(cmm_args, x$dots)

    ## calls

    if(x$clust_fun %in% c('hclust', 'som')) {

      cmm_args$clust_fun <- NULL

    } else if(x$clust_fun == 'dbscan') {

      cmm_args$k <- NULL

      cmm_args$clust_fun <- NULL

    }

    imp_call <- rlang::call2('importance_cluster',
                             !!!purrr::compact(cmm_args),
                             clustering_fun = switch(x$clust_fun,
                                                     hclust = hcluster,
                                                     kmeans = kcluster,
                                                     pam = kcluster,
                                                     dbscan = dbscan_cluster,
                                                     som = som_cluster))

    eval(imp_call)

  }

# END ------
