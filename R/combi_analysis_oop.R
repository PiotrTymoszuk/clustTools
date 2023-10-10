# Specific S3 OOP methods for the combi_analysis class.

# Plotting --------

#' Plot selected features of a combi_analysis object.
#'
#' @description
#' The plotting method for the `combi_analysis` class. Enables
#' plotting of the standard diagnostic plots used for the optimal cluster number
#' determination for the node clustering (dendrogram, WSS- and
#' silhouette curve), results of the reduction analysis,
#' heat map of the distances between the observations and SOM nodes as
#' well as the self-organizing map training process. It is also possible to plot
#' the first two variables of the clustering data frame, an option which is
#' attractive, if the clustering of reduction analysis was performed.
#'
#' @references
#' Wehrens R, Kruisselbrink J. Flexible self-organizing maps in kohonen 3.0.
#' J Stat Softw (2018) 87:1–18. doi:10.18637/jss.v087.i07
#' @references
#' Kassambara A, Mundt F. factoextra: Extract and Visualize the Results
#' of Multivariate Data Analyses. (2020) Available
#' at: https://cran.r-project.org/web/packages/factoextra/index.html
#' @references
#' Galili T. dendextend: an R package for visualizing, adjusting and
#' comparing trees of hierarchical clustering.
#' Bioinformatics (2015) 31:3718–20. doi:10.1093/bioinformatics/btv428
#' @references
#' McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation
#' and Projection for Dimension Reduction. (2018) Available
#' at: https://arxiv.org/abs/1802.03426v3
#' @references
#' Belyadi H, Haghighat A, Nguyen H, Guerin A-J. IOP Conference Series:
#' Earth and Environmental Science Determination of Optimal Epsilon (Eps)
#' Value on DBSCAN Algorithm to Clustering Data on Peatland Hotspots in
#' Sumatra Related content EPS conference comes to London-EPS rewards
#' quasiparticle research-EP. IOP Conf Ser Earth Environ Sci (2016) 31:
#' doi:10.1088/1755-1315/31/1/012012
#' @references
#' Rousseeuw PJ. Silhouettes: A graphical aid to the interpretation and
#' validation of cluster analysis. J Comput Appl Math (1987) 20:53–65.
#' doi:10.1016/0377-0427(87)90125-7
#'
#' @param x a `combi_analysis` object.
#' @param type the type of plots:
#'
#' * `diagnostic` returns a series of diagnostic plots for the SOM construction
#' and clustering of the nodes
#'
#' * `components` plots the results of reduction analysis done with the
#' clustering data or the distance matrix
#' (see: \code{\link{components.clust_analysis}})
#'
#' * `heat_map` plots the distances between observations and nodes as heat maps
#'
#' * `training` plots the mean distance to the SOM winning unit as a function
#' of the iteration number
#'
#' * `data` works only if reduction analysis results were used for clustering
#' and plots the first two components/dimensions.
#'
#' @inheritParams plot.clust_analysis
#'
#' @return a ggplot object or a list of ggplot objects, as specified by the
#' 'type' argument
#'
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

    ## entry control --------

    stopifnot(is_combi_analysis(x))

    type <- match.arg(type[1],
                      c('diagnostic',
                        'components',
                        'heat_map',
                        'training',
                        'data'))

    stopifnot(inherits(cust_theme, 'theme'))
    stopifnot(is.numeric(jitter_width))
    stopifnot(is.numeric(jitter_height))
    stopifnot(is.numeric(point_alpha))

    observation <- NULL

    ## plotting -------

    plot_list <- map(x$clust_analyses,
                     plot.clust_analysis,
                     type = type,
                     cust_theme = cust_theme,
                     jitter_height = jitter_height,
                     jitter_width = jitter_width, ...)

    plot_list <- compact(plot_list)

    if(type == 'training') {

      plot_list$observation <-
        plot_list$observation +
        ggplot2::labs(tag = paste('Iterations: n =',
                                  nrow(x$clust_analyses$observation$clust_obj$changes)))

    }

    ## summary component plots for the final clustering results

    if(type == 'components') {

      obs_red <- components(x$clust_analyses$observation, ...)

      score_tbl <- select(obs_red$component_tbl,
                          observation,
                          dplyr::starts_with('comp'))

      score_tbl <- left_join(score_tbl,
                             x$clust_assignment,
                             by = 'observation')

      ## handling U-matrix requests

      if(all(is.na(score_tbl$clust_id))) return(plot_list)

      sdevs <- var(obs_red)

      if(obs_red$red_fun == 'pca') {

        ax_labs <- map2(c('PC1', 'PC2'),
                        signif(sdevs$perc_var[1:2], 3),
                        ~paste0(.x, ', ', .y, '%'))

      } else {

        ax_labs <- map2(c('Dim 1', 'Dim 2'),
                        signif(sdevs$perc_var[1:2], 3),
                        ~paste0(.x, ', ', .y, '%'))

      }

      plot_list$final <-
        plot_point(data = score_tbl,
                   x_var = 'comp_1',
                   y_var = 'comp_2',
                   fill_var = 'clust_id',
                   plot_title = switch(obs_red$red_fun,
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

    return(plot_list)

  }

#' @rdname plot.combi_analysis
#' @export

  plot.umatrix_analysis <- function(x,
                                    type = c('diagnostic',
                                             'components',
                                             'heat_map',
                                             'training',
                                             'data'),
                                    cust_theme = ggplot2::theme_classic(),
                                    jitter_width = 0,
                                    jitter_height = 0,
                                    point_alpha = 1, ...) {

    ## entry control --------

    stopifnot(is_combi_analysis(x))

    type <- match.arg(type[1],
                      c('diagnostic',
                        'components',
                        'heat_map',
                        'training',
                        'data'))

    stopifnot(inherits(cust_theme, 'theme'))
    stopifnot(is.numeric(jitter_width))
    stopifnot(is.numeric(jitter_height))
    stopifnot(is.numeric(point_alpha))

    observation <- NULL

    extra_args <- list2(...)

    ## most plotting cases are managed by the superclass -----

    if(type %in% c('diagnostic',
                   'heat_map',
                   'training',
                   'data')) {

      return(NextMethod())

    }

    if(!is.null(extra_args$with)) {

      if(extra_args$with != 'data') return(NextMethod())

    }

    ## in case component plots are requested for 'data' -------
    ## plots will be created for each layer

    red_objects <- components(object = x, ...)

    score_tbl <- map(red_objects, extract, 'scores')

    sdevs <- map(red_objects, var)

    if(red_objects[[1]]$red_fun == 'pca') {

      ax_labs <-
        map(sdevs,
            function(layer) map2(c('PC1', 'PC2'),
                                 signif(layer$perc_var[1:2], 3),
                                 ~paste0(.x, ', ', .y, '%')))

    } else {

      ax_labs <-
        map(sdevs,
            function(layer) map2(c('Dim1', 'Dim2'),
                                 signif(layer$perc_var[1:2], 3),
                                 ~paste0(.x, ', ', .y, '%')))

    }

    n_numbers <- map(red_objects, nobs)

    plot_tags <- map(n_numbers,
                     ~paste0('Observations: n = ', .x[['observations']],
                             '\nVariables: n = ', .x[['variables']]))

    pmap(list(data = score_tbl,
              x_lab = map(ax_labs, ~.x[[1]]),
              y_lab = map(ax_labs, ~.x[[2]]),
              plot_tag = plot_tags),
         plot_point,
         x_var = 'comp_1',
         y_var = 'comp_2',
         fill_var = 'clust_id',
         plot_title = switch(red_objects[[1]]$red_fun,
                             pca = 'PCA',
                             mds = 'MDS',
                             umap = 'UMAP'),
         cust_theme = cust_theme,
         jitter_height = jitter_height,
         jitter_width = jitter_width,
         fill_lab = 'Cluster ID',
         point_alpha = point_alpha)

  }

# END ------
