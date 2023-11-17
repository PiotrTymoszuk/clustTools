# Specific S3 OOP for the clust_analysis class.

# Plotting ------

#' Plot selected features of a clust_analysis object.
#'
#' @description
#' The plotting method for the `clust_analysis` class. Enables
#' plotting of the standard diagnostic plots used for the optimal cluster number
#' determination (dendrogram, WSS- and silhouette curve), results of the
#' reduction analysis, heat map of the distances between the observations as
#' well as the self-organizing map training process. It is also possible to plot
#' the first two variables of the clustering data frame, an option which is
#' attractive, if the clustering of reduction analysis was performed.
#'
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
#' @param x a `clust_analysis` object.
#' @param type the type of plots:
#'
#' * 'diagnostic' returns a series of diagnostic plots,
#' for non-SOM clustering those include a dendrogram (hierarchical clustering),
#' WSS and silhouette curve (see: \code{\link{plot_nbclust}}) or the complete
#' output of \code{\link[kohonen]{plot.kohonen}}
#'
#' * 'components' plots the results of reduction analysis done with the
#' clustering data, distance matrix or, for SOM, with U matrix
#' (see: \code{\link{components.clust_analysis}})
#'
#' * 'heat_map' plots the distances between observations as a heat map
#'
#' * 'training' plots the mean distance to the SOM winning unit as a function
#' of the iteration number
#'
#' * 'data' works only if reduction analysis results were used
#' for clustering and plots the first two components/dimensions.
#'
#' @param cust_theme a ggplot theme.
#' @param jitter_width horizontal jittering of the points in the plots.
#' @param jitter_height vertical jittering of the points in the plots.
#' @param point_alpha scatter plot's point alpha.
#' @param ... extra arguments passed to \code{\link{components.clust_analysis}}.
#'
#' @return a `ggplot` object or a list of `ggplot` objects, as specified by the
#' 'type' argument and character of the object.
#'
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

    ## entry control -----

    stopifnot(is_clust_analysis(x))

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
    distance <- NULL
    observation2 <- NULL

    ## plot meta -------

    clust_algorithm <- switch(x$clust_fun,
                              hclust = paste(',', x$hc_method, 'method'),
                              kmeans = '',
                              htk = '',
                              pam = '',
                              som = ', single-layer',
                              supersom = ', multi-layer',
                              dbscan = '')

    clust_method <- switch(x$clust_fun,
                           hclust = 'Hierarchical clustering',
                           kmeans = 'Kmeans clustering',
                           htk = 'HT-KMEANS clustering',
                           pam = 'PAM clustering',
                           som = 'SOM clustering',
                           supersom = 'SOM clustering',
                           dbscan = 'DBSCAN clustering',
                           prediction = 'prediction',
                           supersom_prediction = 'prediction')

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

    ## plots -----------

    if(type == 'data' & x$clust_fun %in% c('supersom', 'supersom_prediction')) {

      warning("The 'data' option is not available for multi-layer SOM.",
              call. = FALSE)

      return(NULL)

    }

    if(type == 'diagnostic') {

      ## diagnostic plots

      plot_list <- list()

      if(x$clust_fun %in% c('hclust', 'kmeans', 'pam')) {

        k <- length(unique(x$clust_assignment$clust_id))

        clust_args <-  switch(x$clust_fun,
                              hclust = list(FUNcluster = factoextra::hcut,
                                            hc_method = x$hc_method),
                              kmeans = list(FUNcluster = stats::kmeans),
                              pam = list(FUNcluster = cluster::pam))

        plot_list[c('wss', 'silhouette')] <-
          map(c('wss', 'silhouette'),
              function(met) call2('plot_nbclust',
                                  data = x$dist_mtx,
                                  method = met,
                                  diss = dist(x, 'distance'),
                                  k = k,
                                  !!!clust_args,
                                  plot_title = 'Optimal cluster number',
                                  plot_subtitle = plot_subtitle,
                                  plot_tag = plot_tag,
                                  cust_theme = cust_theme))

        plot_list <- map(plot_list, purrr::safely(eval))

        plot_list <- map(plot_list, ~.x$result)

        plot_list <- compact(plot_list)

        if(x$clust_fun == 'hclust') {

          plot_list$dendrogram <- plot_dendro(clust_str = x$clust_obj,
                                              k = k,
                                              labels = TRUE,
                                              cust_theme = cust_theme,
                                              plot_tag = plot_tag)

        }

      } else if(x$clust_fun == 'dbscan') {

        plot_list$knn_dist <-
          plot_knn_distance(as.dist(x$dist_mtx),
                            k = x$minPts - 1,
                            eps = x$eps,
                            plot_title = 'kNN distance plot',
                            plot_subtitle = paste('DBSCAN algorithm, eps =',
                                                  x$eps),
                            plot_tag = plot_tag,
                            cust_theme = cust_theme)

      } else if(x$clust_fun %in% c('som', 'supersom')) {

        plot_list <- plot_som(x$clust_obj)

      } else if(x$clust_fun == 'htk') {

        plot_list <- plot_htk(x = x,
                              plot_title = 'Optimal cluster number',
                              plot_subtitle = plot_subtitle,
                              plot_tag = plot_tag,
                              cust_theme = cust_theme)

      } else {

        warning('No training plots available for cluster predictions',
                call. = FALSE)

        return(NULL)

      }

      return(plot_list)

    }

    if(type == 'components') {

      ## MDS, PCA or UMAP

      red_results <- components(x, ...)

      if(is.null(red_results)) return(NULL)

      if(is_red_analysis(red_results)) {

        point_plot <- plot.clust_red(red_results,
                                     type = 'score',
                                     label_clust = TRUE,
                                     cust_theme = cust_theme,
                                     plot_subtitle = plot_subtitle,
                                     jitter_height = jitter_height,
                                     jitter_width = jitter_width,
                                     point_alpha = point_alpha)

        return(point_plot)

      } else {

        point_plots <-
          pmap(list(x = red_results),
               plot.clust_red,
               type = 'score',
               label_clust = TRUE,
               cust_theme = cust_theme,
               plot_subtitle = plot_subtitle,
               jitter_height = jitter_height,
               jitter_width = jitter_width,
               point_alpha = point_alpha)

        return(point_plots)

      }

    }

    if(type == 'heat_map') {

      ## heat map of the distances between the observations

      plotting_tbl <- as.data.frame(x$dist_mtx)

      clust_vars <- colnames(plotting_tbl)

      plotting_tbl <- rownames_to_column(plotting_tbl,
                                         'observation')

      plotting_tbl <- left_join(plotting_tbl,
                                x$clust_assignment[c('observation',
                                                     'clust_id')],
                                by = 'observation')

      plotting_tbl <-
        tidyr::pivot_longer(plotting_tbl,
                            names_to = 'observation2',
                            values_to = 'distance',
                            cols = all_of(clust_vars))

      plotting_tbl <-
        left_join(plotting_tbl,
                  set_names(x$clust_assignment[c('observation',
                                                 'clust_id')],
                            c('observation2',
                              'clust_id2')),
                  by = 'observation2')

      heat_map <-
        ggplot(plotting_tbl,
               aes(x = stats::reorder(observation, distance),
                   y = stats::reorder(observation2, distance),
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

    }

    if(type == 'training') {

      if(!x$clust_fun %in% c('som', 'supersom')) {

        warning('The training plots available only for the SOM cluster analyses', call. = FALSE)

        return(NULL)

      }

      som_training <-
        plot_train_som(kohonen_object = x$clust_obj,
                       plot_title = 'SOM clustering: training',
                       plot_subtitle = plot_subtitle,
                       cust_theme = cust_theme)

      return(som_training)

    }

    if(type == 'data') {

      plotting_tbl <- model.frame(x)

      if(!all(c('comp_1', 'comp_2') %in% names(plotting_tbl))) {

        warning('The data plots are available only, if a red_analysis object was subjected to clustering.',
                call. = FALSE)

        return(NULL)

      }

      ## component variances

      var_list <- select(plotting_tbl, dplyr::starts_with('comp'))

      var_list <- map_dbl(var_list, var)

      perc_var <- var_list/sum(var_list) * 100

      plotting_tbl <- rownames_to_column(plotting_tbl, 'observation')

      plotting_tbl <- left_join(plotting_tbl,
                                x$clust_assignment,
                                by = 'observation')

      point_plot <-
        plot_point(data = plotting_tbl,
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

# Summary quality stats -------

#' Quality control of clustering solutions.
#'
#' @description
#' Computes basic global statistics of quality of a clustering analysis object.
#'
#' @details
#' The statistics retrieved by the `summary()` method are:
#'
#' * _silhouette width_ (`sil_width`)
#'
#' * _fraction of potentially misclassified observations_ with negative
#' silhouette widths (`frac_misclassified`)
#'
#' * _fraction of explained clustering variance_ expressed as the ratio of total
#' between sum of squares to total sum of squares (`frac_var`)
#'
#' * _fraction of preserved nearest neighbors_ (`frac_np`)
#'
#' The statistics are computed with \code{\link{silhouette}}, \code{\link{var}},
#' and \code{\link{np}} methods for the entire clustering structure and not
#' for particular clusters.
#'
#' @references
#' Rousseeuw PJ. Silhouettes: A graphical aid to the interpretation and
#' validation of cluster analysis. J Comput Appl Math (1987) 20:53–65.
#' doi:10.1016/0377-0427(87)90125-7
#'
#' @references
#' Venna J, Kaski S. Neighborhood preservation in nonlinear projection methods:
#' An experimental study. Lect Notes Comput Sci (including Subser Lect Notes
#' Artif Intell Lect Notes Bioinformatics) (2001) 2130:485–491.
#' doi:10.1007/3-540-44668-0_68
#'
#' @return a data frame with columns characterized in Details.
#'
#' @param object a `clust_analysis` or `combi_analysis` object.
#' @param ... extra arguments passed to \code{\link{np}}.
#'
#' @export summary.clust_analysis
#' @export

  summary.clust_analysis <- function(object, ...) {

    stopifnot(is_clust_analysis(object) | is_combi_analysis(object))

    sil_res <- try(silhouette(object), silent = TRUE)

    variance <- try(var(object), silent = TRUE)

    neighbor_preservation <- try(np(object, type = 'final', ...), silent = TRUE)

    stats <-
      tibble(sil_width = tryCatch(mean(sil_res$sil_width, na.rm = TRUE),
                                  error = function(e) NA),
             frac_misclassified = tryCatch(sum(sil_res$sil_width < 0)/nrow(sil_res),
                                           error = function(e) NA),
             frac_var = tryCatch(variance$frac_var,
                                 error = function(e) NA),
             frac_np = tryCatch(mean(neighbor_preservation$frac_np, na.rm = TRUE),
                                error = function(e) NA))

    if(any(map_lgl(stats, is.na))) {

      warning('At least one statistic could not be calculated.',
              call. = FALSE)

    }

    stats

  }

# END ------
