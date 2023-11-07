# S3 methods for the `knb` class

# Summary -------

#' Summary neighborhood preservation statistic for clusters.
#'
#' @description
#' Computes mean, SD, median, interquartile range, 95% range and range of the
#' network preservation statistic for the global clustering structure
#' and particular clusters.
#' Low values of the median or median neighborhood preservation statistic may
#' indicate that the given cluster s poorly separated from other clusters.
#' For objects without cluster assignment (e.g.
#' neighborhood analyses for reduction analysis methods), only global
#' neighborhood statistics are returned.
#' For such objects, low values of the neighborhood preservation statistic
#' suggest that the reduction analysis poorly projects the data point
#' neighborhood into the reduced layout.
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
#' @param object an object of the \code{\link{knb}} class.
#' @param ... extra arguments, currently none.
#'
#' @return  a data frame with numeric statistics for the whole clustering
#' structure (`clust_id` = 'global') and particular clusters.
#'
#' @export summary.knb
#' @export

  summary.knb <- function(object, ...) {

    stopifnot(is_knb(object))

    stat_name <-
      names(object)[!names(object) %in% c('observation',
                                          'clust_id',
                                          'kNN_data',
                                          'kNN_cluster')]

    stat_funs <-
      list(n = function(x) length(x),
           mean = function(x) mean(x, na.rm = TRUE),
           sd = function(x) stats::sd(x, na.rm = TRUE),
           median = function(x) stats::median(x, na.rm = TRUE),
           q025 = function(x) stats::quantile(x, 0.025, na.rm = TRUE),
           q25 = function(x) stats::quantile(x, 0.25, na.rm = TRUE),
           q75 = function(x) stats::quantile(x, 0.75, na.rm = TRUE),
           q975 = function(x) stats::quantile(x, 0.975, na.rm = TRUE),
           min = function(x) min(x, na.rm = TRUE),
           max = function(x) max(x, na.rm = TRUE))

    clust_data <- split(object[[stat_name]], object$clust_id)

    ## global and cluster stats -------

    global_stats <- map_dbl(stat_funs, ~.x(object[[stat_name]]))

    clust_stats <- map(clust_data,
                       function(clust) map(stat_funs, ~.x(clust)))

    stats <- c(list(global = global_stats),
               clust_stats)

    stats <- map(stats, reduce, cbind)

    stats <- map(stats, as.data.frame)

    stats <- map(stats, set_names, names(stat_funs))

    clust_id <- NULL

    stats <-
      map2_dfr(stats, names(stats),
               ~mutate(.x, clust_id = .y))

    stats <-
      mutate(stats,
             clust_id = factor(clust_id, c('global', names(clust_data))))

    dplyr::relocate(as_tibble(stats), clust_id)

  }

# Plotting --------

#' Plot distribution of neighborhood preservation statistic.
#'
#' @description
#' Generates a bar plot of neighborhood preservation statistic values for
#' observations and clusters, similar to a classical silhouette plot (see:
#' \code{\link{plot.sil_extra}}).
#'
#' @param x an object of the \code{\link{knb}} class.
#' @param bar_color color of the bar line.
#' @param bar_fill color of the bars, relevant only for objects without
#' cluster assignment.
#' @param show_stats logical, should the number of observations in the cluster,
#' percentage of negative silhouette widths and average silhouette statistic
#' be shown in the plot? Defaults to TRUE.
#' @param signif_digits significant digits used for rounding of the statistics
#' presented in the plot.
#' @param cust_theme custom ggplot theme.
#' @param ... extra arguments passed to \code{\link[ggplot2]{geom_bar}}.
#'
#' @return a `ggplot` class graphic.
#'
#' @export plot.knb
#' @export

  plot.knb <- function(x,
                       show_stats = TRUE,
                       signif_digits = 2,
                       cust_theme = ggplot2::theme_classic(),
                       bar_color = 'black',
                       bar_fill = 'steelblue', ...) {

    ## entry check ------

    stopifnot(is_knb(x))

    stat_name <-
      names(x)[!names(x) %in% c('observation',
                                'clust_id',
                                'kNN_data',
                                'kNN_cluster')]

    stopifnot(is.logical(show_stats))
    stopifnot(is.numeric(signif_digits))

    signif_digits <- as.integer(signif_digits)

    if(!inherits(cust_theme, 'theme')) {

      stop("'cust_theme' hast to be a vaild ggplot2 'theme' object.",
           call. = FALSE)

    }

    ## labeller ----------

    if(!show_stats) {

      stat_labs <- 'label_value'

    } else {

      plot_lab <- NULL
      n <- NULL

      stats <-
        mutate(summary(x),
               plot_lab = paste0('total: n = ', n,
                                 '\navg = ', signif(mean, signif_digits)),
               plot_lab = paste(clust_id, plot_lab, sep = '\n'))

      stat_labs <- set_names(stats$plot_lab, stats$clust_id)

    }

    ## plotting ------

    clust_id <- NULL
    observation <- NULL

    if(!all(is.na(x$clust_id))) {

      plot_subtitle <- NULL

      if(!all(is.na(x$kNN_data))) {

        plot_subtitle <-
          paste0(x$kNN_data[[1]], '-nearest data point neighbors')

      }

      if(!all(is.na(x$kNN_cluster))) {

        plot_subtitle <- paste0(plot_subtitle, ', ', x$kNN_cluster[[1]],
                                '-nearest cluster neighbors')

      }

      base_plot <- ggplot(x,
                          aes(x = .data[[stat_name]],
                              y = reorder(observation, .data[[stat_name]]),
                              fill = clust_id)) +
        ggplot2::geom_bar(stat = 'identity',
                          color = bar_color, ...) +
        ggplot2::facet_grid(clust_id ~ .,
                            scales = 'free',
                            space = 'free',
                            labeller = ggplot2::as_labeller(stat_labs)) +
        ggplot2::labs(title = 'Neighborhood preservation',
                      subtitle = plot_subtitle,
                      y = 'observation',
                      x = 'Fraction of preserved nearest neighbors',
                      fill = 'Cluster')

    } else {

      if(show_stats) {

        plot_subtitle <-
          stringi::stri_replace_all(stat_labs[['global']],
                                    fixed = '\n',
                                    replacement = ', ')

      }

      base_plot <-
        ggplot(x,
               aes(x = .data[[stat_name]],
                   y = reorder(observation, .data[[stat_name]]))) +
        ggplot2::geom_bar(stat = 'identity',
                          color = bar_color,
                          fill = bar_fill, ...) +
        ggplot2::labs(title = 'Neighborhood preservation',
                      subtitle = plot_subtitle,
                      y = 'observation',
                      x = 'Fraction of preserved nearest neighbors')

    }

    base_plot +
      ggplot2::geom_vline(xintercept = 0,
                          linetype = 'dashed') +
      cust_theme +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())

  }


# END ------
