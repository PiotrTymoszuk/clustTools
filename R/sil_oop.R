# Class-specific S3 OOP methods for 'sil_extra' objects

# Summary -------

#' Summary silhouette width statistic for clusters.
#'
#' @description
#' Computes mean, SD, median, interquartile range, 95% range, range as well as
#' the number and percentage of observations with negative silhouette width.
#' Such observations are likely in an improper cluster.
#'
#' @references
#' Rousseeuw PJ. Silhouettes: A graphical aid to the interpretation and
#' validation of cluster analysis. J Comput Appl Math (1987) 20:53–65.
#' doi:10.1016/0377-0427(87)90125-7
#'
#' @param object an object of the \code{\link{sil_extra}} class.
#' @param ... extra arguments, currently none.
#'
#' @return  a data frame with numeric statistics for the whole clustering
#' structure (`clust_id` = 'global') and particular clusters (mean, SD, median,
#' interquartile range, 95% percentile range, range, number and fraction of
#' potentially misclassified observations with negative silhouette widths).
#'
#' @export summary.sil_extra
#' @export

  summary.sil_extra <- function(object, ...) {

    stopifnot(is_sil_extra(object))

    stat_funs <-
      list(mean = function(x) mean(x, na.rm = TRUE),
           sd = function(x) stats::sd(x, na.rm = TRUE),
           median = function(x) stats::median(x, na.rm = TRUE),
           q025 = function(x) stats::quantile(x, 0.025, na.rm = TRUE),
           q25 = function(x) stats::quantile(x, 0.25, na.rm = TRUE),
           q75 = function(x) stats::quantile(x, 0.75, na.rm = TRUE),
           q975 = function(x) stats::quantile(x, 0.975, na.rm = TRUE),
           min = function(x) min(x, na.rm = TRUE),
           max = function(x) max(x, na.rm = TRUE))

    clust_data <- split(object$sil_width, object$clust_id)

    ## global and cluster stats -------

    global_stats <- map_dbl(stat_funs, ~.x(object$sil_width))

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

    ## n numbers and observations with negative metric values -------

    global_neg <- tibble(n = length(object$sil_width),
                         n_negative = sum(object$sil_width < 0))

    cluster_neg <- map(clust_data,
                       ~tibble(n = length(.x),
                               n_negative = sum(.x < 0)))

    neg <- c(list(global = global_neg),
             cluster_neg)

    neg <- map2_dfr(neg, names(neg),
                    ~mutate(.x,
                            clust_id = .y,
                            frac_misclassified = n_negative/n))

    stats <- left_join(neg, stats, by = 'clust_id')

    stats <-
      mutate(stats,
             clust_id = factor(clust_id, c('global', names(clust_data))))

    dplyr::relocate(as_tibble(stats), clust_id)

  }

# Plotting ------

#' Plots of silhouette statistics.
#'
#' @description
#' Generates a classical bar plot of silhouette width distribution in clusters.
#'
#' @references
#' Rousseeuw PJ. Silhouettes: A graphical aid to the interpretation and
#' validation of cluster analysis. J Comput Appl Math (1987) 20:53–65.
#' doi:10.1016/0377-0427(87)90125-7
#'
#' @param x an object of the \code{\link{sil_extra}} class.
#' @param fill_by defines the color coding of the bar fill color.
#' For `cluster`, the bars are colored after cluster assignment of the
#' observations (default).
#' For `neighbor`, the bar color codes for the nearest neighbor cluster.
#' For `value`, the bar color codes for the silhouette width.
#' For `sign`, the bar color represents the sign of the silhouette width.
#' @param bar_color color of the bar line.
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
#' @export plot.sil_extra
#' @export

  plot.sil_extra <- function(x,
                             fill_by = c('cluster', 'neighbor', 'value', 'sign'),
                             show_stats = TRUE,
                             signif_digits = 2,
                             cust_theme = ggplot2::theme_classic(),
                             bar_color = 'black', ...) {

    ## entry control -----

    stopifnot(is_sil_extra(x))

    fill_by <- match.arg(fill_by[1],
                         c('cluster', 'neighbor', 'value', 'sign'))

    if(!inherits(cust_theme, 'theme')) {

      stop("'cust_theme' hast to be a vaild ggplot2 'theme' object.",
           call. = FALSE)

    }

    x <- mutate(x,
                sil_sign = ifelse(sil_width > 0,
                                  'positive',
                                  ifelse(sil_width < 0,
                                         'negative', 'zero')),
                sil_sign = factor(sil_sign, c('positive', 'zero', 'negative')))

    stopifnot(is.logical(show_stats))
    stopifnot(is.numeric(signif_digits))

    signif_digits <- as.integer(signif_digits)

    sil_width <- NULL
    sil_sign <- NULL
    n <- NULL
    frac_misclassified <- NULL
    clust_id <- NULL
    observation <- NULL
    neighbor_id <- NULL

    ## fill scales -------

    fill_scales <-
      switch(fill_by,
             value = ggplot2::scale_fill_gradient2(low = 'steelblue',
                                                   mid = 'white',
                                                   high = 'firebrick',
                                                   midpoint = mean(range(x$sil_width)),
                                                   name = 'Silhouette width'),
             sign = ggplot2::scale_fill_manual(values = c(negative = 'steelblue',
                                                          zero = 'white',
                                                          positive = 'firebrick'),
                                               name = 'Silhouette sign'))

    ## labeller ----------

    if(!show_stats) {

      stat_labs <- 'label_value'

    } else {

      plot_lab <- NULL

      stats <-
        mutate(summary(x),
               plot_lab = paste0('total: n = ', n,
                                 '\nnegative fraction = ',
                                 signif(frac_misclassified, signif_digits),
                                 '\navg(s) = ', signif(mean, signif_digits)),
               plot_lab = paste(clust_id, plot_lab, sep = '\n'))

      stat_labs <- set_names(stats$plot_lab, stats$clust_id)

    }

    ## plots -------

    bar_plot <-
      switch(fill_by,
             cluster = ggplot(x,
                              aes(x = sil_width,
                                  y = reorder(observation, sil_width),
                                  fill = clust_id)),
             neighbor = ggplot(x,
                               aes(x = sil_width,
                                   y = reorder(observation, sil_width),
                                   fill = neighbor_id)),
             value = ggplot(x,
                            aes(x = sil_width,
                                y = reorder(observation, sil_width),
                                fill = sil_width)),
             sign = ggplot(x,
                           aes(x = sil_width,
                               y = reorder(observation, sil_width),
                               fill = sil_sign)))

    bar_plot <- bar_plot +
      ggplot2::geom_bar(stat = 'identity',
                        color = bar_color, ...) +
      ggplot2::geom_vline(xintercept = 0,
                          linetype = 'dashed') +
      cust_theme +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank()) +
      ggplot2::facet_grid(clust_id ~ .,
                          scales = 'free',
                          space = 'free',
                          labeller = ggplot2::as_labeller(stat_labs)) +
      ggplot2::labs(title = 'Silhouette width',
                    y = 'observation',
                    x = 'Silhouette width')

    if(fill_by %in% c('cluster', 'neighbor')) {

      bar_plot <- bar_plot +
        ggplot2::labs(fill = if(fill_by == 'cluster') 'Cluster' else 'Nearest\nneighbor\ncluster')

    } else {

      bar_plot <- bar_plot +
        fill_scales

    }

    bar_plot

  }

# END ------
