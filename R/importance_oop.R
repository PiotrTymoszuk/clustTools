# S3 OOP methods for the importance class.

# Plotting ------

#' Plot feature importance as a scatter or bar plot.
#'
#' @description Generates a bar, scatter or box plot with the importance
#' statistic for the clustering variables.
#' The importance statistic is the difference in the
#' clustering variance fraction between the original clustering structure and
#' the clustering objects with the given variable reshuffled randomly.
#' @return a ggplot bar or scatter plot.
#' @param x an `importance` object.
#' @param type type of the plot: scatter or bar. Defaults to scatter.
#' This parameter is silently ignored, if evaluation of the importance was done
#' in multiple iterations
#' (e.g. `n_iter` set to > 1 in \code{\link{impact.clust_analysis}}).
#' In such cases, a box plot of importance metrics obtained in algorithm
#' iteration is generated.
#' @param fill_color fill color for the bars or boxes.
#' @param point_color size of the points, refers only to scatter and box plots.
#' @param point_size size of the points, refers only to scatter and box plots.
#' @param point_alpha alpha of the points, refers only to box plots.
#' @param point_wjitter width of the data point jittering, refers only to
#' box plots.
#' @param point_hjitter height of the data point jittering, refers only to
#' box plots.
#' @param box_alpha alpha of the boxes, refers only to box plots.
#' @param label logical, should the points be labeled with the importance
#' stat value?
#' @param txt_size label text size.
#' @param signif_digits significant digits for rounding of the statistic value.
#' @param plot_title plot title.
#' @param plot_subtitle plot subtitle.
#' @param plot_tag plot tag.
#' @param cust_theme a ggplot theme.
#' @param ... extra arguments, currently none.
#' @export plot.importance
#' @export

  plot.importance <- function(x,
                              type = c('scatter', 'bar'),
                              fill_color = 'cornsilk3',
                              point_color = fill_color,
                              point_size = 2,
                              point_alpha = 0.5,
                              point_wjitter = 0,
                              point_hjitter = 0.1,
                              box_alpha = 0.25,
                              label = TRUE,
                              txt_size = 2.75,
                              signif_digits = 2,
                              plot_title = NULL,
                              plot_subtitle = NULL,
                              plot_tag = NULL,
                              cust_theme = ggplot2::theme_classic(), ...) {

    ## entry control -------

    stopifnot(is_importance(x))

    type <- match.arg(type[1], c('scatter', 'bar'))

    stopifnot(is.logical(label))
    stopifnot(inherits(cust_theme, 'theme'))

    variable <- NULL
    frac_diff <- NULL

    ## plotting: single iteration ---------

    if(!'run' %in% names(x)) {

      base_plot <-
        ggplot(filter(x, variable != 'data'),
               aes(x = frac_diff,
                   y = reorder(variable, frac_diff))) +
        cust_theme +
        ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
        ggplot2::labs(title = plot_title,
                      subtitle = plot_subtitle,
                      tag = plot_tag,
                      x = expression(Delta*' clustering variance'))

      if(type == 'scatter') {

        base_plot <- base_plot +
          ggplot2::geom_point(size = point_size,
                              shape = 21,
                              fill = point_color)

      } else {

        base_plot <- base_plot +
          ggplot2::geom_bar(stat = 'identity',
                            color = 'black',
                            fill = fill_color)

      }

      if(label) {

        base_plot <- base_plot +
          ggplot2::geom_text(aes(label = signif(frac_diff,
                                                signif_digits)),
                             size = txt_size,
                             vjust = if(type == 'bar') 0.5 else -0.9,
                             hjust = if(type == 'bar') 1.2 else 0.6)

      }

      return(base_plot)

    }

    ## multiple iterations: a box plot with single runs visualized as points ------

    ggplot(filter(x, variable != 'data'),
           aes(x = frac_diff,
               y = reorder(variable, frac_diff))) +
      ggplot2::geom_boxplot(alpha = box_alpha,
                            fill = fill_color,
                            color = 'black',
                            outlier.color = NA) +
      ggplot2::geom_point(shape = 21,
                          size = point_size,
                          alpha = point_alpha,
                          color = 'black',
                          fill = point_color,
                          position = ggplot2::position_jitter(width = point_wjitter,
                                                              height = point_hjitter)) +
      cust_theme +
      ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
      ggplot2::labs(title = plot_title,
                    subtitle = plot_subtitle,
                    tag = plot_tag,
                    x = expression(Delta*' clustering variance'))

  }

# Summary --------

#' Importance statistic summary.
#'
#' @description
#' If the permutation importance analysis for clustering variables was done
#' in multiple iterations
#' (e.g. `n_iter` set to > 1 in \code{\link{impact.clust_analysis}}),
#' number of iterations, mean, SD, median, interquartile range
#' and range of the difference in
#' clustering variance for each clustering variable is computed.
#' Otherwise, a plain data frame with importance statistics is returned.
#'
#' @param object an `importance` class object.
#' @param ... extra arguments, currently none.
#'
#' @return a data frame with importance metrics.
#'
#' @export summary.importance
#' @export

  summary.importance <- function(object, ...) {

    stopifnot(is_importance(object))

    if(!'run' %in% names(object)) return(object)

    object <- filter(object, variable != 'data')

    n_iter <- length(unique(object$run))

    variable <- NULL

    stats <- dplyr::group_by(object, variable)

    frac_diff <- NULL
    mean_imp <- NULL
    sd_imp <- NULL
    median_imp <- NULL
    q25_imp <- NULL
    q75_imp <- NULL
    min_imp <- NULL
    max_imp <- NULL

    stats <-
      dplyr::summarise(stats,
                       mean_imp = mean(frac_diff, na.rm = TRUE),
                       sd_imp = stats::sd(frac_diff, na.rm = TRUE),
                       median_imp = stats::median(frac_diff, na.rm = TRUE),
                       q25_imp = stats::quantile(frac_diff, 0.25, na.rm = TRUE),
                       q75_imp = stats::quantile(frac_diff, 0.75, na.rm = TRUE),
                       min_imp = min(frac_diff, na.rm = TRUE),
                       max_imp = max(frac_diff, na.rm = TRUE))

    stats <- set_names(stats,
                       stringi::stri_replace(names(stats),
                                             fixed = '_imp', replacement = ''))

    mutate(stats,
           n_iter = n_iter)

  }

# END -----
