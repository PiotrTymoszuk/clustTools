# S3 OOP interface for the `tuner` class

# Extraction -------

#' Extract cluster assignment predictions.
#'
#' @description
#' The function extracts the `clust_analysis` object, quality statistics or
#' the best combination of the tuning parameters.
#'
#' @param x a `tuner` object.
#' @param type feature to be extracted from the object:
#'
#' * `clust_object` or `analysis` returns the `clust_analysis` object generated
#' with the best combination of the tuning parameters
#'
#' * `stats` extracts quality statistics for combinations of the tuning
#' parameters
#'
#' * `criteria` returns a data frame with criteria of selection of the best
#' combination of the tuning parameters
#'
#' * `best_tune` extracts a data frame with the best values of the tuning
#' parameters
#'
#' @param ... extra arguments, currently none.
#'
#' @return a `clust_analysis` object.
#'
#' @export

  extract.tuner <- function(x,
                            type = c('clust_object',
                                     'analysis',
                                     'stats',
                                     'criteria',
                                     'best_tune'), ...) {

    stopifnot(is_tuner(x))

    type <- match.arg(type[1],
                      c('clust_object',
                        'analysis',
                        'stats',
                        'criteria',
                        'best_tune'))

    switch(type,
           clust_object = x$analysis,
           analysis = x$analysis,
           stats = x$stats,
           criteria = x$tune_criteria,
           best_tune = x$best_tune)

  }

# Summary -------

#' Summary of quality statistics.
#'
#' @description
#' The `summary()` method called for `tuner` class objects extracts cluster
#' assignment prediction statistics for subsequent combinations of the tuning
#' parameters from the object.
#'
#' @param object a `tuner` object.
#' @param ... extra arguments, currently none.
#'
#' @return a data frame with the following columns:
#'
#' * `sil_width`: silhouette width
#'
#' * `frac_misclassified`: fraction of observations with negative silhouette
#' widths suggestive of misclassification
#'
#' * `frac_var`: fraction of explained clustering variance
#'
#' * `frac_np`: fraction of preserved nearest neighbors
#'
#' * columns named after names of the tuning parameters and containing their
#' values
#'
#' @export

  summary.tuner <- function(object, ...) object$stats

# Plotting -------

#' Plot cluster quality statistic values for tuning.
#'
#' @description
#' The function plots quality statistic values (such as explained clustering
#' variance, silhouette width or mean neighborhood preservation) for
#' combinations of the tuning parameters.
#'
#' @param x a `tuner` object.
#' @param cust_theme a custom `ggplot` theme.
#' @param line_alpha alpha of the lines, applies only to plots of
#' regularization paths.
#' @param point_size size of the data points.
#' @param ... extra arguments, currently none.
#'
#' @return a list of `ggplot` graphic objects.
#'
#' @export plot.tuner
#' @export

  plot.tuner <- function(x,
                         cust_theme = ggplot2::theme_classic(),
                         line_alpha = 1,
                         point_size = 2, ...) {

    ## entry control ------

    stopifnot(is_tuner(x))

    if(!inherits(cust_theme, 'theme')) {

      stop("'cust_theme' has to be a valid ggplot theme.",
           call. = FALSE)

    }

    ## variable declaration, compatibility with CHECK ------

    n_active_vars <- NULL
    value <- NULL
    statistic <- NULL
    labs <- NULL
    lambda <- NULL
    distance <- NULL
    variable <- NULL

    ## plot data and metadata --------

    stat_vars <- c('sil_width', 'frac_misclassified', 'frac_var', 'frac_np')

    stat_titles <-
      set_names(c('Cluster separation',
                  'Potential misclassification',
                  'Explained clustering variance',
                  'Neighborhood preservation'),
                stat_vars)

    stat_labels <-
      set_names(c('mean silhouette width',
                  'fraction of negative silhouette widths',
                  'fraction of explained clustering variance',
                  'fraction of preserved nearest neighbors',
                  'fraction of active variables'),
                c(stat_vars, 'n_active_vars'))

    stat_colors <-
      set_names(c('coral3',
                  'steelblue3',
                  'darkolivegreen',
                  'plum4',
                  'firebrick'),
                c(stat_vars, 'n_active_vars'))

    data_label <- switch(x$dataset,
                         train = 'training',
                         cv = 'CV')

    criteria <- paste(stat_labels[names(x$tune_criteria)], collapse = ', ')

    best_label <- map_dbl(x$best_tune, signif, 2)

    best_label <- map2_chr(names(best_label), best_label,
                           paste, sep = ' = ')

    best_label <- paste(best_label, collapse = ', ')

    plot_subtitle <- paste0('dataset: ', data_label,
                            '\ncriteria: ', criteria,
                            '\nbest tune:', best_label)

    ## plotting table -------

    stats <- summary(x)

    if(x$fun == 'tune_htk') {

      stat_vars <- c(stat_vars, 'n_active_vars')

      stats <- mutate(stats,
                      n_active_vars = n_active_vars/length(x$clust_vars))

    }

    stats <-
      tidyr::pivot_longer(stats,
                          cols = dplyr::all_of(stat_vars),
                          names_to = 'statistic',
                          values_to = 'value')

    ## plotting: base plots --------

    if(x$fun %in% c('prediter', 'tune_htk')) {

      base_var <- x$tune_params[1]

      x_intercept <- x$best_tune[[base_var]][1]

      stat_plot <-
        ggplot(stats,
               aes(x = .data[[base_var]],
                   y = value,
                   color = statistic)) +
        ggplot2::geom_path() +
        ggplot2::geom_point(shape = 16,
                            size = point_size) +
        ggplot2::geom_vline(xintercept = x_intercept,
                            linetype = 'dashed') +
        ggplot2::scale_color_manual(values = stat_colors,
                                    labels = stat_labels,
                                    name = 'statistic') +
        cust_theme +
        labs(title = paste('Tuning of', base_var),
             subtitle = plot_subtitle,
             x = base_var,
             y = 'statistic value')

    }

    ## plotting: regularization paths ------

    if(x$fun %in% c('tune_htk')) {

      x_intercept <- x$best_tune[['lambda']][1]

      reg_data <-
        tidyr::pivot_longer(x$center_distances,
                            cols = dplyr::all_of(x$clust_vars),
                            names_to = 'variable',
                            values_to = 'distance')

      reg_plot <-
        ggplot(reg_data,
               aes(x = lambda,
                   y = distance,
                   color = variable)) +
        ggplot2::geom_line(alpha = line_alpha) +
        ggplot2::geom_vline(xintercept = x_intercept,
                            linetype = 'dashed') +
        cust_theme +
        ggplot2::labs(title = 'Regularization paths',
                      subtitle = plot_subtitle,
                      x = 'lambda',
                      y = 'cluster center - mean distance')

    }

    ## output ----

    if(x$fun == 'prediter') {

      return(stat_plot)

    } else if(x$fun == 'tune_htk') {

      return(list(statistics = stat_plot,
                  regularization = reg_plot))

    }

  }

# END ------
