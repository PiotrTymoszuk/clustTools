# Provides specific methods for the red_analysis class.

# variance and summary -----

#' Variance and summary for a red_analysis object.
#'
#' @description Variance associated with the
#' components and statistic summary for `red_analysis` class objects.
#'
#' @param object a `red_analysis` object.
#' @param x a `red_analysis` object.
#' @param ... extra arguments, currently none.
#'
#' @return `var()` returns a data frame with components'
#' variances, `summary()` returns a set of summary statistic specific for
#' the wrapped dimensionality reduction function.
#'
#' @export

  var.red_analysis <- function(x, ...) {

    stopifnot(is_red_analysis(x))

    extract.red_analysis(x, type = 'sdev')

  }

#' @rdname var.red_analysis
#' @export

  summary.red_analysis <- function(object, ...) {

    stopifnot(is_red_analysis(object))

    summary(object$red_obj)

  }

# plotting ----

#' Plot features of a red_analysis object.
#'
#' @description Plots the component table, loadings table - in both cases the
#' first two components/dimensions in form of scatter plots -  or generates
#' a scree plot of the variance percentages associated with
#' the components/dimensions.
#' @details The loadings table plot is available only for the PCA `red_analysis`
#' objects.
#' @param x a `red_analysis` object, created with \code{\link{reduce_data}}.
#' @param type plot type:
#' 'component_tbl' or 'score' present the scores for particular observations in
#' a scatter plot.
#' 'loadings' plot the variable PCA loadings as a scatter plot.
#' 'scree' plots the percentage of component's variances as a line plot.
#' @param label_points logical, should the variable names be displayed in the
#' plot? Valid only for the PCA loadings plot.
#' @param cust_theme a ggplot plot theme.
#' @param segment_color color of the lines presented in the PCA loading plot.
#' @param ... extra arguments passed to \code{\link{plot_point}}.
#' @return a ggplot object.
#' @export plot.red_analysis
#' @export

  plot.red_analysis <- function(x,
                                type = c('component_tbl',
                                         'scores',
                                         'loadings',
                                         'scree'),
                                label_points = TRUE,
                                cust_theme = ggplot2::theme_classic(),
                                segment_color = 'steelblue', ...) {

    ## entry control --------

    stopifnot(is_red_analysis(x))
    stopifnot(inherits(cust_theme, 'theme'))
    stopifnot(is.logical(label_points))

    type <- match.arg(type[1],
                      c('component_tbl', 'scores', 'loadings', 'scree'))

    if(x$red_fun %in% c('mds', 'umap')) {

      if((type == 'loadings')) {

        warning('Loadings are not implemented for the MDS and UMAP dimansionality reduction algotithms',
                call. = FALSE)

        return(NULL)

      }

    }

    component <- NULL
    perc_var <- NULL
    line_group <- NULL

    ## plot meta -------

    plot_n <- nobs(x)

    plot_tag <- paste0('\nObservations: n = ',
                       plot_n$observations,
                       '\nVariables: n = ',
                       plot_n$variables)

    sdevs <- extract(x, 'sdev')

    if(x$red_fun %in% c('mds', 'umap', 'fa')) {

      ax_labs <- map2(c('Dim 1', 'Dim 2'),
                      signif(sdevs$perc_var[1:2], 3),
                      ~paste0(.x, ', ', .y, '%'))

    } else {

      ax_labs <- map2(c('PC1', 'PC2'),
                      signif(sdevs$perc_var[1:2], 3),
                      ~paste0(.x, ', ', .y, '%'))

    }

    ## plotting

    if(type %in% c('component_tbl', 'scores')) {

      plot_data <- extract(x, 'scores')

      if(is.null(plot_data)) {

        warning('No component table/scores available.', call. = FALSE)

        return(NULL)

      }

      point_plot <- plot_point(data = plot_data,
                               x_var = 'comp_1',
                               y_var = 'comp_2',
                               fill_var = NULL,
                               label_var = NULL,
                               plot_title = switch(x$red_fun,
                                                   mds = 'MDS dimensions',
                                                   umap = 'UMAP dimensions',
                                                   pca = 'PCA scores',
                                                   fa = 'FA scores'),
                               plot_tag = plot_tag,
                               x_lab = ax_labs[[1]],
                               y_lab = ax_labs[[2]],
                               cust_theme = cust_theme, ...)

      return(point_plot)

    }

    if(type == 'loadings') {

      point_plot <-
        plot_point(data = extract(x, 'loadings'),
                   x_var = 'comp_1',
                   y_var = 'comp_2',
                   fill_var = NULL,
                   label_var = if(label_points) 'variable' else NULL,
                   plot_title = switch(x$red_fun,
                                       pca = 'PCA loadings',
                                       fa = 'FA loadings'),
                   plot_tag = plot_tag,
                   x_lab = ax_labs[[1]],
                   y_lab = ax_labs[[2]],
                   cust_theme = cust_theme,
                   show_segments = TRUE,
                   segment_color = segment_color, ...)

      return(point_plot)

    }

    if(type == 'scree'){

      sdevs <- mutate(sdevs,
                      line_group = 'gr1')

      scree_plot <-
        ggplot(sdevs,
               aes(x = component,
                   y = perc_var,
                   group = line_group)) +
        ggplot2::geom_line(color = segment_color) +
        cust_theme +
        ggplot2::labs(title = switch(x$red_fun,
                                     mds = 'MDS variance',
                                     umap = 'UMAP variance',
                                     pca = 'PCA variance'),
                      y = '% total variance',
                      x = switch(x$red_fun,
                                 mds = 'Dimension',
                                 umap = 'Dimension',
                                 pca = 'Principal component'),
                      tag = plot_tag) +
        ggplot2::scale_x_continuous(breaks = sdevs$component)

      return(scree_plot)

    }

  }

# END ------
