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

#' Plot features of a `red_analysis` object.
#'
#' @description
#' Plots the component table, loadings table - in both cases the
#' first two components/dimensions in form of scatter plots -  or generates
#' a scree plot of the variance percentages associated with
#' the components/dimensions.
#'
#' @details
#' The loadings table plot is available only for the PCA and factor
#' analysis `red_analysis` objects.
#' For `red_analysis` objects created with `clust_analysis` objects, scatter
#' plots of the  scores/components table/layout can convey
#' the cluster assignment information coded by the point or bar color
#' (`label_clust = TRUE`).
#'
#' @param x a `red_analysis` object, created with \code{\link{reduce_data}} or
#' `components()` called for clustering analyses.
#' @param type plot type:
#'
#' * 'component_tbl' or 'score' present the scores (layout) for particular
#' observations in a scatter plot. If `label_clust = TRUE`, point color codes
#' for the optional cluster assignment information.
#'
#' * 'loadings' plot the variable PCA loadings as a scatter plot.
#'
#' * 'scree' plots the percentage of component's variances as a line plot.
#'
#' * 'neighborhood' plots fractions of nearest neighbors of the data points
#' preserved in the neighborhood in the reduced layout. See: \code{\link{np}}.
#'
#' @param label_points logical, should the variable names be displayed in the
#' plot? Valid only for the PCA loadings plot.
#' @param label_clust logical, should the cluster assignment (if available) be
#' coded be the point color?
#' @param cust_theme a ggplot plot theme.
#' @param segment_color color of the lines presented in the PCA loading plot.
#' @param ... extra arguments passed to \code{\link{plot_point}}
#' ('component_tbl', 'score', 'loadings') or \code{\link{np}} ('neighborhood').
#'
#' @return a `ggplot` object.
#'
#' @export plot.red_analysis
#' @export

  plot.red_analysis <- function(x,
                                type = c('component_tbl',
                                         'scores',
                                         'loadings',
                                         'scree',
                                         'neighborhood'),
                                label_points = TRUE,
                                label_clust = FALSE,
                                cust_theme = ggplot2::theme_classic(),
                                segment_color = 'steelblue', ...) {

    ## entry control --------

    stopifnot(is_red_analysis(x))
    stopifnot(inherits(cust_theme, 'theme'))
    stopifnot(is.logical(label_points))

    type <- match.arg(type[1],
                      c('component_tbl', 'scores',
                        'loadings', 'scree',
                        'neighborhood'))

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

    ## neighborhood preservation plots --------

    if(type == 'neighborhood') {

      np_res <- np(x, ...)

      return(plot(np_res))

    }

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

#' @rdname plot.red_analysis
#' @export

  plot.clust_red <- function(x,
                             type = c('component_tbl',
                                      'scores',
                                      'loadings',
                                      'scree',
                                      'neighborhood'),
                             label_points = TRUE,
                             label_clust = TRUE,
                             cust_theme = ggplot2::theme_classic(),
                             segment_color = 'steelblue', ...) {

    ## entry control is done be the superclass method -------

    stopifnot(is.logical(label_clust))

    type <- match.arg(type[1],
                      c('component_tbl',
                        'scores',
                        'loadings',
                        'scree',
                        'neighborhood'))

    ## calls for the superclass method for plots without cluster assignment ------

    if(!label_clust) return(NextMethod())

    if(!type %in% c('component_tbl', 'scores', 'neighborhood')) {

      return(NextMethod())

    }

    ## numbers of observations and variables ------

    plot_n <- nobs(x)

    plot_tag <- paste0('\nObservations: n = ',
                       plot_n$observations,
                       '\nVariables: n = ',
                       plot_n$variables)

    ## scatter plots with the cluster assignment information ------

    if(type %in% c('scores', 'component_tbl')) {

      plot_data <- extract(x, 'scores')

      if(!all(c('comp_1', 'comp_2') %in% names(plot_data))) {

        stop(paste('Atempt to plot a 1-dimensional reduction analysis result.',
                   'Adjust kdim?'),
             call. = FALSE)

      }

      sdevs <- var(x)

      if(x$red_fun == 'pca') {

        ax_labs <- map2(c('PC1', 'PC2'),
                        signif(sdevs$perc_var[1:2], 3),
                        ~paste0(.x, ', ', .y, '%'))

      } else {

        ax_labs <- map2(c('Dim 1', 'Dim 2'),
                        signif(sdevs$perc_var[1:2], 3),
                        ~paste0(.x, ', ', .y, '%'))

      }

      point_plot <-
        plot_point(data = plot_data,
                   x_var = 'comp_1',
                   y_var = 'comp_2',
                   fill_var = 'clust_id',
                   plot_title = switch(x$red_fun ,
                                       pca = 'PCA',
                                       mds = 'MDS',
                                       umap = 'UMAP'),
                   plot_tag = plot_tag,
                   x_lab = ax_labs[[1]],
                   y_lab = ax_labs[[2]],
                   cust_theme = cust_theme,
                   fill_lab = 'Cluster ID', ...)

      return(point_plot)

    }

  }


# prediction ------

#' Project new data onto a reduction analysis layout.
#'
#' @description
#' Predicts reduction analysis scores for a new piece of data with a reduction
#' algorithm-specific methodology.
#'
#' @details
#' Currently implemented only for PCA and UMAP reduction analysis objects.
#' The method employs internally the function `predict()` from the
#' packages `stats` and `umap`. For reduction analysis objects
#' created with other methods, `NULL` is returned with a warning.
#' Of note, specifically for UMAP reduction analysis objects, the method will
#' work only for few basic distances implemented by the `umap` package by
#' default (Euclidean, Manhattan and cosine).
#'
#' @return an object of the \code{\link{red_analysis}} class.
#'
#' @param object a `red_analysis` object, see Details.
#' @param newdata a numeric data frame or a numeric matrix with the new
#' data set.
#' @param ... extra arguments, currently none.
#'
#' @export predict.red_analysis
#' @export

   predict.red_analysis <- function(object, newdata, ...) {

     stopifnot(is_red_analysis(object))

     check_numeric(newdata)

     if(!object$red_fun %in% c('umap', 'pca')) {

       warning(paste("'predict()' is implemented only for redction analysis",
                     "objects created with the PCA and UMAP algorithms."),
               call. = FALSE)

       return(NULL)

     }

     data <- enexpr(newdata)

     ## predictions -------

     preds <- predict(object$red_obj, newdata)

     preds <- set_names(as.data.frame(preds),
                        paste0('comp_', 1:ncol(preds)))

     preds <- rownames_to_column(preds, 'observation')

     preds <- list(data = quo(!!data),
                   red_obj = NULL,
                   red_fun = 'predicton',
                   dist_method = object$dist_method,
                   component_tbl = preds,
                   loadings = NULL)

     red_analysis(preds)

   }

# END ------
