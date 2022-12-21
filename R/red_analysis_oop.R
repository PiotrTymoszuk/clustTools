# Provides methods for the red_analysis class.

# general extractor -----

#' Extract features of a red_analysis object.
#'
#' @description A general extractor method for accessing properties and features
#' of a red_analysis object.
#' @param x a red_analysis object, created with \code{\link{reduce_data}}.
#' @param type the feature name:
#' 'component_tbl' or 'scores' return the component pr score tables for the
#' observations, 'loadings' retrieves the table of variable loadings (PCA only),
#' 'data' returns the input data frame, 'sdev' returns standard deviations
#' associated with the components (PCA only), 'object' returns the PCA or UMAP
#' object.
#' @param ... extra arguments, currently none.
#' @return the requested feature/property.
#' @export extract.red_analysis
#' @export

  extract.red_analysis <- function(x,
                                   type = c('component_tbl',
                                            'scores',
                                            'loadings',
                                            'data',
                                            'sdev',
                                            'object'), ...) {

    stopifnot(all(class(x) == 'red_analysis'))

    type <- match.arg(type[1],
                      c('component_tbl',
                        'scores',
                        'loadings',
                        'data',
                        'sdev',
                        'object'))

    if(type != 'sdev') {

      switch(type,
             component_tbl = x$component_tbl,
             scores = x$component_tbl,
             loadings = x$loadings,
             data = rlang::eval_tidy(x$data),
             object = x$red_obj)

    } else if(x$red_fun == 'pca') {

      tibble::tibble(component = 1:length(x$red_obj$sdev),
                     sdev = x$red_obj$sdev,
                     perc_sdev = x$red_obj$sdev/sum(x$red_obj$sdev) * 100,
                     var = x$red_obj$sdev^2,
                     perc_var = x$red_obj$sdev^2/sum(x$red_obj$sdev^2) * 100)

    } else if(x$red_fun == 'fa') {

      p <- nrow(x$red_obj$loadings)

      sum_sq <- colSums(x$red_obj$loadings^2)

      tibble::tibble(component = 1:ncol(x$red_obj$loadings),
                     sdev = sqrt(sum_sq/p),
                     perc_sdev = sqrt(sum_sq/p) * 100,
                     var = sum_sq/p,
                     perc_var = sum_sq/p * 100)

    } else {

      score_tbl <- x$component_tbl

      dims <- stringi::stri_extract(names(score_tbl), regex = 'comp_\\d+')

      dims <- dims[!is.na(dims)]

      variances <- purrr::map_dbl(score_tbl[dims], var)

      std_devs <- purrr::map_dbl(score_tbl[dims], sd)

      tibble::tibble(component = 1:length(dims),
                     sdev = std_devs,
                     perc_sdev = std_devs/sum(std_devs)*100,
                     var = variances,
                     perc_var = variances/sum(variances)*100)

    }

  }

# appearance and class testing -----

#' Test for the red_analysis class.
#'
#' @description Tests whwtrer the object is a red_analysis class.
#' @param x an object to test.
#' @return a logical value.
#' @export

  is_red_analysis <- function(x) {

    all(class(x) == 'red_analysis')

  }

#' Print a red_analysis object.
#'
#' @description Prints the red_analysis object.
#' @param x a red_analysis object.
#' @param ... extra arguments, currently none.
#' @export

  print.red_analysis <- function(x, ...) {

    stopifnot(all(class(x) == 'red_analysis'))

    cat(paste0(toupper(x$red_fun), ' reduction analysis object.'))
    cat('\nComponents:\n')
    print(tibble::as_tibble(x$component_tbl))


  }

# data, observation number -----

#' Model frame/data of a red_analysis object.
#'
#' @description Extracts the data frame used for construction of a red_analysis
#' object.
#' @param formula a red_analysis object.
#' @param ... extra arguments, currently none.
#' @export

  model.frame.red_analysis <- function(formula, ...) {

    stopifnot(all(class(formula) == 'red_analysis'))

    rlang::eval_tidy(formula$data)

  }

#' Get the number of observations for a red_analysis object.
#'
#' @description Returns the number of complete observations used for
#' construction of a red_analysis object.
#' @param object a red_analysis object.
#' @param ... extra arguments, currently none.
#' @export

  nobs.red_analysis <- function(object, ...) {

    stopifnot(all(class(object) == 'red_analysis'))

    clustTools:::get_data_dim(rlang::eval_tidy(object$data))

  }

#' component variance.
#'
#' @description Computes the standard variations and variances associated with
#' the dimensions/component of a red_analysis object.
#' @param x a red_analysis object.
#' @param ... extra arguments, currently none.
#' @export

  var.red_analysis <- function(x, ...) {

    stopifnot(all(class(x) == 'red_analysis'))

    clustTools::extract(x, type = 'sdev')

  }

#' Summary of a red_analysis object.
#'
#' @description Summary method for the red_analysis class. Prints the summary
#' of the reduction analysis object.
#' @param object a red_analysis object.
#' @param ... extra arguments, currently none.
#' @export

  summary.red_analysis <- function(object, ...) {

    stopifnot(all(class(object) == 'red_analysis'))

    summary(object$red_obj)

  }

# plotting ----

#' Plot features of a red_analysis object.
#'
#' @description Plots the component table, loadings table - in both cases the
#' first two components/dimensions in form of scatter plots -  or generates
#' a scree plot of the variance percentages associated with
#' the components/dimensions.
#' @details The loadings table plot is available only for the PCA red_analysis
#' objects.
#' @param x a red_analysis object, created with \code{\link{reduce_data}}.
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

    ## entry control

    stopifnot(all(class(x) == 'red_analysis'))
    stopifnot(any(class(cust_theme) == 'theme'))
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

    ## plot meta

    plot_n <- clustTools::nobs(x)

    plot_tag <- paste0('\nObservations: n = ',
                       plot_n$observations,
                       '\nVariables: n = ',
                       plot_n$variables)

    sdevs <- clustTools::extract(x, 'sdev')

    if(x$red_fun %in% c('mds', 'umap', 'fa')) {

      ax_labs <- purrr::map2(c('Dim 1', 'Dim 2'),
                             signif(sdevs$perc_var[1:2], 3),
                             ~paste0(.x, ', ', .y, '%'))

    } else {

      ax_labs <- purrr::map2(c('PC1', 'PC2'),
                             signif(sdevs$perc_var[1:2], 3),
                             ~paste0(.x, ', ', .y, '%'))

    }

    ## plotting

    if(type %in% c('component_tbl', 'scores')) {

      plot_data <- clustTools::extract(x, 'scores')

      if(is.null(plot_data)) {

        warning('No component table/scores available.', call. = FALSE)

        return(NULL)

      }

      clustTools::plot_point(data = plot_data,
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

    } else if(type == 'loadings') {

      clustTools::plot_point(data = clustTools::extract(x, 'loadings'),
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

    } else {

      sdevs <- dplyr::mutate(sdevs,
                             line_group = 'gr1')

      ggplot2::ggplot(sdevs,
                      ggplot2::aes(x = component,
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

    }

  }
