# Heat maps of clustering variable levels.

#' Plot levels of clustering features in a heat map.
#'
#' @description Generates a heat map of the clustering features, cluster
#' assignment is indicated by the plot faceting.
#'
#' @details
#' `plot_clust_hm()` is a S3 generic function.
#' Note that it is not possible to visualize clustering variable levels for
#' `clust_analysis` objects generated with user-provided dissimilarity matrices.
#' In such cases, `NULL` is returned with a warning.
#'
#' @param x_object a `clust_analysis` or `combi_analysis` object, specifies
#' clustering of the observations.
#' @param y_object a `clust_analysis` or `combi_analysis` object, specifies
#' clustering of the features, an optional parameter. Ignored in case of
#' multi-layer SOM.
#' @param line_color color of the line around heat map tiles.
#' @param plot_title plot title. If `NULL`, the plots generated for multi-layer
#' SOM analyses will be named after the data layers.
#' @param plot_subtitle plot subtitle.
#' @param x_lab x axis title.
#' @param fill_lab fill scale title.
#' @param cust_theme a ggplot theme.
#' @param discrete_fill logical, force a discrete fill scale?
#' @param ... extra arguments passed to methods.
#'
#' @return a `ggplot` object (single-layer analysis) or a list of
#' `ggplot` objects (multi-layer cases).
#'
#' @export

  plot_clust_hm <- function(x_object, ...) UseMethod('plot_clust_hm')

#' @rdname plot_clust_hm
#' @export plot_clust_hm.clust_analysis
#' @export

  plot_clust_hm.clust_analysis <- function(x_object,
                                           y_object = NULL,
                                           line_color = NA,
                                           plot_title = NULL,
                                           plot_subtitle = NULL,
                                           x_lab = 'Sample',
                                           fill_lab = 'Feature level',
                                           cust_theme = ggplot2::theme_classic(),
                                           discrete_fill = FALSE, ...) {

    ## entry control -----

    stopifnot(is_clust_analysis(x_object))

    if(!is.null(y_object)) {

      if(!is_clust_analysis(y_object) & !is_combi_analysis(y_object)) {

        stop(paste("'y_object' has to be of clust_analysis',
                   'or combi_analysis class."),
             call. = FALSE)

      }

      if(is_multi_layer(y_object)) {

        stop("Multi-layer SOM 'y_objects' are not implemented.",
             call. = FALSE)

      }

    }

    stopifnot(inherits(cust_theme, 'theme'))
    stopifnot(is.logical(discrete_fill))

    value <- NULL
    feature <- NULL
    sample_node <- NULL

    ## observation numbers ------

    n_numbers <- ngroups(x_object)

    plot_tag <-
      map2_chr(n_numbers[[1]], n_numbers[[2]],
               paste, sep = ': n = ')

    plot_tag <- paste(plot_tag, collapse = ', ')

    ## plotting -----

    if(!is_multi_layer(x_object)) {

      hm_plot <- ft_hm_single(x_object, y_object, line_color = line_color)

      if(!discrete_fill) {

        hm_plot <- hm_plot +
          ggplot2::scale_fill_gradient2(low = 'steelblue',
                                        mid = 'black',
                                        high = 'firebrick')

      }

      hm_plot <- hm_plot +
        cust_theme +
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                       plot.tag.position = 'bottom') +
        ggplot2::labs(title = plot_title,
                      subtitle = plot_subtitle,
                      x = x_lab,
                      fill = fill_lab,
                      tag = plot_tag)

    } else {

      hm_plot <- ft_hm_multi(x_object, line_color = line_color)

      if(!discrete_fill) {

        hm_plot <-
          map(hm_plot,
              ~.x +
                ggplot2::scale_fill_gradient2(low = 'steelblue',
                                              mid = 'black',
                                              high = 'firebrick'))

      }

      hm_plot <-
        map(hm_plot,
            ~.x +
              cust_theme +
              ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                             plot.tag.position = 'bottom') +
              ggplot2::labs(subtitle = plot_subtitle,
                            x = x_lab,
                            fill = fill_lab,
                            tag = plot_tag))

      if(is.null(plot_title)) {

        hm_plot <-
          map2(hm_plot, names(hm_plot),
               ~.x +
                 labs(title = .y))

      }

    }

    return(hm_plot)

  }

#' @rdname plot_clust_hm
#' @export

  plot_clust_hm.min_analysis <- function(x_object, ...) {

    warning(paste('It is not possible to visualize levels of clustering',
                  'features for cluster analysis objects generated with',
                  'user-provided dissimilarity matrices.'),
            call. = FALSE)

    return(NULL)

  }

#' @rdname plot_clust_hm
#' @export plot_clust_hm.combi_analysis
#' @export

  plot_clust_hm.combi_analysis <- function(x_object,
                                           y_object = NULL,
                                           line_color = NA,
                                           plot_title = NULL,
                                           plot_subtitle = NULL,
                                           x_lab = 'Sample',
                                           fill_lab = 'Feature level',
                                           cust_theme = ggplot2::theme_classic(),
                                           discrete_fill = FALSE, ...) {

    ## entry control -------

    stopifnot(is_combi_analysis(x_object))

    if(!is.null(y_object)) {

      if(!is_clust_analysis(y_object) & !is_combi_analysis(y_object)) {

        stop(paste("'y_object' has to be of clust_analysis',
                   'or combi_analysis class."),
             call. = FALSE)

      }

      if(is_multi_layer(y_object)) {

        stop("Multi-layer SOM 'y_objects' are not implemented.",
             call. = FALSE)

      }

    }

    stopifnot(inherits(cust_theme, 'theme'))
    stopifnot(is.logical(discrete_fill))

    value <- NULL
    feature <- NULL
    sample_node <- NULL

    ## observation numbers ------

    n_numbers <- ngroups(x_object)$final

    plot_tag <-
      map2_chr(n_numbers[[1]], n_numbers[[2]],
               paste, sep = ': n = ')

    plot_tag <- paste(plot_tag, collapse = ', ')

    ## plotting -------

    hm_plot <- ft_hm_single(x_object, y_object, line_color = line_color)

    if(!discrete_fill) {

      hm_plot <- hm_plot +
        ggplot2::scale_fill_gradient2(low = 'steelblue',
                                      mid = 'black',
                                      high = 'firebrick')

    }

    hm_plot <- hm_plot +
      cust_theme +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                     plot.tag.position = 'bottom') +
      ggplot2::labs(title = plot_title,
                    subtitle = plot_subtitle,
                    x = x_lab,
                    fill = fill_lab,
                    tag = plot_tag)


    hm_plot

  }

#' @rdname plot_clust_hm
#' @export

  plot_clust_hm.umatrix_analysis <- function(x_object,
                                             line_color = NA,
                                             plot_title = NULL,
                                             plot_subtitle = NULL,
                                             x_lab = 'Sample',
                                             fill_lab = 'Feature level',
                                             cust_theme = ggplot2::theme_classic(),
                                             discrete_fill = FALSE, ...) {

    ## entry check ------

    stopifnot(is_combi_analysis(x_object))
    stopifnot(is_umatrix_analysis(x_object))

    stopifnot(inherits(cust_theme, 'theme'))

    stopifnot(is.logical(discrete_fill))

    ## n numbers -------

    n_numbers <- ngroups(x_object)$final

    plot_tag <-
      map2_chr(n_numbers[[1]], n_numbers[[2]],
               paste, sep = ': n = ')

    plot_tag <- paste(plot_tag, collapse = ', ')

    ## plotting -------

    plot_lst <- ft_hm_multi(x_object = x_object,
                            line_color = line_color,
                            discrete_fill = discrete_fill)

    if(!discrete_fill) {

      plot_lst <-
        map(plot_lst,
            ~.x +
              scale_fill_gradient2(low = 'steelblue',
                                   mid = 'black',
                                   high = 'firebrick'))

    }

    map2(plot_lst, names(plot_lst),
         ~.x +
           ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                          plot.tag.position = 'bottom') +
           ggplot2::labs(title = .y,
                         subtitle = plot_subtitle,
                         tag = plot_tag,
                         x = x_lab,
                         fill = fill_lab))

  }

# END ------
