# Additional plotting functions.

#' Plot levels of clustering features in a heat map.
#'
#' @description Generates a heat map of the clustering features, cluster
#' assignment is indicated by the plot faceting.
#' @param x_object a `clust_analysis` or `combi_analysis` object, specifies
#' clustering of the observations.
#' @param y_object a `clust_analysis` or `combi_analysis` object, specifies
#' clustering of the features. Optional.
#' @param plot_title plot title.
#' @param plot_subtitle plot subtitle.
#' @param x_lab x axis title.
#' @param fill_lab fill scale title.
#' @param cust_theme a ggplot theme.
#' @param discrete_fill logical, force a discrete fill scale?
#' @return a ggplot graph.
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom purrr map2_chr
#' @importFrom rlang set_names
#' @export

  plot_clust_hm <- function(x_object,
                            y_object = NULL,
                            plot_title = NULL,
                            plot_subtitle = NULL,
                            x_lab = 'Sample',
                            fill_lab = 'Feature level',
                            cust_theme = ggplot2::theme_classic(),
                            discrete_fill = FALSE) {

    ## entry control ------

    if(!class(x_object) %in% c('combi_analysis', 'clust_analysis')) {

      stop('The sample cluster object has to be of clust_analysis or combi_analysis class',
           call. = FALSE)

    }

    if(!is.null(y_object)) {

      if(!class(y_object) %in% c('combi_analysis', 'clust_analysis')) {

        stop('The sample cluster object has to be of clust_analysis or combi_analysis class',
             call. = FALSE)

      }

    }

    stopifnot(inherits(cust_theme, 'theme'))
    stopifnot(is.logical(discrete_fill))

    value <- NULL
    feature <- NULL
    sample_node <- NULL

    ## assignment list ---------

    if('node' %in% names(x_object$clust_assignment)) {

      x_ass_names <- c('sample', 'sample_node', 'sample_clust')

    } else {

      x_ass_names <- c('sample', 'sample_clust')

    }

    if(!is.null(y_object)) {

      y_ass_names <- switch(class(y_object)[1],
                            clust_analysis = c('feature',
                                               'feature_clust'),
                            combi_analysis = c('feature',
                                               'feature_node',
                                               'feature_clust'))

      cmm_assignment <- map(list(x_object, y_object),
                            ~.x$clust_assignment)

      cmm_assignment <- map2(cmm_assignment,
                             list(x_ass_names, y_ass_names),
                             set_names)

      cmm_assignment <- set_names(cmm_assignment,
                                  c('sample', 'feature'))

    } else {

      cmm_assignment <- list()

      cmm_assignment$sample <- set_names(x_object$clust_assignment,
                                         x_ass_names)

    }

    ## data in a long format

    data <- switch(class(x_object),
                   clust_analysis = model.frame(x_object),
                   combi_analysis = model.frame(x_object)[[1]])

    data <- as_tibble(data)

    features <- names(data)

    data <- mutate(data,
                   sample = cmm_assignment$sample$sample)

    data <-
      tidyr::pivot_longer(data,
                          cols = all_of(features),
                          names_to = 'feature',
                          values_to = 'value')

    ## joining the data table with the cluster assignment information

    data <- left_join(data,
                      cmm_assignment$sample,
                      by = 'sample')

    if(!is.null(y_object)) {

      if(all(!features %in% cmm_assignment$feature$feature)) {

        stop(paste('Unable to retrieve the cluster assignment information ',
                   'or the clustering variables from the feature_clust_object.'),
             call. = FALSE)

      }

      data <- left_join(data,
                        cmm_assignment$feature,
                        by = 'feature')

    }

    ## plot tag with the sample n numbers

    n_tag <- ngroups(x_object)

    if(is_combi_analysis(x_object)) n_tag <- n_tag$final

    n_tag <- map2_chr(n_tag$clust_id,
                      n_tag$n,
                      ~paste0(.x, ': n = ', .y))

    n_tag <- paste(n_tag, collapse = ', ')

    ## plotting

    if(is_clust_analysis(x_object)) {

      base_plot <-
        ggplot(data,
               aes(x = reorder(sample,
                               as.numeric(factor(value))),
                   y = reorder(feature,
                               as.numeric(factor(value))),
                   fill = if(discrete_fill) factor(value) else value))

    } else {

      base_plot <-
        ggplot(data,
               aes(x = reorder(sample,
                               as.numeric(sample_node)),
                   y = reorder(feature,
                               as.numeric(factor(value))),
                   fill = if(discrete_fill) factor(value) else value))

    }

    if(!is.null(y_object)) {

      base_plot <- base_plot +
        ggplot2::facet_grid(feature_clust ~ sample_clust,
                            scales = 'free',
                            space = 'free')

    } else {

      base_plot <- base_plot +
        ggplot2::facet_grid(. ~ sample_clust,
                            scales = 'free',
                            space = 'free')

    }

    base_plot <- base_plot +
      ggplot2::geom_tile() +
      cust_theme +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     axis.line = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank()) +
      ggplot2::labs(title = plot_title,
                    subtitle = plot_subtitle,
                    tag = n_tag,
                    x = x_lab)

    if(!discrete_fill){

      base_plot <- base_plot +
        ggplot2::scale_fill_gradient2(low = 'steelblue',
                                      mid = 'black',
                                      high = 'firebrick',
                                      midpoint = mean(range(data$value)),
                                      name = fill_lab)

    }

    base_plot

  }

# END ------
