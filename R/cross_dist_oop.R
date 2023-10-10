# S3 OOP for the cross-dist class.
# Summary -------

#' Summary of cross-distances between clusters.
#'
#' @description
#' Computes summary statistics of homologous or heterologous cross-distances
#' between the clusters.
#'
#' @param object a `cross_dist` object.
#' @param ... extra arguments, currently none.
#'
#' @return a data frame with mean, SD, median, interquartile range, 95% range
#' and range of cross-distances between the cluster pairs.
#'
#' @export summary.cross_dist
#' @export

  summary.cross_dist <- function(object, ...) {

    stopifnot(is_cross_dist(object))

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

    stats <- map(object, function(obj) map_dbl(stat_funs, ~.x(obj)))

    comparison <- NULL
    n_pairs <- NULL
    type <- NULL
    distance <- NULL
    cluster1 <- NULL
    cluster2 <- NULL

    stats <- as_tibble(reduce(stats, rbind))

    stats <-
      mutate(stats,
             comparison = names(object),
             n_pairs = map_dbl(object, length),
             type = attr(object, 'type'),
             distance = attr(object, 'dist_method'),
             cluster1 = stringi::stri_split_fixed(comparison,
                                                  pattern = ' vs ',
                                                  simplify = TRUE)[, 1],
             cluster2 = stringi::stri_split_fixed(comparison,
                                                  pattern = ' vs ',
                                                  simplify = TRUE)[, 2],
             cluster1 = factor(cluster1, attr(object, 'x_levels')),
             cluster2 = factor(cluster2, attr(object, 'x_levels')))

    select(stats,
           cluster1, cluster2,
           type, distance, n_pairs,
           all_of(names(stat_funs)))

  }

# Plotting --------

#' Plots of cross-distances.
#'
#' @description
#' Visualizes pairwise cross-distances as heat maps for observation pairs,
#' heat maps of average cross-distances between the clusters or histograms.
#'
#' @return
#' a `ggplot` graphic, whose elements like themes or fill scales can be easily
#' modified by the user.
#'
#'
#' @param x a `cross_distance` class object.
#' @param type type of the plot:
#' `heat_map` (default) generates a heat map of
#' homologous or heterologous cross-distances for observation pairs with
#' mean distances and 95% ranges of distances,
#' `mean` plots mean distances with 95% ranges as a heat map, and
#' `histogram` generates a faceted panel of cross-distance histograms
#' (for heterologous distances, x object clusters are represented by horizontal
#' facets, y object clusters are represented by vertical facets) .
#' @param reorder logical: should distances in the heat maps be ordered
#' by mean distance? Defaults to FALSE.
#' @param upper should the upper half of the distance heat map be plotted?
#' Defaults to TRUE.
#' @param signif_digits significant digits for mean distances and distance
#' ranges presented in the heat map.
#' @param line_color color of the line around the tile, used only if `type`
#' is set to 'mean'.
#' @param show_txt logical, should the mean distance be presented in the plot?
#' @param txt_size of the mean distance text.
#' @param labeller a \code{\link[ggplot2]{labeller}} object to provide
#' customized labels of the facets of the histogram panel.
#' @param cust_theme a custom ggplot theme.
#' @param ... extra arguments, such as color or number of bins, passed to
#' \code{\link[ggplot2]{geom_histogram}}.
#'
#' @export plot.cross_dist
#' @export

  plot.cross_dist <- function(x,
                              type = c('heat_map', 'mean', 'histogram'),
                              reorder = FALSE,
                              upper = TRUE,
                              signif_digits = 2,
                              line_color = 'black',
                              show_txt = TRUE,
                              txt_size = 2.75,
                              labeller = NULL,
                              cust_theme = ggplot2::theme_classic(), ...) {

    ## entry control ---------

    stopifnot(is_cross_dist(x))

    type <- match.arg(type[1], c('heat_map', 'mean', 'histogram'))

    if(!inherits(cust_theme, 'theme')) {

      stop("'cust_theme' must be a valid ggplot2 theme object.",
           call. = FALSE)

    }

    stopifnot(is.logical(reorder))
    stopifnot(is.logical(upper))
    stopifnot(is.logical(show_txt))

    ## plotting data --------

    if(type %in% c('heat_map', 'histogram')) {

      plot_data <- map(x, as.data.frame)

      y_cols <- map(plot_data, names)

      x_obs <- NULL
      y_obs <- NULL
      distance <- NULL

      plot_data <- map(plot_data, rownames_to_column, 'x_obs')

      plot_data <-
        map2(plot_data, y_cols,
             ~tidyr::pivot_longer(.x,
                                  cols = all_of(.y),
                                  names_to = 'y_obs',
                                  values_to = 'distance'))

      comparison <- NULL
      cluster1 <- NULL
      cluster2 <- NULL
      cluster_prov <- NULL
      obs_prov <- NULL

      plot_data <-
        map2_dfr(plot_data, names(plot_data),
                 ~mutate(.x, comparison = .y))

      plot_data <-
        mutate(plot_data,
               cluster1 = stringi::stri_split_fixed(comparison,
                                                    pattern = ' vs ',
                                                    simplify = TRUE)[, 1],
               cluster2 = stringi::stri_split_fixed(comparison,
                                                    pattern = ' vs ',
                                                    simplify = TRUE)[, 2])

      if(upper & attr(x, 'type') == 'homologous') {

        upper_data <- filter(plot_data,
                             cluster1 != cluster2)

        upper_data <-
          mutate(upper_data,
                 cluster_prov = cluster1,
                 cluster1 = cluster2,
                 cluster2 = cluster_prov,
                 obs_prov = x_obs,
                 x_obs = y_obs,
                 y_obs = obs_prov,
                 x_obs = stringi::stri_replace(x_obs,
                                               regex = '^y_',
                                               replacement = 'x_'),
                 y_obs = stringi::stri_replace(y_obs,
                                               regex = '^x_',
                                               replacement = 'y_'))

        upper_data <- select(upper_data, - cluster_prov, - obs_prov)

        plot_data <- rbind(plot_data, upper_data)

      }

      plot_data <-
        mutate(plot_data,
               cluster1 = factor(cluster1, attr(x, 'x_levels')),
               cluster2 = factor(cluster2, attr(x, 'y_levels')))

    } else {

      cluster1 <- NULL
      cluster2 <- NULL

      plot_data <-
        mutate(summary(x),
               cluster1 = as.character(cluster1),
               cluster2 = as.character(cluster2))

      if(upper & attr(x, 'type') == 'homologous') {

        upper_data <- filter(plot_data,
                             cluster1 != cluster2)

        prov_cluster <- NULL

        upper_data <-
          mutate(upper_data,
                 prov_cluster = cluster1,
                 cluster1 = cluster2,
                 cluster2 = prov_cluster)

        upper_data <- select(upper_data, - prov_cluster)

        plot_data <- rbind(plot_data, upper_data)

      }

      plot_lab <- NULL
      q025 <- NULL
      q975 <- NULL

      plot_data <-
        mutate(plot_data,
               plot_lab = paste0(signif(mean, signif_digits),
                                 '\n[', signif(q025, signif_digits),
                                 ' - ', signif(q975, signif_digits), ']'),
               cluster1 = factor(cluster1, attr(x, 'x_levels')),
               cluster2 = factor(cluster2, attr(x, 'y_levels')))

    }

    ## plotting heat maps for observation pairs -------

    if(type == 'heat_map') {

      if(reorder) {

        heat_map <- ggplot(plot_data,
                           aes(x = x_obs,
                               y = y_obs,
                               fill = distance))

      } else {

        heat_map <- ggplot(plot_data,
                           aes(x = reorder(x_obs, distance),
                               y = reorder(y_obs, distance),
                               fill = distance))

      }

      heat_map <- heat_map +
        ggplot2::geom_tile() +
        ggplot2::facet_grid(cluster2 ~ cluster1,
                            scales = 'free',
                            space = 'free') +
        ggplot2::scale_fill_gradient2(low = 'firebrick',
                                      mid = 'white',
                                      high = 'steelblue',
                                      name = paste0(attr(x, 'dist_method'),
                                                    '\ndistance'),
                                      midpoint = mean(range(plot_data$distance,
                                                            na.rm = TRUE))) +
        cust_theme +
        ggplot2::theme(axis.text = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       axis.line = ggplot2::element_blank()) +
        ggplot2::labs(title = 'Pairwise cross-distances',
                      subtitle = paste0(attr(x, 'type'),
                                        ' ',
                                        attr(x, 'dist_method'),
                                        ' pairwise cross-distances'),
                      x = 'x object observations',
                      y = if(attr(x, 'type') == 'homologous') 'x object observations' else 'y object observations')

      return(heat_map)

    }

    ## plotting heat maps for mean cross-distances -------

    if(type == 'mean') {

      heat_map <-
        ggplot(plot_data,
               aes(x = cluster1,
                   y = cluster2,
                   fill = mean)) +
        ggplot2::geom_tile(color = line_color) +
        ggplot2::geom_text(aes(label = plot_lab),
                           color = if(show_txt) 'black' else NA,
                           size = txt_size,
                           hjust = 0.5,
                           vjust = 0.5) +
        ggplot2::scale_fill_gradient2(low = 'firebrick',
                                      mid = 'white',
                                      high = 'steelblue',
                                      name = paste0(attr(x, 'dist_method'),
                                                    '\ndistance'),
                                      midpoint = mean(range(plot_data$mean,
                                                            na.rm = TRUE))) +
        cust_theme +
        ggplot2::theme(axis.line = ggplot2::element_blank()) +
        ggplot2::labs(title = 'Pairwise cross-distances',
                      subtitle = paste0(attr(x, 'type'),
                                        ' ',
                                        attr(x, 'dist_method'),
                                        ' pairwise cross-distances'),
                      x = 'x object clusters',
                      y = if(attr(x, 'type') == 'homologous') 'x object clusters' else 'y object clusters')

      return(heat_map)

    }

    ## plotting histogram panels -------

    if(type == 'histogram') {

      if(is.null(labeller)) {

        labeller <-
          ggplot2::labeller(.rows = set_names(paste('y:', attr(x, 'y_levels')),
                                              attr(x, 'y_levels')),
                            .cols = set_names(paste('x:', attr(x, 'x_levels')),
                                              attr(x, 'x_levels')))

      }

      hist_panel <-
        ggplot(plot_data,
               aes(x = distance,
                   fill = cluster2)) +
        ggplot2::geom_histogram(...) +
        ggplot2::facet_grid(cluster2 ~ cluster1,
                            labeller = labeller) +
        cust_theme +
        ggplot2::labs(title = 'Pairwise cross-distance distribution',
                      x = paste0(attr(x, 'type'),
                                 ' ',
                                 attr(x, 'dist_method'),
                                 ' pairwise cross-distances'),
                      y = 'observation count')

      return(hist_panel)

    }

  }

# END -----
