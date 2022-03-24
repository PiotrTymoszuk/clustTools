# S3 OOP methods for the importance class.

# Plotting ------

#' Plot feature importance as a scatter or bar plot.
#'
#' @description Generates a bar plot with the importance statistic for
#' the clustering variables. The importance statistic is the difference in the
#' clustering variance fraction between the original clustering structure and
#' the clustering objects with the given variable reshuffled randomly.
#' @return a ggplot bar or scatter plot.
#' @param x an importance object.
#' @param type type of the plot: scatter or bar. Defaults to scatter.
#' @param fill_color fill color for the bars or points.
#' @param label logical, should the points be labeled with the importance
#' stat value?
#' @param txt_size label text size.
#' @param signif_digits significant digits for rounding of the statistic value.
#' @param plot_title plot title.
#' @param plot_subtitle plot subtitle.
#' @param plot_tag plot tag.
#' @param cust_theme a ggplot theme.
#' @export plot.importance
#' @export

  plot.importance <- function(x,
                              type = c('scatter', 'bar'),
                              fill_color = 'cornsilk3',
                              label = TRUE,
                              txt_size = 2.75,
                              signif_digits = 2,
                              plot_title = NULL,
                              plot_subtitle = NULL,
                              plot_tag = NULL,
                              cust_theme = ggplot2::theme_classic()) {

    ## entry control

    stopifnot(any(class(x) == 'importance'))

    type <- match.arg(type[1], c('scatter', 'bar'))

    stopifnot(is.logical(label))
    stopifnot(any(class(cust_theme) == 'theme'))

    ## plotting

    base_plot <- ggplot2::ggplot(dplyr::filter(x, variable != 'data'),
                                 ggplot2::aes(x = frac_diff,
                                              y = reorder(variable, frac_diff))) +
      cust_theme +
      ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
      ggplot2::labs(title = plot_title,
                    subtitle = plot_subtitle,
                    tag = plot_tag,
                    x = expression(Delta*' clustering variance'))

    if(type == 'scatter') {

      base_plot <- base_plot +
        ggplot2::geom_point(size = 2,
                            shape = 21,
                            fill = fill_color)

    } else {

      base_plot <- base_plot +
        ggplot2::geom_bar(stat = 'identity',
                          color = 'black',
                          fill = fill_color)

    }

    if(label) {

      base_plot <- base_plot +
        ggplot2::geom_text(ggplot2::aes(label = signif(frac_diff,
                                                       signif_digits)),
                           size = txt_size,
                           vjust = if(type == 'bar') 0.5 else -0.9,
                           hjust = if(type == 'bar') 1.2 else 0.6)

    }

    base_plot

  }

# END -----
