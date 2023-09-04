# Utilities for plotting and QC

# Distance calculation -----

#' Calculate distance.
#'
#' @description
#' The `calculate_dist()` function calculates distances
#' between observations (rows) of a data
#' frame or matrix. See \code{\link{get_kernel_info}} for a vector of available
#' distance measures.
#'
#' @details
#' Provides a handy wrapper for \code{\link[philentropy]{distance}}.
#' The smc (simple matching coefficient) distance is calculated with the
#' \code{\link[nomclust]{sm}} function.
#'
#' @references
#' Drost H-G. Philentropy: Information Theory and Distance Quantification
#' with R. J Open Source Softw (2018) 3:765. doi:10.21105/joss.00765
#' @references
#' Boriah S, Chandola V, Kumar V. Similarity measures for categorical data:
#' A comparative evaluation. in Society for Industrial and
#' Applied Mathematics - 8th SIAM International Conference on Data Mining
#' 2008, Proceedings in Applied Mathematics 130, 243–254.
#' doi:10.1137/1.9781611972788.22
#' @references
#' Sulc Z, Cibulkova J, Rezankova H. nomclust: Hierarchical Cluster Analysis
#' of Nominal Data. (2021)
#' Available at: https://cran.r-project.org/package=nomclust
#'
#' @return `calculate_dist()`: a matrix with the distance statistics;
#' `get_kernel_info()`: a vector with names of available distance measures.
#'
#' @param data a data frame or matrix. Row names are preserved as
#' observation IDs.
#' @param method the name of dissimilarity measure, see:
#' \code{\link{get_kernel_info}} for available distances.
#'
#' @export

  calculate_dist <- function(data, method) {

    if(!is.data.frame(data) & !is.matrix(data)) {

      stop('Provide a valid data.frame, tibble or matrix as input data',
           call. = FALSE)

    }

    if(!method %in% get_kernel_info()) {

      stop('The requested distance measure is not available.',
           call. = FALSE)

    }

    av_distances <- get_kernel_info()

    if(method == 'smc') {

      data <- as.data.frame(data)

      dist_mtx <- as.matrix(nomclust::sm(data))

    } else {

      if(method == 'sumofsquares') method <- 'squared_euclidean'

      data <- as.matrix(data)

      dist_mtx <- philentropy::distance(data,
                                        method = method,
                                        use.row.names = TRUE,
                                        mute.message = TRUE)

      if(method %in% c('cosine', 'ruzicka')) {

        ### handling the similarity coefficients

        dist_mtx <- 1 - dist_mtx

      }

    }

    dist_mtx

  }

#' @rdname calculate_dist
#' @export

  get_kernel_info <- function() {

    av_distances <- somKernels::get_kernel_info()

    av_distances[av_distances != 'BrayCurtis']

  }

# Specific plotting helpers -----

#' Plot the mean distance to k-nearest neighbors.
#'
#' @description
#' Plots the sorted (ascending) distances to k-nearest neighbors
#' (kNN) for each observation in the provided dissimilarity object.
#'
#' @details
#' Internally, the mean kNN distances are calculated with the
#' \code{\link[dbscan]{kNNdist}} function.
#'
#' @references
#' Hahsler M, Piekenbrock M, Doran D. Dbscan: Fast density-based clustering
#' with R. J Stat Softw (2019) 91:1–30. doi:10.18637/jss.v091.i01
#' @references
#' Belyadi H, Haghighat A, Nguyen H, Guerin A-J. IOP Conference Series:
#' Earth and Environmental Science Determination of Optimal Epsilon (Eps)
#' Value on DBSCAN Algorithm to Clustering Data on Peatland Hotspots in
#' Sumatra Related content EPS conference comes to London-EPS rewards
#' quasiparticle research-EP. IOP Conf Ser Earth Environ Sci (2016) 31:
#' doi:10.1088/1755-1315/31/1/012012
#'
#' @return A ggplot object.
#'
#' @param diss_obj a dissimilarity object (e.g. 'dist' class).
#' @param k the k number of the nearest neighbors.
#' @param eps the distance to be presented in the plot as a horizontal dashed
#' line. If NULL, the line is hidden.
#' @param plot_title plot title.
#' @param plot_subtitle plot subtitle.
#' @param plot_tag plot tag.
#' @param cust_theme custom plot theme, a ggplot2 theme object.
#'
#' @export

  plot_knn_distance <- function(diss_obj,
                                k,
                                eps = NULL,
                                plot_title = NULL,
                                plot_subtitle = NULL,
                                plot_tag = NULL,
                                cust_theme = ggplot2::theme_classic()) {

    stopifnot(inherits(diss_obj, 'dist'))
    stopifnot(inherits(cust_theme, 'theme'))

    point <- NULL
    knn_dist <- NULL

    sort_distances <- dbscan::kNNdist(x = diss_obj,
                                      k = k,
                                      all = FALSE)

    sort_distances <- sort(sort_distances)

    sort_distances <-
      tibble(sample = names(sort_distances),
             knn_dist = unname(sort_distances),
             point = 1:length(sort_distances))

    dist_plot <- ggplot(data = sort_distances,
                        aes(x = point,
                            y = knn_dist)) +
      ggplot2::geom_area(color = 'cornsilk4',
                         fill = 'cornsilk2') +
      cust_theme +
      ggplot2::labs(x = 'Sample',
                    y = paste(k, 'NN distance', sep = '-'),
                    title = plot_title,
                    subtitle = plot_subtitle,
                    tag = plot_tag)

    if(!is.null(eps)) {

      dist_plot <- dist_plot +
        ggplot2::geom_hline(yintercept = eps,
                            linetype = 'dashed',
                            color = 'black')

    }

    dist_plot

  }

#' Plot a dendrogram.
#'
#' @description
#' Plots a dendrogram given a clustering object generated by
#' \code{\link[stats]{hclust}}.
#'
#' @details
#' The dendrogram structure is generated with the
#' \code{\link[stats]{as.dendrogram}} function and graphical layout provided
#' by \code{\link[dendextend]{color_branches}} and
#' \code{\link[dendextend]{set}}.
#'
#' @references
#' Galili T. dendextend: an R package for visualizing, adjusting and
#' comparing trees of hierarchical clustering. Bioinformatics (2015) 31:3718–20.
#' doi:10.1093/bioinformatics/btv428
#'
#' @return A ggplot object.
#'
#' @param clust_str an object of the 'hclust' class.
#' @param k an integer, the cluster number.
#' @param labels logical, should observation labels be presented in the x axis?
#' @param cluster_colors colors of the cluster branches, a vector of the length
#' k + 1. The last color codes for the connecor branches.
#' @param cluster_labels cluster names, a text vector of the lenght k.
#' @param cluster_leg_title cluster legend title.
#' @param plot_title plot title.
#' @param plot_subtitle plot subtitle.
#' @param plot_tag plot tag.
#' @param y_lab y axis title.
#' @param cust_theme custom plot theme, a ggplot2 theme object.
#' @param ... extra arguments, currently none.
#'
#' @export

  plot_dendro <- function(clust_str,
                          k,
                          labels = TRUE,
                          cluster_colors = NULL,
                          cluster_labels = paste0('Cluster #', 1:k),
                          cluster_leg_title = 'Cluster',
                          plot_title = NULL,
                          plot_subtitle = NULL,
                          plot_tag = NULL,
                          y_lab = NULL,
                          cust_theme = ggplot2:: theme_classic(), ...) {

    ## entry control

    stopifnot(inherits(clust_str, 'hclust'))
    stopifnot(is.logical(labels))
    stopifnot(inherits(cust_theme, 'theme'))

    k <- as.integer(k)

    ## dendrogram data

    dendro_data <- stats::as.dendrogram(clust_str)

    dendro_data <- dendextend::color_branches(dend = dendro_data, k = k)

    dendro_data <- dendextend::set(dend = dendro_data,
                                   'branches_lwd', 0.5)

    dendro_data <- dendextend::set(dend = dendro_data,
                                   'labels_cex', 0.5)

    ## plotting

    dendro_plot <- ggplot2::ggplot(data = dendro_data,
                                   labels = labels) +
      cust_theme +
      ggplot2::theme(axis.line.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank()) +
      ggplot2::labs(title = plot_title,
                    subtitle = plot_subtitle,
                    tag = plot_tag,
                    y = y_lab)

    if(all(c(!is.null(cluster_colors),
             !is.null(cluster_labels),
             !is.null(cluster_leg_title)))) {

      suppressMessages(dendro_plot <- dendro_plot +
                         ggplot2::scale_color_manual(values = cluster_colors,
                                                     labels = cluster_labels,
                                                     name = cluster_leg_title) +
                         ggplot2::guides(color = ggplot2::guide_legend()))

    }

    dendro_plot

  }

#' Plot WSS curve and silhouette statistic values as a function of cluster
#' number.
#'
#' @description
#' Plots the values of the total within-cluster sum-of-squares and
#' silhouette statistic as a function of the cluster number.
#'
#' @details
#' Takes a distance matrix (e.g. the \code{\link{get_kernel_info}}
#' output) and a clustering function, for the details, see:
#' \code{\link[factoextra]{fviz_nbclust}}.
#'
#' @return a ggplot object.
#'
#' @references
#' Kassambara A, Mundt F. factoextra: Extract and Visualize the Results
#' of Multivariate Data Analyses. (2020) Available at:
#' https://cran.r-project.org/web/packages/factoextra/index.html
#' @references
#' Rousseeuw PJ. Silhouettes: A graphical aid to the interpretation and
#' validation of cluster analysis. J Comput Appl Math (1987) 20:53–65.
#' doi:10.1016/0377-0427(87)90125-7
#'
#' @param data a numeric matrix with the distances or a data frame.
#' @param k an integer, the cluster number.
#' @param FUNcluster a clustering function. See:
#' \code{\link[factoextra]{fviz_nbclust}} for details.
#' @param method a statistic to be plotted. See:
#' \code{\link[factoextra]{fviz_nbclust}} for details.
#' @param plot_title plot title.
#' @param plot_subtitle plot subtitle.
#' @param plot_tag plot tag.
#' @param cust_theme custom plot theme, a ggplot2 theme object.
#' @param ... extra arguments passed to \code{\link[factoextra]{fviz_nbclust}}.
#'
#' @export

  plot_nbclust <- function(data,
                           k,
                           FUNcluster = NULL,
                           method = c('wss', 'silhouette', 'gap_stat'),
                           plot_title = NULL,
                           plot_subtitle = NULL,
                           plot_tag = NULL,
                           cust_theme = ggplot2::theme_classic(), ...) {

    ## entry control

    check_numeric(data)

    stopifnot(inherits(cust_theme, 'theme'))

    k <- as.integer(k)

    method <- match.arg(method[1],
                        c('silhouette', 'wss', 'gap_stat'))

    ## plotting

    factoextra::fviz_nbclust(data,
                             FUNcluster = FUNcluster,
                             method = method, ...) +
      ggplot2::geom_vline(xintercept = k,
                          linetype = 'dashed',
                          color = 'coral3') +
      cust_theme +
      ggplot2::labs(title = plot_title,
                    subtitle = plot_subtitle,
                    tag = plot_tag)

  }

#' Plot diagnostic plots for the self-organizing map.
#'
#' @description
#' Generates a set of diagnostic plots for the 'kohonen' class
#' object as specified by \code{\link[kohonen]{plot.kohonen}}.
#'
#' @details The plots are wrapped as a list of (non-editable) ggplot objects.
#'
#' @references
#' Wehrens R, Kruisselbrink J. Flexible self-organizing maps in kohonen 3.0.
#' J Stat Softw (2018) 87:1–18. doi:10.18637/jss.v087.i07
#'
#' @return  a list of non-editable ggplot objects.
#'
#' @param kohonen_object a 'kohonen' class object.
#' See: \code{\link[kohonen]{som}} for details.
#'
#' @export

  plot_som <- function(kohonen_object) {

    ## generates diagnostic plots for the test cohort

    stopifnot(inherits(kohonen_object, 'kohonen'))

    plot_exprs <- rlang::enexpr(kohonen_object)

    plot_methods <- c('change',
                      'codes',
                      'counts',
                      'mapping',
                      'dist.neighbours',
                      'quality')

    plot_exprs_lst <- purrr::map(plot_methods,
                                 function(x) rlang::expr(~plot(!!plot_exprs,
                                                               !!x)))

    plot_exprs_lst <- rlang::set_names(plot_exprs_lst,
                                       plot_methods)

    purrr::map(plot_exprs_lst,
               function(x)
                 cowplot::ggdraw(cowplot::as_grob(rlang::eval_tidy(x))))

  }

#' Visualize the SOM training process.
#'
#' @description
#' Plots the mean distance to the neuron/winning unit as a
#' function of the iteration number.
#'
#' @references
#' Wehrens R, Kruisselbrink J. Flexible self-organizing maps in kohonen 3.0.
#' J Stat Softw (2018) 87:1–18. doi:10.18637/jss.v087.i07
#'
#' @return a ggplot object.
#'
#' @param kohonen_object a 'kohonen' class object.
#' @param plot_title plot title.
#' @param plot_subtitle plot subtitle.
#' @param cust_theme custom plot theme, a ggplot2 theme object.
#' @param ... extra arguments, currently none specified.
#'
#' @export

  plot_train_som <- function(kohonen_object,
                             plot_title = NULL,
                             plot_subtitle = NULL,
                             cust_theme = ggplot2::theme_classic(), ...) {

    ## plots the distance to the winning unit (neuron) during the training process
    ## a nicer plot using ggplot(). May pass additional arguments to the smoothing function

    stopifnot(inherits(kohonen_object, 'kohonen'))
    stopifnot(inherits(cust_theme, 'theme'))

    Iteration <- NULL

    change_data <- as.data.frame(kohonen_object$changes)

    change_data <- rlang::set_names(change_data, 'dist')

    change_data <- tibble::rownames_to_column(change_data, 'Iteration')

    ggplot2::ggplot(change_data,
                    ggplot2::aes(x = as.numeric(Iteration),
                                 y = dist)) +
      ggplot2::geom_point(shape = 16,
                          size = 1,
                          color = 'gray60') +
      ggplot2::geom_smooth(color = 'steelblue',
                           method = 'loess', ...) +
      cust_theme +
      ggplot2::labs(x = 'Iteration',
                    y = 'Mean distance\nto the winning unit',
                    title = plot_title,
                    subtitle = plot_subtitle,
                    tag = paste0('\nIterations: n = ',
                                 nrow(change_data)))

  }

# General plotting functions -------

#' Generate a custom scatter ggplot.
#'
#' @description
#' Generates a simple scatter ggplot.
#'
#' @return a ggplot object.
#'
#' @param data a data frame.
#' @param x_var the name of the variable to be presented in the x axis.
#' @param y_var the name of the variable to be presented in the y axis.
#' @param fill_var optional, the name of the variable coded by the point fill.
#' If NULL, the point fill is specified by the point_color argument.
#' @param label_var optional, the name of the variable to be presented in the
#' point labels. If NULL, no point labels are displayed.
#' @param plot_title plot title.
#' @param plot_subtitle plot subtitle.
#' @param plot_tag plot tag.
#' @param x_lab x axis title.
#' @param y_lab y axis title.
#' @param fill_lab fill legend title.
#' @param cust_theme custom plot theme, a ggplot2 theme object.
#' @param point_color point fill color.
#' @param point_alpha point alpha.
#' @param show_segments logical, should lines connecting the (0,0) point with
#' the plot point be displayed?
#' @param segment_color color of the connecting lines.
#' @param segment_alpha alpha of the connecting lines.
#' @param label_color color of the text labels.
#' @param txt_color color of the text presented in the labels.
#' @param txt_size size of the text presented in the labels.
#' @param txt_type type of the displayed text: either as geom_text or geom_label
#' @param jitter_width horizontal jittering of the points.
#' @param jitter_height vertical jittering of the points.
#'
#' @export

  plot_point <- function(data,
                         x_var,
                         y_var,
                         fill_var = NULL,
                         label_var = NULL,
                         plot_title = NULL,
                         plot_subtitle = NULL,
                         plot_tag = NULL,
                         x_lab = x_var,
                         y_lab = y_var,
                         fill_lab = NULL,
                         cust_theme = ggplot2::theme_classic(),
                         point_color = 'steelblue',
                         point_alpha = 1,
                         show_segments = FALSE,
                         segment_color = 'steelblue',
                         segment_alpha = 1,
                         label_color = point_color,
                         txt_color = 'black',
                         txt_size = 2.5,
                         txt_type = c('label', 'text'),
                         jitter_width = 0,
                         jitter_height = 0) {

    ## entry control

    stopifnot(is.data.frame(data))
    stopifnot(is.logical(show_segments))

    if(!x_var %in% names(data) | !y_var %in% names(data)) {

      stop('x_var or y_var absent from the data frame.', call. = FALSE)

    }

    if(!is.null(fill_var)) {

      if(!fill_var %in% names(data)) {

        stop('fill_var absent from the data frame.', call. = TRUE)

      }

    }

    if(!is.null(label_var)) {

      if(!label_var %in% names(data)) {

        stop('label_var absent from the data frame.', call. = TRUE)

      }

    }

    ## plotting

    if(is.null(fill_var)) {

      pplot <- ggplot(data = data,
                      aes(x = .data[[x_var]],
                          y = .data[[y_var]]))

    } else {

      pplot <- ggplot(data = data,
                      aes(x = .data[[x_var]],
                          y = .data[[y_var]],
                          fill = .data[[fill_var]]))

    }

    if(show_segments) {

      pplot <- pplot +
        ggplot2::geom_segment(ggplot2::aes(x = 0,
                                           y = 0,
                                           xend = .data[[x_var]],
                                           yend = .data[[y_var]]),
                              color = segment_color,
                              alpha = segment_alpha)

    }

    if(is.null(fill_var)) {

      pplot <- pplot +
        ggplot2::geom_point(shape = 21,
                            size = 2,
                            fill = point_color,
                            alpha = point_alpha,
                            position = ggplot2::position_jitter(width = jitter_width,
                                                                height = jitter_height))

    } else {

      pplot <- pplot +
        ggplot2::geom_point(shape = 21,
                            size = 2,
                            alpha = point_alpha,
                            position = ggplot2::position_jitter(width = jitter_width,
                                                                height = jitter_height))

    }

    if(!is.null(label_var)) {

      txt_type <- match.arg(txt_type[1], c('label', 'text'))

      if(txt_type == 'label') {

        pplot <- pplot +
          ggrepel::geom_label_repel(aes(label = .data[[label_var]]),
                                    size = txt_size,
                                    fill = label_color,
                                    color = txt_color,
                                    box.padding = 0.1,
                                    label.padding = 0.1)

      } else {

        pplot <- pplot +
          ggrepel::geom_text_repel(aes(label = .data[[label_var]]),
                                   size = txt_size,
                                   color = txt_color,
                                   box.padding = 0.1)

      }

    }

    pplot +
      cust_theme +
      ggplot2::labs(title = plot_title,
                    subtitle = plot_subtitle,
                    tag = plot_tag,
                    x = x_lab,
                    y = y_lab,
                    fill = fill_lab)


  }

# Calculation helpers ----

#' Calculate clustering sum of squares.
#'
#' @description
#' Calculates total, within cluster and between cluster
#' sum of squares (ss).
#'
#' @details
#' The calculation method is independent of the clustering method.
#'
#' @return a list with the values of within-cluster ss for the particular
#' clusters, total within-cluster ss, total ss, total between-cluster ss as well
#' as the ratio of between-cluster ss to total ss, interpreted as the fraction
#' of 'explained' clustering variance.
#'
#' @param dist_mtx a numeric matrix with the distances.
#'
#' @param assignment a data frame with the variable 'clust_id' specifying the
#' assignment of the observations to the clusters.

  get_sum_sq <- function(dist_mtx, assignment) {

    ## entry control

    stopifnot(is.matrix(dist_mtx))
    stopifnot(is.data.frame(assignment))

    ## within sum of squares

    x_ss <- stats::aggregate(dist_mtx,
                             by = assignment[, 'clust_id'],
                             FUN = function(x) sum(scale(x, scale = FALSE)^2,
                                                   na.rm = TRUE))

    wss <- rowSums(x_ss[, -1])

    total_wss <- sum(x_ss[, -1])

    ## total sum of squares

    total_ss <- sum(scale(dist_mtx, scale = FALSE)^2, na.rm = TRUE)

    ## output

    list(wss = wss,
         total_wss = total_wss,
         total_ss = total_ss,
         between_ss = total_ss - total_wss,
         frac_var = (total_ss - total_wss)/total_ss)

  }

#' Get dimensions of a data frame or matrix.
#'
#' @description
#' Gets the number of observations and variables of
#' the given object.
#'
#' @return a list with the requested statistics.
#'
#' @param data a data frame or matrix.

  get_data_dim <- function(data) {

    stopifnot(is.data.frame(data) | is.matrix(data))

    list(observations = nrow(data),
         variables = ncol(data))

  }

#' Find the most frequently occurring element of a vector.
#'
#' @description
#' Finds the element of a vector with the highers number of
#' occurrences.
#'
#' @details Ties may be resolved at random (resolve_ties = TRUE), otherwise,
#' if a tie exists, the alphabetically first element is returned.
#'
#' @return the most frequent element.
#'
#' @param vector a vector.
#' @param resolve_ties logical, should the ties be resolved at random?

  vote_simple <- function(vector, resolve_ties = FALSE) {

    stopifnot(is.vector(vector))
    stopifnot(is.logical(resolve_ties))

    voting_res <- sort(table(vector), decreasing = TRUE)

    best_res <- names(voting_res[voting_res == max(voting_res)])

    if(!resolve_ties) return(best_res[1])

    if(length(best_res) == 1) {

      return(best_res[1])

    } else {

      best_index <- sample(1:length(best_res), size = 1)

      return(best_res[best_index])

    }

  }

#' Find the most frequently occurring element with distance weighting.
#'
#' @description
#' Finds the element of a vector with the highers number of
#' occurrences. The voting is distance weighted by the given kernel function.
#'
#' @references the most frequent element.
#'
#' @param vector a vector.
#' @param dist_vec a numeric vector with the distance values.
#' @param kernel_fun a kernel function.

  vote_kernel <- function(vector,
                          dist_vec,
                          kernel_fun = function(x) 1/x) {

    ## entry control -------

    stopifnot(is.vector(vector))
    stopifnot(is.vector(dist_vec))
    stopifnot(is.numeric(dist_vec))
    stopifnot(is.function(kernel_fun))

    raw_votes <- NULL
    weighted_votes <- NULL
    vote_sum <- NULL

    ## voting

    vote_tbl <- data.frame(raw_votes = vector,
                           weighted_votes = kernel_fun(dist_vec))

    vote_tbl <- filter(vote_tbl,
                       complete.cases(vote_tbl))

    vote_sums <- dplyr::group_by(vote_tbl, raw_votes)

    vote_sums <- dplyr::summarise(vote_sums,
                                  vote_sum = sum(weighted_votes))

    vote_sums <- dplyr::arrange(vote_sums, -vote_sum)

    vote_sums[['raw_votes']][1]

  }

# Entry control -----

#' Check for a numeric data frame or a matrix.
#'
#' @description
#' Checks if an object is a numeric data frame or a matrix.
#'
#' @param object an object.
#'
#' @return none. Throws exceptions if the object is not a numeric data frame or
#' a matrix.

  check_numeric <- function(object) {

    err_msg <- 'A numeric data frame is required.'

    if(is.data.frame(object)) {

      classes <- map_lgl(object, is.numeric)

      if(any(!classes)) {

        stop(err_msg, call. = FALSE)

      }

    } else if(is.matrix(object)){

      if(!is.numeric(object)) {

        stop(err_msg, call. = FALSE)

      }

    } else {

      stop(err_msg, call. = FALSE)

    }

  }

# END -----
