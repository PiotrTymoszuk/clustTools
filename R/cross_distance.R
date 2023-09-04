# Tools to compute cross-distances between data frames or clusters

# Cross-distance between data frames -----

#' Compute cross-distances between two data frames
#'
#' @description
#' Computes cross-distances between two data frames with the same variable sets.
#'
#' @details
#' Distances (for available distances, see: \code{\link{get_kernel_info}}) are
#' computed in a pair-wise manner employing \code{\link[philentropy]{distance}}.
#' Preserves row names.
#' If a single data frame is provided, pairwise observations
#' between the observations are computed with \code{\link{calculate_dist}}.
#' If a single `clust_analysis` or `combi_analysis` object is provided,
#' cross-distances between the clusters within the object are computed.
#' `cross_distance()` is a S3 generic function.
#'
#' @references
#' Drost H-G. Philentropy: Information Theory and Distance Quantification
#' with R. J Open Source Softw (2018) 3:765. doi:10.21105/joss.00765
#' @references
#' Sulc Z, Cibulkova J, Rezankova H. nomclust: Hierarchical Cluster
#' Analysis of Nominal Data. (2021)
#' Available at: https://cran.r-project.org/package=nomclust
#'
#' @param x a data frame, `clust_analysis` or `combi_analysis` object.
#' @param y an object like `x` or NULL (default).
#' @param method distance metric name as specified by
#' \code{\link{get_kernel_info}}.
#' For `clust_analysis` or `combi_analysis` instances.
#' If `method` is set to NULL,
#' the metric name is extracted from the object (distance between observations
#' for `clust_analysis` and `combi_analysis`, not between the SOM nodes).
#' @param ... extra arguments passed to methods.
#'
#' @return For data frames: a matrix with pairwise distances,
#' observations of the `x` data frame
#' are present in rows, observations of the `y` data frame are presented
#' in columns.
#' For `clust_analysis` and `combi_analysis` results:
#' a list of cross-distance matrices of class \code{\link{cross_dist}}
#' with defined \code{\link{summary.cross_dist}} and
#' \code{\link{plot.cross_dist}} methods.
#'
#' @export

  cross_distance <- function(x, ...) {

    UseMethod('cross_distance')

  }

#' @rdname cross_distance
#' @export cross_distance.data.frame
#' @export

  cross_distance.data.frame <- function(x,
                                        y = NULL,
                                        method = 'euclidean', ...) {

    ## entry control -------

    if(!is.data.frame(x)) {

      stop("'x' has to be a data frame.",
           call. = FALSE)

    }

    check_numeric(x)

    if(!is.null(y)) {

      if(!is.data.frame(y)) {

        stop("'y' has to be a data frame.",
             call. = FALSE)

      }

      check_numeric(y)

    }

    if(!method %in% get_kernel_info()) {

      stop(paste("Unsupported distance method. Please refer",
                 "to 'get_kernel_info()'"),
           call. = FALSE)

    }

    ## no y: pairwise distances -------

    if(is.null(y)) return(calculate_dist(x, method = method))

    ## common variables -------

    cmm_variables <- intersect(names(x), names(y))

    if(length(cmm_variables) == 0) {

      stop('No variables are shared between the inputa data frames.',
           call. = FALSE)

    }

    x <- x[cmm_variables]
    y <- y[cmm_variables]

    ## a common computation data frame and distance matrix ---------

    x <- as.data.frame(x)
    y <- as.data.frame(y)

    x <- set_rownames(x, row_names = paste0('x_', rownames(x)))
    y <- set_rownames(y, row_names = paste0('y_', rownames(y)))

    cmm_data <- rbind(x, y)

    cmm_dist <- calculate_dist(cmm_data, method = method)

    cmm_dist[rownames(x), rownames(y)]

  }

#' @rdname cross_distance
#' @export cross_distance.clust_analysis
#' @export

  cross_distance.clust_analysis <- function(x,
                                            y = NULL,
                                            method = NULL, ...) {

    ## entry control -------

    if(!is_clust_analysis(x) & !is_combi_analysis(x)) {

      stop(paste("'x' has to be an instance of the 'clust_analysis'",
                 "or 'combi_analysis' class."),
           call. = FALSE)

    }

    if(!is.null(y)) {

      if(!is_clust_analysis(y) & !is_combi_analysis((y))) {

        stop(paste("'y' has to be an instance of the 'clust_analysis'",
                   "or 'combi_analysis' class."),
             call. = FALSE)

      }

    }

    if(is.null(method)) {

      if(is_clust_analysis(x)) {

        method <- x$dist_method

      } else {

        method <- x$clust_analyses$observation$dist_method

      }

    }

    if(!method %in% get_kernel_info()) {

      stop(paste("Unsupported distance metric.",
                 "Please refer to 'get_kernel_info()"),
           call. = FALSE)

    }

    ## no y provided: cross-distances between the clusters -------

    if(is.null(y)) {

      clust_ass <- extract(x, 'assignment')

      clust_ass <- split(clust_ass$observation, clust_ass$clust_id)

      if(is_clust_analysis(x)) {

        clust_data <- model.frame(x)

      } else {

        clust_data <- model.frame(x)$observation

      }

      clust_data <-
        map(clust_ass,
            ~clust_data[.x, ])

      pairs <- utils::combn(names(clust_data), m = 2, simplify = FALSE)

      identities <- map(names(clust_data), ~c(.x, .x))

      pairs <- c(pairs, identities)

      pair_names <- map_chr(pairs, paste, collapse = ' vs ')

      pairs <- set_names(pairs, pair_names)

      dist_lst <-
        map(pairs,
            ~cross_distance(clust_data[[.x[1]]],
                            clust_data[[.x[2]]],
                            method = method))

      return(cross_dist(dist_lst,
                        type = 'homologous',
                        method = method,
                        x_levels = names(clust_ass),
                        y_levels = names(clust_ass)))

    }

    ## y provided: cross-distances between the clusters of the objects -------

    clust_list <- list(x = x,
                       y = y)

    ### cluster assignment and cluster pairs

    clust_ass <- map(clust_list, extract, 'assignment')

    clust_ass <- map(clust_ass,
                     ~split(.x$observation, .x$clust_id))

    pairs <-
      map(names(clust_ass$x),
          function(x) map(names(clust_ass$y),
                          ~c(x, .x)))

    pairs <- unlist(pairs, recursive = FALSE)

    pair_names <- map(pairs, paste, collapse = ' vs ')

    pairs <- set_names(pairs, pair_names)

    ### clustering data, reduced to common variables

    clust_combi <- map(clust_list, is_combi_analysis)

    clust_data <- map(clust_list, model.frame)

    clust_data <-
      map2(clust_data, clust_combi,
           function(x, y) if(y) x$observation else x)

    cmm_vars <- map(clust_data, names)

    cmm_vars <- reduce(cmm_vars, intersect)

    if(length(cmm_vars) == 0) {

      stop("No common variables for 'x' and 'y'.", call. = FALSE)

    }

    clust_data <- map(clust_data, ~.x[cmm_vars])

    clust_data <- map2(clust_data, clust_ass,
                       function(data, clust) map(clust, ~data[.x, ]))


    dist_lst <- map(pairs,
                    ~cross_distance(clust_data$x[[.x[1]]],
                                    clust_data$y[[.x[2]]],
                                    method = method))

    return(cross_dist(dist_lst,
                      type = 'heterologous',
                      method = method,
                      x_levels = names(clust_ass[['x']]),
                      y_levels = names(clust_ass[['y']])))

  }

#' @rdname cross_distance
#' @export cross_distance.combi_analysis
#' @export

  cross_distance.combi_analysis <- function(x,
                                            y = NULL,
                                            method = NULL, ...) {

    cross_distance.clust_analysis(x = x, y = y, method = method)

  }

# END ------
