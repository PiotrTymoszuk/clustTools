# Tools to compute cross-distances between data frames or clusters

# Cross-distance between data frames -----

#' Compute cross-distances between two data frames
#'
#' @description
#' Computes cross-distances between two data frames with the same variable sets.
#'
#' @details
#' `cross_distance()` is a S3 generic function.
#' Distances (for available distances, see: \code{\link{get_kernel_info}}) are
#' computed in a pair-wise manner employing \code{\link[philentropy]{distance}}.
#' Preserves row names.
#' If a single data frame is provided, pairwise observations
#' between the observations are computed with \code{\link{calculate_dist}}.
#' If a single `clust_analysis` or `combi_analysis` object is provided,
#' cross-distances between the clusters within the object are computed -
#' so called 'homologous' distances, as opposed to 'heterologous' distances
#' computed in a pair-wise manner between clusters of two clustering analysis
#' objects.
#' Note: it is not possible to compute heterologous distances is cases of
#' clustering analyses done with an user-provided dissimilarity matrix
#' (subclass `min_analysis` of `clust_analysis` parent class). In such cases,
#' the `method` argument is ignored as well.
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
#' for `clust_analysis` and `combi_analysis`, not between the SOM nodes). For
#' multi-layer SOM and their prediction, the distance methods are extracted
#' from the clustering objects.
#' @param .parallel logical, should the operation be run in parallel?
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

    if(!is.matrix(cmm_dist)) {

      warning('Distance calculation failed: NAs or not enough data points?',
              call. = FALSE)

      return(NA)

    }

    cmm_dist[rownames(x), rownames(y), drop = FALSE]

  }

#' @rdname cross_distance
#' @export cross_distance.clust_analysis
#' @export

  cross_distance.clust_analysis <- function(x,
                                            y = NULL,
                                            method = NULL,
                                            .parallel = FALSE, ...) {

    if(.parallel) future::plan('multisession')

    on.exit(future::plan('sequential'))

    multi_types <- c('supersom', 'supersom_prediction')

    if(is.null(y)) {

      if(!is_multi_layer(x)) {

        return(cross_single_homolog(x, method))

      } else {

        return(cross_multi_homolog(x))

      }

    }

    if(!is.null(y)) {

      if(!is_multi_layer(x) & !is_multi_layer(y)) {

        return(cross_single_heterolog(x, y, method))

      }

      warn_txt <- "Incompatible analysis types: SOM and non-SOM clustering"

      if(is_multi_layer(x) & !is_multi_layer(y)) {

        warning(warn_txt, call. = FALSE)

        return(NULL)

      }

      if(!is_multi_layer(x) & is_multi_layer(y)) {

        warning(warn_txt, call. = FALSE)

        return(NULL)

      }

      return(cross_multi_heterolog(x, y))

    }

  }

#' @rdname cross_distance
#' @export cross_distance.min_analysis
#' @export

  cross_distance.min_analysis <- function(x, ...) {

    ## the function extracts pair-wise distances within and
    ## between the clusters from the general distance matrix provided
    ## within the analysis object

    ## clust assignment and distances -------

    general_dist <- as.matrix(dist(x))

    clust_assignment <- extract(x, 'assignment')

    clust_assignment <-
      split(clust_assignment$observation, clust_assignment$clust_id)

    ## cluster pairs --------

    pairs <- utils::combn(names(clust_assignment), m = 2, simplify = FALSE)

    identities <- map(names(clust_assignment), ~c(.x, .x))

    pairs <- c(pairs, identities)

    pair_names <- map_chr(pairs, paste, collapse = ' vs ')

    pairs <- set_names(pairs, pair_names)

    pairs <-
      map(pairs,
          ~list(clust_assignment[[.x[1]]],
                clust_assignment[[.x[2]]]))

    ## cross-distance objects

    dist_lst <-
      map(pairs,
          ~general_dist[.x[[1]], .x[[2]]], drop = FALSE)

    cross_dist(dist_lst,
               type = 'homologous',
               method = 'custom',
               x_levels = names(clust_assignment),
               y_levels = names(clust_assignment))

  }

#' @rdname cross_distance
#' @export cross_distance.combi_analysis
#' @export

  cross_distance.combi_analysis <- function(x,
                                            y = NULL,
                                            method = NULL,
                                            .parallel = FALSE, ...) {

    if(.parallel) future::plan('multisession')

    on.exit(future::plan('sequential'))

    if(is.null(y)) {

      return(cross_single_homolog(x, method))

    } else {

      return(cross_single_heterolog(x, y, method))

    }

  }

#' @rdname cross_distance
#' @export

  cross_distance.umatrix_analysis <- function(x,
                                              y = NULL,
                                              method = NULL,
                                              .parallel = FALSE, ...) {

    if(.parallel) future::plan('multisession')

    on.exit(future::plan('sequential'))

    if(is.null(y)) {

      return(cross_multi_homolog(x))

    } else {

      if(!is_multi_layer(y)) {

        warning(paste("'y' has to be an analysis or prediction objects",
                      "employing a multi-layer SOM."),
                call. = FALSE)

        return(NULL)

      }

      return(cross_multi_heterolog(x, y))

    }

  }
# END ------
