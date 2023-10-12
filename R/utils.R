# Utilities for plotting and QC

# Distance calculation -----

#' Calculate distance.
#'
#' @description
#' The `calculate_dist()` function calculates distances
#' between observations (rows) of a data
#' frame or matrix.
#' `calculate_weighted_dist()` computes a matrix of weighted distances
#' between observations (rows) for a list of numeric data frames or matrices.
#' See \code{\link{get_kernel_info}} for a vector of available
#' distance measures.
#'
#' @details
#' `calculate_dist()` and `calculate_weighted_dist()` provide handy wrappers for
#' \code{\link[philentropy]{distance}}.
#' The smc (simple matching coefficient) distance is calculated with the
#' \code{\link[nomclust]{sm}} function. Similarity coefficients returned by
#' \code{\link[philentropy]{distance}} (methods: cosine, ruzicka, intersection,
#' inner_product, harmonic_mean, hassebrook, fidelity) are handled with the
#' formula `dist = 1 - simil`.
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
#' @param data a numeric data frame or matrix (`calculate_dist()`) or a list
#' of such objects (for `calculate_weighted_dist()`). Row names are preserved as
#' observation IDs.
#' @param method the name of dissimilarity measure (for `calculate_dist()`)
#' or a vector of distance nammes (`calculate_weighted_dist()`). See:
#' \code{\link{get_kernel_info}} for available distances.
#' @param weights a numeric vector of weights.
#' @param FUN a function used to integrate the weighted distances.
#' @param ... extra arguments passed to `FUN`.
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

      simil <- c('cosine', 'ruzicka', 'intersection',
                 'inner_product', 'harmonic_mean', 'hassebrook',
                 'fidelity')

      if(method %in% simil) {

        ### handling the similarity coefficients

        dist_mtx <- 1 - dist_mtx

      }

    }

    dist_mtx

  }

#' @rdname calculate_dist
#' @export

  calculate_weighted_dist <- function(data,
                                      method = 'euclidean',
                                      weights = 1,
                                      FUN = function(x, y) x + y, ...) {

    ## entry control --------

    err_txt <- "'data' has to be a list of numeric data frames or matrices."

    if(!is.list(data)) {

      stop(err_txt, call. = FALSE)

    }

    class_check <- map_lgl(data, function(x) is.data.frame(x) | is.matrix(x))

    if(any(!class_check)) {

      stop(err_txt, call. = FALSE)

    }

    if(length(method) != 1 & length(method) != length(data)) {

      stop("The length of 'method' has to match the lenght of 'data'.",
           call. = FALSE)

    }

    av_distances <- get_kernel_info()

    if(any(!method %in% av_distances)) {

      stop('At least one unsupported distance method.',
           call. = FALSE)

    }

    stopifnot(is.numeric(weights))

    if(length(weights) != 1 & length(weights) != length(data)) {

      stop("The length of 'weights' has to match the lenght of 'data'.",
           call. = FALSE)

    }

    if(!rlang::is_function(sum)) {

      stop("'FUN' has to be a valid R function.",
           call. = FALSE)

    }

    ## distance matrix -------

    dist_mtx <- map2(data, method, calculate_dist)

    dist_mtx <- map2(dist_mtx, weights, `*`)

    reduce(dist_mtx, FUN, ...)

  }

#' @rdname calculate_dist
#' @export

  get_kernel_info <- function() {

    av_distances <- somKernels::get_kernel_info()

    av_distances[av_distances != 'BrayCurtis']

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

# Cross-validation folds -------

#' Create cross-validation folds for clustering.
#'
#' @description
#' Creates cross-validation folds in the single-layer (data frame or matrix)
#' or multi-layer clustering data.
#'
#' @details
#' The folds are generated with \code{\link[caret]{createFolds}}.
#' Intended for internal use.
#'
#' @param data the input data: a data frame, matrix or a list of data frames
#' or matrices.
#' @param k number of cross-validaiton folds.
#' @param seed an integer or `NULL` which specifies the seed for random
#' number generator.
#'
#' @return a list with two elements: `train` and `test` with the data subsets
#' intended for development of the clustering object and predictions,
#' respectively.

  create_clust_folds <- function(data, k, seed = 1234) {

    ## entry control is taken over by an upstream function

    stopifnot(is.numeric(k))
    k <- as.integer(k)

    set.seed(seed)

    ## single layer data -------

    if(is.data.frame(data) | is.matrix(data)) {

      fold_ids <- caret::createFolds(1:nrow(data), k = k)

      fold_set <- map(fold_ids,
                      ~list(train = data[-.x, ], test = data[.x, ]))

      return(purrr::transpose(fold_set))

    }

    ## multi-layer data ------

    fold_ids <- caret::createFolds(1:nrow(data[[1]]), k = k)

    if(is.null(names(data))) {

      data <- set_names(data,
                        paste0('layer_', 1:length(data)))

    }

    fold_set <- list()

    for(i in names(data)) {

      fold_set[['train']][[i]] <-
        map(fold_ids, ~data[[i]][-.x, ])

      fold_set[['test']][[i]] <-
        map(fold_ids, ~data[[i]][.x, ])

    }

    map(fold_set, purrr::transpose)

  }

# Permutation of clustering variables ------

#' Permute clustering variables.
#'
#' @description
#' Creates a series of clustering data sets with consecutive variables
#' reshuffled by random (permuted). Works for single- and multiple-layer data.
#'
#' @return a list of data frames or data lists.
#' The first element `data` stores the genuine data set.
#' Each subsequent element is named after the permuted variable.
#'
#' @param data a data frame, matrix or a list of such objects.

  permute_clust_data <- function(data) {

    ## handling single data frames or matrices -------

    if(is.data.frame(data) | is.matrix(data)) {

      if(is.null(colnames(data))) {

        colnames(data) <- paste0('V_', 1:ncol(data))

      }

      data <- as.data.frame(data)

      cl_feats <- names(data)

      noised_set <-
        map(cl_feats,
            ~mutate(data,
                    !!.x := sample(.data[[.x]],
                                   size = nrow(data),
                                   replace = TRUE)))

      noised_set <- set_names(noised_set, cl_feats)

      return(c(list(data = data), noised_set))

    }

    ## lists: formatting --------

    if(is.null(names(data))) {

      data <- set_names(data,
                        paste0('layer_', 1:length(data)))

    }

    cl_feats <- list()

    for(i in names(data)) {

      if(is.null(colnames(data[[i]]))) {

        colnames(data[[i]]) <- paste0(i, '_V', 1:ncol(data[[i]]))

      }

      data[[i]] <- data.frame(data[[i]])

      cl_feats[[i]] <- names(data[[i]])

    }

    all_vars <- unlist(cl_feats)

    if(any(duplicated(all_vars))) {

      stop('Duplicated variable names in the multi-layer data are not allowed.',
           call. = FALSE)

    }

    ## lists: permuting ------

    noised_set <- list()

    for(i in all_vars) {

      noised_set[[i]] <- data

      for(j in names(noised_set[[i]])) {

        if(!i %in% names(noised_set[[i]][[j]])) next

        noised_set[[i]][[j]] <-
          mutate(noised_set[[i]][[j]],
                 !!i := sample(.data[[i]],
                               size = nrow(noised_set[[i]][[j]]),
                               replace = TRUE))

      }
    }

    c(list(data = data), noised_set)

  }

# Multi-layer analysis ------

#' Check for multi-layer analysis.
#'
#' @description
#' Checks if the `clust_analysis` or `combi_analysis` object uses
#' multi-layer SOM.
#'
#' @param x an instance of the `clust_analysis` or `combi_analysis` class.
#'
#' @return a logical value

  is_multi_layer <- function(x) {

    if(is_clust_analysis(x)) {

      return(x$clust_fun %in% c('supersom', 'supersom_prediction'))

    } else {

      return(is_umatrix_analysis(x))

    }

  }

# END -----
