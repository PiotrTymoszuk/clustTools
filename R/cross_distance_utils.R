# Cross_distance calculation helpers

#' Cross-distances calculation helpers.
#'
#' @description
#' The functions compute homologous and heterologous distance between the
#' clusters of a single or two clustering analysis objects, respectively.
#'
#' @return an instance of the `cross_dist` class.
#'
#' @details
#' Designed solely for internal use. `cross_single_homolog()` and
#' `cross_single_heterolog()` compute cross-distances for single layer data.
#' `cross_multi_homolog()` and `cross_multi_heterolog()` tackle multi-layer SOM.
#' In multi-layer analyses, the cross-distance matrix is a weighted sum of
#' cross-distances between the corresponding data layers. The weights are
#' extracted from the SOM analysis object.
#'
#' @param x a `clust_analysis` or `combi_analysis` object.
#' @param y a `clust_analysis` or `combi_analysis` object.
#' @param method name of the distance metric. If NULL, it will be extracted
#' from the `x` object.

  cross_single_homolog <- function(x, method) {

    ## entry control -------

    err_txt <- paste("'x' has to be an instance of the 'clust_analysis'",
                     "or 'combi_analysis' class.")

    if(!is_clust_analysis(x) & !is_combi_analysis(x)) {

      stop(err_txt, call. = FALSE)

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

    ## distance computation -------

    clust_ass <- extract(x, 'assignment')

    clust_ass <- split(clust_ass$observation, clust_ass$clust_id)

    if(is_clust_analysis(x)) {

      clust_data <- model.frame(x)

    } else {

      clust_data <- model.frame(x)$observation

    }

    clust_data <-
      map(clust_ass,
          ~clust_data[.x, , drop = FALSE])

    clust_data <- map(clust_data, as.data.frame)

    pairs <- utils::combn(names(clust_data), m = 2, simplify = FALSE)

    identities <- map(names(clust_data), ~c(.x, .x))

    pairs <- c(pairs, identities)

    pair_names <- map_chr(pairs, paste, collapse = ' vs ')

    pairs <- set_names(pairs, pair_names)

    dist_lst <-
      furrr::future_map(pairs,
                        ~cross_distance(clust_data[[.x[1]]],
                                        clust_data[[.x[2]]],
                                        method = method),
                        .options = furrr::furrr_options(seed = TRUE))

    cross_dist(dist_lst,
               type = 'homologous',
               method = method,
               x_levels = names(clust_ass),
               y_levels = names(clust_ass))

  }

#' @rdname cross_single_homolog

  cross_single_heterolog <- function(x, y, method) {

    ## entry control -------

    err_txt <- paste("'x' has to be an instance of the 'clust_analysis'",
                     "or 'combi_analysis' class.")

    if(!is_clust_analysis(x) & !is_combi_analysis(x)) {

      stop(err_txt, call. = FALSE)

    }

    if(!is.null(y)) {

      if(!is_clust_analysis(y) & !is_combi_analysis((y))) {

        stop(err_txt, call. = FALSE)

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

    ## distance computation ------

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

    clust_data <- map(clust_data, as.data.frame)

    cmm_vars <- map(clust_data, names)

    cmm_vars <- reduce(cmm_vars, intersect)

    if(length(cmm_vars) == 0) {

      stop("No common variables for 'x' and 'y'.", call. = FALSE)

    }

    clust_data <- map(clust_data, ~.x[cmm_vars])

    clust_data <- map2(clust_data, clust_ass,
                       function(data, clust) map(clust,
                                                 ~data[.x, , drop = FALSE]))


    dist_lst <-
      furrr::future_map(pairs,
                        ~cross_distance(clust_data$x[[.x[1]]],
                                        clust_data$y[[.x[2]]],
                                        method = method),
                        .options = furrr::furrr_options(seed = TRUE))

    cross_dist(dist_lst,
               type = 'heterologous',
               method = method,
               x_levels = names(clust_ass[['x']]),
               y_levels = names(clust_ass[['y']]))

  }

#' @rdname cross_single_homolog

  cross_multi_homolog <- function(x) {

    ## entry control ------

    err_txt <- paste("'x' has to be an instance of the 'clust_analysis'",
                     "or 'combi_analysis' class.")

    if(!is_clust_analysis(x)) {

      if(!is_combi_analysis(x)) {

        stop(err_txt, call. = FALSE)
      }

    }

    ## distance metrics, weights ------

    if(is_clust_analysis(x)) {

      if(x$clust_fun == 'supersom') {

        distance_names <- x$clust_obj$dist.fcts

        weights <-
          x$clust_obj$user.weights * x$clust_obj$distance.weights

      } else {

        distance_names <- x$dots$dist.fcts

        weights <-
          x$dots$user.weights * x$dots$distance.weights

      }

    } else {

      distance_names <- x$clust_analyses$observation$clust_obj$dist.fcts

      weights <-
        x$clust_analyses$observation$clust_obj$user.weights *
        x$clust_analyses$observation$clust_obj$distance.weights

    }

    ## cluster pairs and data layers -------

    clust_ass <- extract(x, 'assignment')

    clust_ass <- split(clust_ass$observation, clust_ass$clust_id)

    if(is_clust_analysis(x)) {

      data_layers <- model.frame(x)

    } else {

      data_layers <- model.frame(x)$observation

    }

    data_layers <- map(data_layers, data.frame)

    for(i in seq_along(data_layers)) {

      data_layers[[i]] <-
        map(clust_ass,
            ~data_layers[[i]][.x, , drop = FALSE])

    }

    pairs <- utils::combn(names(clust_ass), m = 2, simplify = FALSE)

    identities <- map(names(clust_ass), ~c(.x, .x))

    pairs <- c(pairs, identities)

    pair_names <- map_chr(pairs, paste, collapse = ' vs ')

    pairs <- set_names(pairs, pair_names)

    ## distance computation -------

    dist_lst <- list()

    for(i in seq_along(data_layers)) {

      dist_lst[[i]] <-
        furrr::future_map(pairs,
                          ~cross_distance(data_layers[[i]][[.x[1]]],
                                          data_layers[[i]][[.x[2]]],
                                          method = distance_names[[i]]),
                          .options = furrr::furrr_options(seed = TRUE))

    }

    if(!is.null(names(data_layers))) {

      dist_lst <- set_names(dist_lst, names(data_layers))

    } else {

      dist_lst <- set_names(data_layers,
                            paste0('layer_', 1:length(dist_lst)))

    }

    dist_lst <- purrr::transpose(dist_lst)

    ## weighted sum ------

    for(i in seq_along(dist_lst)) {

      dist_lst[[i]] <-
        map2(dist_lst[[i]], weights, `*`)

      dist_lst[[i]] <- reduce(dist_lst[[i]], `+`)

    }

    ## removal of computation failures

    dist_lst <- map(dist_lst, function(x) if(all(is.na(x))) NULL else x)

    dist_lst <- compact(dist_lst)

    cross_dist(dist_lst,
               type = 'homologous',
               method = 'weighted_som',
               x_levels = names(clust_ass),
               y_levels = names(clust_ass))

  }

#' @rdname cross_single_homolog

  cross_multi_heterolog <- function(x, y) {

    # entry control -----

    err_txt <- paste("'x' and 'y' have to be instances of the 'clust_analysis'",
                     "or 'combi_analysis' class.")

    if(!is_clust_analysis(x)) {

      if(!is_combi_analysis(x)) {

        stop(err_txt, call. = FALSE)

      }

    }

    if(!is_clust_analysis(y)) {

      if(!is_combi_analysis(y)) {

        stop(err_txt, call. = FALSE)

      }

    }

    ## data layer consistency and names ------

    clust_list <- list(x = x, y = y)

    data_layers <-
      map(clust_list,
          function(x) if(is_clust_analysis(x)) model.frame(x) else model.frame(x)$observation)

    if(length(data_layers[[1]]) != length(data_layers[[2]])) {

      stop("Sizes of data layers in 'x' and 'y' differ.",
           call. = FALSE)

    }

    for(i in seq_along(data_layers)) {

      if(is.null(names(data_layers[[i]]))) {

        data_layers[[i]] <-
          set_names(data_layers[[i]],
                    paste0('layer_', 1:length(data_layers[[i]])))

      }

      data_layers[[i]] <- map(data_layers[[i]], as.data.frame)

    }

    cmm_names <-
      map2(data_layers[['x']],
           data_layers[['y']],
           ~intersect(names(.x), names(.y)))

    if(any(map_dbl(cmm_names, length) == 0)) {

      stop(paste("No common variables in at least one data layer.",
                 "Sure you are comparing the right objects?"),
           call. = FALSE)

    }

    for(i in seq_along(data_layers)) {

      data_layers[[i]] <-
        map2(data_layers[[i]], cmm_names, ~.x[.y])

    }

    ## distance metrics and weights -------

    if(is_clust_analysis(x)) {

      if(clust_list[['x']]$clust_fun == 'supersom') {

        distance_names <- clust_list[['x']]$clust_obj$dist.fcts

        weights <-
          clust_list[['x']]$clust_obj$user.weights * clust_list[['x']]$clust_obj$distance.weights

      } else {

        distance_names <- clust_list[['x']]$dots$dist.fcts

        weights <-
          clust_list[['x']]$dots$user.weights * clust_list[['x']]$dots$distance.weights

      }

    } else {

      distance_names <- x$clust_analyses$observation$clust_obj$dist.fcts

      weights <-
        x$clust_analyses$observation$clust_obj$user.weights *
        x$clust_analyses$observation$clust_obj$distance.weights

    }

    ## cluster pairs ------

    clust_ass <- map(clust_list, extract, 'assignment')

    clust_ass <- map(clust_ass,
                     ~split(.x$observation, .x$clust_id))

    for(i in seq_along(data_layers[[1]])) {

      data_layers[['x']][[i]] <-
        map(clust_ass[['x']], ~data_layers[['x']][[i]][.x, , drop = FALSE])

      data_layers[['y']][[i]] <-
        map(clust_ass[['y']], ~data_layers[['y']][[i]][.x, , drop = FALSE])

    }

    pairs <-
      map(names(clust_ass[['x']]),
          function(x) map(names(clust_ass[['y']]), ~c(x, .x)))

    pairs <- unlist(pairs, recursive = FALSE)

    pair_names <- map(pairs, paste, collapse = ' vs ')

    pairs <- set_names(pairs, pair_names)

    ## cross-distance matrices for the layers ---------

    dist_lst <- list()

    for(i in seq_along(data_layers[[1]])) {

      dist_lst[[i]] <-
        furrr::future_map(pairs,
                          ~cross_distance(data_layers$x[[i]][[.x[1]]],
                                          data_layers$y[[i]][[.x[2]]],
                                          method = distance_names[[i]]),
                          .options = furrr::furrr_options(seed = TRUE))

    }

    dist_lst <- set_names(dist_lst, paste0('layer_', 1:length(dist_lst)))

    dist_lst <- purrr::transpose(dist_lst)

    ## weighted distances and output -------

    for(i in seq_along(dist_lst)) {

      dist_lst[[i]] <-
        map2(dist_lst[[i]], weights, `*`)

      dist_lst[[i]] <- reduce(dist_lst[[i]], `+`)

    }

    ## removal of computation failures

    dist_lst <- map(dist_lst, function(x) if(all(is.na(x))) NULL else x)

    dist_lst <- compact(dist_lst)

    cross_dist(dist_lst,
               type = 'heterologous',
               method = 'weighted_som',
               x_levels = names(clust_ass[['x']]),
               y_levels = names(clust_ass[['y']]))

  }

# END ------
