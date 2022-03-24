# Determines importance of clustering features

#' Determine clustering feature importance.
#'
#' @description Determines importance of specific clustering variables by
#' comparing the fraction of 'explained' clustering variance of the input
#' clustering object and the object generated with the variable
#' re-shuffled randomly.
#' @param data a numeric data frame, matrix or a red_analysis object. If a
#' red_analysis object is provided as the data argument, the observation
#' component/score table is subjected to clustering.
#' @param clustering_fun clustering function. Should return a
#' clust_analysis object.
#' @param seed initial setting of the random number generator.
#' @param .parallel logical, should the CV be run in parallel? Experimental.
#' @param ... extra arguments passed to the clustering_fun.
#' @return a data frame with the values of sum of squares and the clustering
#' variances.
#' @export

  importance_cluster <- function(data,
                                 clustering_fun = kcluster,
                                 seed = 1234,
                                 .parallel = FALSE, ...) {

    ## entry control

    if(all(class(data) == 'red_analysis')) {

      data <- tibble::column_to_rownames(data$component_tbl,
                                         'observation')

    }

    clustTools:::check_numeric(data)

    stopifnot(is.logical(.parallel))

    cl_feats <- names(data)

    ## noised set and the unnoised data

    noised_set <- purrr::map(cl_feats,
                             function(x) dplyr::mutate(data,
                                                       !!rlang::ensym(x) := sample(.data[[x]],
                                                                            size = nrow(data),
                                                                            replace = TRUE)))

    noised_set <- rlang::set_names(noised_set, cl_feats)

    noised_set <- c(list(data = data),
                    noised_set)

    ## generating the clustering objects, calculating the variances

    if(.parallel) {

      future::plan('multisession')

      var_lst <-  furrr::future_map(noised_set,
                                    function(x) clustering_fun(data = x, ...),
                                    .options = furrr::furrr_options(seed = TRUE,
                                                                    packages = c('tidyverse',
                                                                                 'rlang',
                                                                                 'cluster',
                                                                                 'kohonen',
                                                                                 'dbscan',
                                                                                 'Rcpp',
                                                                                 'somKernels',
                                                                                 'stringi',
                                                                                 'nomclust',
                                                                                 'clustTools')))

      future::plan('sequential')

    } else {

      var_lst <- purrr::map(noised_set,
                            function(x) clustering_fun(data = x, ...))

    }

    var_lst <- purrr::map_dfr(var_lst,
                              ~tibble::as_tibble(var(.x)[c('total_wss',
                                                           'total_ss',
                                                           'between_ss',
                                                           'frac_var')]))

    clustTools::importance(dplyr::mutate(var_lst,
                                         variable = c('data', names(data)),
                                         frac_diff = var_lst$frac_var[1] - frac_var))

  }

# END -----
