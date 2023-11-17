# Determines importance of clustering features

# Permutation importance --------

#' Determine clustering feature importance.
#'
#' @description
#' Determines importance of specific clustering variables by
#' comparing the fraction of 'explained' clustering variance of the input
#' clustering object and the object generated with the variable
#' re-shuffled randomly - i.e. so called 'permutation' importance.
#'
#' @references
#' Breiman L. Random forests. Mach Learn (2001) 45:5–32.
#' doi:10.1023/A:1010933404324
#'
#' @param data a numeric data frame, matrix or a red_analysis object. If a
#' red_analysis object is provided as the data argument, the observation
#' component/score table is subjected to clustering.
#' @param clustering_fun clustering function. Should return a
#' `clust_analysis` object.
#' @param seed initial setting of the random number generator.
#' @param .parallel logical, should the CV be run in parallel?
#' @param ... extra arguments passed to the clustering_fun.
#'
#' @return a data frame with the values of sum of squares and the clustering
#' variances.

  importance_cluster <- function(data,
                                 clustering_fun = kcluster,
                                 seed = 1234,
                                 .parallel = FALSE, ...) {

    ## entry control ------

    set.seed(seed)

    stopifnot(is.logical(.parallel))

    ## data control ------

    if(is.data.frame(data) | is.matrix(data)) {

      if(is_red_analysis(data)) {

        data <- column_to_rownames(data$component_tbl,
                                   'observation')

      }

      check_numeric(data)

    } else if(is.list(data)) {

      purrr::walk(data, check_numeric)

      if(is.null(names(data))) {

        data <- set_names(data, paste0('layer_', 1:length(data)))

      }

    } else {

      stop('Unsupported data format.', call. = FALSE)

    }

    ## noised set and the unnoised data ------

    noised_set <- permute_clust_data(data)

    if(.parallel) {

      future::plan('multisession')

      exports <- c('tidyverse',
                   'rlang',
                   'cluster',
                   'kohonen',
                   'dbscan',
                   'Rcpp',
                   'somKernels',
                   'stringi',
                   'nomclust',
                   'clustTools')

      var_lst <-
        furrr::future_map(noised_set,
                          function(x) clustering_fun(data = x, ...),
                          .options = furrr::furrr_options(seed = TRUE,
                                                          packages = exports))

      future::plan('sequential')

    } else {

      var_lst <- map(noised_set,
                     function(x) clustering_fun(data = x, ...))

    }

    var_lst <- map_dfr(var_lst,
                       ~as_tibble(var(.x)[c('total_wss',
                                            'total_ss',
                                            'between_ss',
                                            'frac_var')]))

    frac_var <- NULL

    importance(mutate(var_lst,
                      variable = names(noised_set),
                      frac_diff = var_lst$frac_var[1] - frac_var))

  }

# S3 generics and methods ------

#' Permutation importance of clustering features.
#'
#' @description
#' Determines importance of specific clustering variables by
#' comparing the fraction of 'explained' clustering variance of the input
#' clustering object and the object generated with the variable
#' re-shuffled randomly - so called 'permutation' importance.
#'
#' @details
#' `impact()` is a S3 generic function.
#' The permutation importance algorithm is 'blind' or agnostic to the
#' clustering procedure.
#' Note that it is not possible to compute clustering feature importance
#' for clustering analyses done with an user-provided dissimilarity objects
#' (subclass `min_analysis` of `clust_analysis`). In such cases, `NULL` is
#' returned with a warning.
#'
#' @references
#' Breiman L. Random forests. Mach Learn (2001) 45:5–32.
#' doi:10.1023/A:1010933404324
#'
#' @param x a `clust_analysis` object.
#' @param n_iter number of iterations, 1 by default.
#' If the arguments is larger that 1, the function is run multiple times,
#' which may help at testing variable importance in a more objective way
#' for different permutations.
#' @param seed initial setting of the random number generator.
#' @param .parallel logical, should the CV be run in parallel? Experimental.
#' @param ... extra arguments, currently none.
#'
#' @return a data frame of class \code{\link{importance}} with the defined
#' \code{\link{plot.importance}} and \code{\link{summary.importance}} methods.
#'
#' @export

  impact <- function(x, ...) UseMethod('impact')

#' @rdname impact
#' @export

  impact.clust_analysis <- function(x,
                                    n_iter = 1,
                                    seed = 1234,
                                    .parallel = FALSE, ...) {

    ## entry control --------

    start_time <- Sys.time()

    stopifnot(is_clust_analysis(x))
    stopifnot(is.numeric(n_iter))

    n_iter <- as.integer(n_iter)

    run <- NULL

    if(x$clust_fun != 'supersom') {

      distance_names <- x$dist_method

    } else {

      distance_names <- x$clust_obj$dist.fcts

    }

    clustering_fun <-
      switch(x$clust_fun,
             hclust = hcluster,
             kmeans = kcluster,
             htk = htk_cluster,
             pam = kcluster,
             dbscan = dbscan_cluster,
             som = som_cluster,
             supersom = som_cluster)

    ## parallel backend and benchmarking -------

    message(paste('Permutation importance testing with',
                  n_iter, 'iterations'))

    if(.parallel) future::plan('multisession')

    on.exit(message(paste('Elapsed:', Sys.time() - start_time)),
            add = TRUE)

    on.exit(future::plan('sequential'))

    ## common parameters --------

    cmm_list <-
      list2(data = model.frame(x),
            distance_method = distance_names,
            .parallel = .parallel,
            k = nrow(ngroups(x)),
            hc_method = x$hc_method,
            clust_fun = x$clust_fun,
            lambdas = x$lambdas,
            xdim = x$grid$xdim,
            ydim = x$grid$ydim,
            topo = x$grid$topo,
            eps = x$eps,
            minPts = x$minPts,
            neighbourhood.fct = as.character(x$grid$neighbourhood.fct),
            toroidal = x$grid$toroidal)

    if(n_iter == 1) {

      cmm_args <- c(cmm_list, list2(seed = seed))

      cmm_args <- c(cmm_args, x$dots)

    } else {

      set.seed(seed)

      seeds <- sample(1:1e6, size = n_iter, replace = FALSE)

      seeds <- set_names(seeds, paste0('run_', 1:n_iter))

      cmm_args <-
        map(seeds, ~c(cmm_list, list2(seed = .x)))

      cmm_args <- map(cmm_args, ~c(.x, x$dots))

    }

    ## calls: no iterations ---------

    if(n_iter == 1) {

      if(x$clust_fun %in% c('hclust', 'som', 'supersom', 'htk')) {

        cmm_args$clust_fun <- NULL

      } else if(x$clust_fun == 'dbscan') {

        cmm_args$k <- NULL

        cmm_args$clust_fun <- NULL

      }

      imp_call <-
        call2('importance_cluster',
              !!!compact(cmm_args),
              clustering_fun = clustering_fun)

      return(eval(imp_call))

    }

    ## calls: multiple iterations ----------

    if(x$clust_fun %in% c('hclust', 'som', 'supersom', 'htk')) {

      for(i in names(cmm_args)) {

        cmm_args[[i]]$clust_fun <- NULL

      }

    } else if(x$clust_fun == 'dbscan') {

      for(i in names(cmm_args)) {

        cmm_args[[i]]$k <- NULL

        cmm_args[[i]]$clust_fun <- NULL

      }

    }

    imp_call <-
      map(cmm_args,
          ~call2('importance_cluster',
                 !!!compact(.x),
                 clustering_fun = clustering_fun))

    res <-
      suppressMessages(furrr::future_map(imp_call,
                                         eval,
                                         .options = furrr::furrr_options(seed = TRUE)))

    res <-
      map2_dfr(res, names(res),
               ~mutate(.x, run = .y))

    return(res)

  }

#' @rdname impact
#' @export

  impact.min_analysis <- function(x, ...) {

    warning(paste("It is not possible to estimate variable importance",
                  "for analyses done with an user-provided",
                  "dissimilarity matrix."),
            call. = FALSE)

    NULL

  }

#' @rdname impact
#' @export

  impact.combi_analysis <- function(x,
                                    n_iter = 1,
                                    seed = 1234,
                                    .parallel = FALSE, ...) {

    ## entry control -------

    start_time <- Sys.time()

    stopifnot(is_combi_analysis(x))
    stopifnot(is.numeric(n_iter))

    n_iter <- as.integer(n_iter)

    if(is_umatrix_analysis(x)) {

      clustering_fun <- multi_cluster

      distance_som <- NULL

      distance_method <- x$clust_analyses$observation$clust_obj$dist.fcts

    } else {

      clustering_fun <- combi_cluster

      distance_som <- x$clust_analyses$observation$dist_method

      distance_method <- NULL

    }

    node_clust_fun = switch(x$clust_analyses$node$clust_fun,
                            hclust = hcluster,
                            kmeans = kcluster,
                            htk = htk_cluster,
                            pam = kcluster,
                            dbscan = dbscan_cluster,
                            som = som_cluster)

    ## parallel backend and benchmarking --------

    message(paste('Permutation importance testing with',
                  n_iter, 'iterations'))

    if(.parallel) future::plan('multisession')

    on.exit(message(paste('Elapsed:', Sys.time() - start_time)),
            add = TRUE)

    on.exit(future::plan('sequential'))

    ## common parameters -------

    cmm_list <-
      list2(data = model.frame(x)$observation,
            clustering_fun = clustering_fun,
            distance_som = distance_som,
            distance_method = distance_method,
            xdim = x$clust_analyses$observation$grid$xdim,
            ydim = x$clust_analyses$observation$grid$ydim,
            topo = x$clust_analyses$observation$grid$topo,
            neighbourhood.fct = as.character(x$clust_analyses$observation$grid$neighbourhood.fct),
            toroidal = x$clust_analyses$observation$grid$toroidal,
            rlen = nrow(x$clust_analyses$observation$clust_obj$changes),
            node_clust_fun = node_clust_fun,
            .parallel = .parallel)

    if(n_iter == 1) {

      args <- c(cmm_list, list(seed = seed))

      args <- c(args, x$dots)

      args <- compact(args)

    } else {

      set.seed(seed)

      seeds <- sample(1:1e6, size = n_iter, replace = FALSE)

      seeds <- set_names(seeds, paste0('run_', 1:n_iter))

      cmm_list$.parallel <- FALSE

      args <-
        map(seeds, ~c(cmm_list, list(seed = .x)))

      args <- map(args, ~c(.x, x$dots))

      args <- map(args, compact)

    }

    ## calls and evaluation: one iteration -------

    if(n_iter == 1) {

      imp_call <- call2('importance_cluster',
                        !!!args)

      return(eval(imp_call))

    }

    ## calls and evaluation: multiple iterations

    imp_call <-
      map(args, ~call2('importance_cluster', !!!.x))

    exports <- c('tidyverse',
                 'rlang',
                 'cluster',
                 'kohonen',
                 'dbscan',
                 'Rcpp',
                 'somKernels',
                 'stringi',
                 'nomclust',
                 'clustTools')

    res <-
      suppressMessages(furrr::future_map(imp_call,
                                         eval,
                                         .options = furrr::furrr_options(seed = TRUE,
                                                                         packages = exports)))

    res <-
      map2_dfr(res, names(res),
               ~mutate(.x, run = .y))

    return(res)

  }

# END -----
