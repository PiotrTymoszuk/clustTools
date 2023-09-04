# Determines importance of clustering features

# The engine --------

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
#'
#' @export

  importance_cluster <- function(data,
                                 clustering_fun = kcluster,
                                 seed = 1234,
                                 .parallel = FALSE, ...) {

    ## entry control ------

    set.seed(seed)

    if(is_red_analysis(data)) {

      data <- column_to_rownames(data$component_tbl,
                                 'observation')

    }

    check_numeric(data)

    stopifnot(is.logical(.parallel))

    cl_feats <- names(data)

    ## noised set and the unnoised data ------

    noised_set <-
      purrr::map(cl_feats,
                 function(x) mutate(data,
                                    !!rlang::ensym(x) := sample(.data[[x]],
                                                                size = nrow(data),
                                                                replace = TRUE)))

    noised_set <- set_names(noised_set, cl_feats)

    noised_set <- c(list(data = data),
                    noised_set)

    ## generating the clustering objects, calculating the variances

    if(.parallel) {

      future::plan('multisession')

      var_lst <-
        furrr::future_map(noised_set,
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
                      variable = c('data', names(data)),
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
#' The permutation importance algorithm is 'blind' or agnostic to the
#' clustering procedure.
#' `impact()` is a S3 generic function.
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
#' @export

  impact <- function(x, ...) UseMethod('impact')

#' @rdname impact
#' @export impact.clust_analysis
#' @export

  impact.clust_analysis <- function(x,
                                    n_iter = 1,
                                    seed = 1234,
                                    .parallel = FALSE, ...) {

    ## entry control --------

    stopifnot(is_clust_analysis(x))
    stopifnot(is.numeric(n_iter))

    n_iter <- as.integer(n_iter)

    run <- NULL

    ## common parameters --------

    if(n_iter == 1) {

      cmm_args <-
        list2(data = eval_tidy(x$data),
              distance_method = x$dist_method,
              seed = seed,
              .parallel = .parallel,
              k = nrow(ngroups(x)),
              hc_method = x$hc_method,
              clust_fun = x$clust_fun,
              xdim = x$grid$xdim,
              ydim = x$grid$ydim,
              topo = x$grid$topo,
              eps = x$eps,
              minPts = x$minPts,
              neighbourhood.fct = as.character(x$grid$neighbourhood.fct),
              toroidal = x$grid$toroidal)

      cmm_args <- c(cmm_args, x$dots)

    } else {

      set.seed(seed)

      seeds <- sample(1:1e6, size = n_iter, replace = FALSE)

      seeds <- set_names(seeds, paste0('run_', 1:n_iter))

      cmm_args <-
        map(seeds,
            ~list2(data = eval_tidy(x$data),
                   distance_method = x$dist_method,
                   seed = .x,
                   .parallel = FALSE,
                   k = nrow(ngroups(x)),
                   hc_method = x$hc_method,
                   clust_fun = x$clust_fun,
                   xdim = x$grid$xdim,
                   ydim = x$grid$ydim,
                   topo = x$grid$topo,
                   eps = x$eps,
                   minPts = x$minPts,
                   neighbourhood.fct = as.character(x$grid$neighbourhood.fct),
                   toroidal = x$grid$toroidal))

      cmm_args <- map(cmm_args, ~c(.x, x$dots))

    }

    ## calls: no iterations ---------

    if(n_iter == 1) {

      if(x$clust_fun %in% c('hclust', 'som')) {

        cmm_args$clust_fun <- NULL

      } else if(x$clust_fun == 'dbscan') {

        cmm_args$k <- NULL

        cmm_args$clust_fun <- NULL

      }

      imp_call <-
        call2('importance_cluster',
              !!!compact(cmm_args),
              clustering_fun = switch(x$clust_fun,
                                      hclust = hcluster,
                                      kmeans = kcluster,
                                      pam = kcluster,
                                      dbscan = dbscan_cluster,
                                      som = som_cluster))

      return(eval(imp_call))

    }

    ## calls: multiple iterations ----------

    if(x$clust_fun %in% c('hclust', 'som')) {

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
                 clustering_fun = switch(x$clust_fun,
                                         hclust = hcluster,
                                         kmeans = kcluster,
                                         pam = kcluster,
                                         dbscan = dbscan_cluster,
                                         som = som_cluster)))

    if(.parallel) {

      future::plan('multisession')

      res <-
        furrr::future_map(imp_call,
                          eval,
                          .options = furrr::furrr_options(seed = TRUE))

      future::plan('sequential')

    } else {

      res <- map(imp_call, eval)

    }

    res <-
      map2_dfr(res, names(res),
               ~mutate(.x, run = .y))

    return(res)

  }

#' @rdname impact
#' @export impact.combi_analysis
#' @export

  impact.combi_analysis <- function(x,
                                    n_iter = 1,
                                    seed = 1234,
                                    .parallel = FALSE, ...) {

    ## entry control -------

    stopifnot(is_combi_analysis(x))
    stopifnot(is.numeric(n_iter))

    n_iter <- as.integer(n_iter)

    ## common parameters -------

    if(n_iter == 1) {

      args <-
        list2(data = eval_tidy(x$clust_analyses$observation$data),
              clustering_fun = combi_cluster,
              distance_som = x$clust_analyses$observation$dist_method,
              xdim = x$clust_analyses$observation$grid$xdim,
              ydim = x$clust_analyses$observation$grid$ydim,
              topo = x$clust_analyses$observation$grid$topo,
              neighbourhood.fct = as.character(x$clust_analyses$observation$grid$neighbourhood.fct),
              toroidal = x$clust_analyses$observation$grid$toroidal,
              rlen = nrow(x$clust_analyses$observation$clust_obj$changes),
              node_clust_fun = switch(x$clust_analyses$node$clust_fun,
                                      hclust = hcluster,
                                      kmeans = kcluster,
                                      pam = kcluster,
                                      dbscan = dbscan_cluster,
                                      som = som_cluster),
              seed = seed,
              .parallel = .parallel)

      args <- c(args, x$dots)

    } else {

      set.seed(seed)

      seeds <- sample(1:1e6, size = n_iter, replace = FALSE)

      seeds <- set_names(seeds, paste0('run_', 1:n_iter))

      args <-
        map(seeds,
            ~list2(data = eval_tidy(x$clust_analyses$observation$data),
                   clustering_fun = combi_cluster,
                   distance_som = x$clust_analyses$observation$dist_method,
                   xdim = x$clust_analyses$observation$grid$xdim,
                   ydim = x$clust_analyses$observation$grid$ydim,
                   topo = x$clust_analyses$observation$grid$topo,
                   neighbourhood.fct = as.character(x$clust_analyses$observation$grid$neighbourhood.fct),
                   toroidal = x$clust_analyses$observation$grid$toroidal,
                   rlen = nrow(x$clust_analyses$observation$clust_obj$changes),
                   node_clust_fun = switch(x$clust_analyses$node$clust_fun,
                                           hclust = hcluster,
                                           kmeans = kcluster,
                                           pam = kcluster,
                                           dbscan = dbscan_cluster,
                                           som = som_cluster),
                   seed = .x,
                   .parallel = FALSE))

      args <- map(args, ~c(.x, x$dots))

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

    if(.parallel) {

      future::plan('multisession')

      res <-
        furrr::future_map(imp_call,
                          eval,
                          .options = furrr::furrr_options(seed = TRUE))

      future::plan('sequential')

    } else {

      res <- map(imp_call, eval)

    }

    res <-
      map2_dfr(res, names(res),
               ~mutate(.x, run = .y))

    return(res)

  }

# END -----
