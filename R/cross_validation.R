# Stability of the clustering structure by cross-validation.

# The engine -----------

#' Cross-validate a clustering algorithm.
#'
#' @description
#' Checks the stability of a clustering solution by
#' cross-validation (CV) and the classification error as a measure of
#' the cluster stability.
#'
#' @details
#' By principle, similar to cross-validation of any machine learning
#' multi-level classifier. The training portion of a CV split is used to develop
#' of a cluster structure and the projection on the test portion is accomplished
#' by k-nearest neighbor (kNN) label propagation algorithm. For its
#' implementation details, see: \code{\link{propagate}}.
#' The fold are generated with \code{\link[caret]{createFolds}}.
#'
#' @references
#' Lange T, Roth V, Braun ML, Buhmann JM. Stability-based validation of
#' clustering solutions. Neural Comput (2004) 16:1299–1323.
#' doi:10.1162/089976604773717621
#' @references
#' Leng M, Wang J, Cheng J, Zhou H, Chen X. Adaptive semi-supervised
#' clustering algorithm with label propagation. J Softw Eng (2014) 8:14–22.
#' doi:10.3923/jse.2014.14.22
#' @references
#' Kuhn M. Building predictive models in R using the caret package.
#' J Stat Softw (2008) 28:1–26. doi:10.18637/jss.v028.i05
#'
#' @param data a numeric data frame, matrix or a `red_analysis` object. If a
#' `red_analysis` object is provided as the data argument, the observation
#' component/score table is subjected to clustering.
#' @param nfolds number of CV folds.
#' @param kNN number of the nearest neighbors.
#' @param simple_vote logical, should classical unweighted k-NN classification
#' be applied? If FALSE, distance-weighted k-NN is used with the provided kernel
#' function.
#' @param resolve_ties logical, should the ties be resolved at random? Applies
#' only to the simple unweighted voting algorithm.
#' @param kernel_fun kernel function transforming the distance into weight.
#' @param clustering_fun clustering function. Should return a
#' `clust_analysis` object.
#' @param seed initial setting of the random number generator.
#' @param .parallel logical, should the CV be run in parallel?
#' @param ... extra arguments passed to the clustering_fun.
#'
#' @return a list containing the global \code{\link{clust_analysis}}
#' object, projection (prediction) results and
#' prediction summary for each fold and a prediction
#' summary for the whole CV.
#'
#' @export

  cv_cluster <- function(data,
                         nfolds = 5,
                         kNN = 5,
                         simple_vote = TRUE,
                         resolve_ties = FALSE,
                         kernel_fun = function(x) 1/x,
                         clustering_fun = clustTools::kcluster,
                         seed = 1234,
                         .parallel = FALSE, ...) {

    ## entry control --------

    if(is_red_analysis(data)) {

      data <- column_to_rownames(data$component_tbl,
                                 'observation')

    }

    check_numeric(data)

    nfolds <- as.integer(nfolds)

    kNN <- as.integer(kNN)

    stopifnot(is.logical(simple_vote))
    stopifnot(is.logical(resolve_ties))
    stopifnot(is.logical(.parallel))

    stopifnot(is.function(kernel_fun))

    global_clust <- NULL
    fold_clust <- NULL
    correct <- NULL

    ## tests stability of the clustering algorithm by cross validation --------

    start_time <- Sys.time()
    message(paste('CV: ', nfolds, 'folds'))
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))

    set.seed(seed = seed)

    ## data split generation

    set.seed(seed)

    fold_ids <- caret::createFolds(1:nrow(data), k = nfolds)

    fold_set <- map(fold_ids, ~list(train = data[-.x, ],
                                    test = data[.x, ]))

    fold_set <- purrr::transpose(fold_set)

    ## defining the global classifier

    glob_classif <- clustering_fun(data = data, ...)

    glob_assign <- glob_classif$clust_assignment[c('observation', 'clust_id')]

    glob_assign <- set_names(glob_assign,
                             c('observation', 'global_clust'))

    ## creating the training set classifiers

    if(.parallel) {

      future::plan('multisession')

      train_classif <-
        furrr::future_map(fold_set$train ,
                          function(x) clustering_fun(data = x, ...),
                          .options = furrr::furrr_options(seed = TRUE,
                                                          packages = c('tidyverse',
                                                                       'rlang',
                                                                       'cluster',
                                                                       'kohonen',
                                                                       'dbscan',
                                                                       'Rcpp',
                                                                       'somKernels',
                                                                       'clustTools',
                                                                       'stringi',
                                                                       'nomclust')))

     future:: plan('sequential')

    } else {

      train_classif <- map(fold_set$train,
                           function(x) clustering_fun(data = x, ...))

    }

    ## obtaining the predictions and comparing with the global classifier

    test_preds <- pmap(list(object = train_classif,
                            newdata = fold_set$test),
                       purrr::safely(predict),
                       type = 'propagation',
                       kNN = kNN,
                       simple_vote = simple_vote,
                       resolve_ties = resolve_ties,
                       kernel_fun = kernel_fun)

    test_preds <- map(test_preds, ~.x$result)

    test_preds <- compact(test_preds)

    test_preds <-
      map(test_preds,
          ~set_names(.x$clust_assignment[c('observation', 'clust_id')],
                     c('observation', 'fold_clust')))

    test_preds <- map(test_preds,
                      ~left_join(.x,
                                 glob_assign,
                                 by = 'observation'))

    test_preds <-
      map(test_preds,
          mutate,
          correct = as.character(global_clust) == as.character(fold_clust))

    ## summary stats

    test_stats <-
      map(test_preds,
          dplyr::summarise,
          corr_rate = mean(as.numeric(correct),
                           na.rm = TRUE))

    test_stats <- map2_dfr(test_stats,
                           names(test_stats),
                           ~mutate(.x,
                                   err_rate = 1 - corr_rate,
                                   fold = .y))

    test_stats <- test_stats[c('fold', 'corr_rate', 'err_rate')]

    bca_err <- coxed::bca(test_stats$err_rate,
                          conf.level = 0.95)

    test_summary <-
      tibble(mean_error = mean(test_stats$err_rate,
                               na.rm = TRUE),
             lower_ci = bca_err[1],
             upper_ci = bca_err[2])

    test_preds <-
      map2_dfr(test_preds,
               names(test_preds),
               ~mutate(.x, fold = .y))

    ## output

    list(clust_analysis_object = glob_classif,
         predictions = test_preds,
         fold_stats = test_stats,
         summary = test_summary)

  }

# S3 generics and methods --------

#' Cross-validate the clustering analysis object.
#'
#' @description
#' Checks the stability of a clustering solution by
#' cross-validation (CV) and the classification error in CV folds
#' as a measure of the cluster stability.
#'
#' @details
#' By principle similar to cross-validation of any machine learning
#' multi-class classifier. The training portion of a CV split is used to develop
#' of a cluster structure and the projection on the test portion is accomplished
#' by k-nearest neighbor (kNN) label propagation algorithm. For its
#' implementation details, see: \code{\link{propagate}}.
#' For `combi_analysis` objects, assignment of the observations in the CV folds
#' is done for the 'top' assignment of the observations to the clusters;
#' nodes are ignored!
#' The fold are generated with \code{\link[caret]{createFolds}}.
#' `cv()` is a S3 generic function.
#'
#' @references
#' Lange T, Roth V, Braun ML, Buhmann JM. Stability-based validation of
#' clustering solutions. Neural Comput (2004) 16:1299–1323.
#' doi:10.1162/089976604773717621
#' @references
#' Leng M, Wang J, Cheng J, Zhou H, Chen X. Adaptive semi-supervised
#' clustering algorithm with label propagation. J Softw Eng (2014) 8:14–22.
#' doi:10.3923/jse.2014.14.22
#' @references
#' Kuhn M. Building predictive models in R using the caret package.
#' J Stat Softw (2008) 28:1–26. doi:10.18637/jss.v028.i05
#'
#' @param x an object.
#' @param nfolds number of CV folds.
#' @param kNN number of the nearest neighbors.
#' @param simple_vote logical, should classical unweighted k-NN classification
#' be applied? If FALSE, distance-weighted k-NN is used with the provided kernel
#' function.
#' @param resolve_ties logical, should the ties be resolved at random? Applies
#' only to the simple unweighted voting algorithm.
#' @param kernel_fun kernel function transforming the distance into weight.
#' @param seed initial setting of the random number generator.
#' @param .parallel logical, should the CV be run in parallel?
#' @param ... extra arguments, currently none.
#'
#' @return a list containing the global \code{\link{clust_analysis}} object,
#' projection (prediction) results and prediction summary
#' for each fold and a prediction
#' summary for the whole CV.
#'
#' @export

  cv <- function(x, ...) UseMethod('cv')

#' @rdname cv
#' @export cv.clust_analysis
#' @export

  cv.clust_analysis <- function(x,
                                nfolds = 5,
                                kNN = 5,
                                simple_vote = TRUE,
                                resolve_ties = FALSE,
                                kernel_fun = function(x) 1/x,
                                seed = 1234,
                                .parallel = FALSE, ...) {

    ## entry control ---------

    stopifnot(is_clust_analysis(x))

    nfolds <- as.integer(nfolds)

    kNN <- as.integer(kNN)

    stopifnot(is.logical(simple_vote))
    stopifnot(is.logical(resolve_ties))
    stopifnot(is.logical(.parallel))

    stopifnot(is.function(kernel_fun))

    ## common parameters --------

    cmm_args <-
      list2(data = eval_tidy(x$data),
            nfolds = nfolds,
            kNN = kNN,
            simple_vote = simple_vote,
            resolve_ties = resolve_ties,
            kernel_fun = kernel_fun,
            distance_method = x$dist_method,
            seed = seed,
            .parallel = .parallel,
            k = nrow(ngroups(x)),
            hc_method = x$hc_method,
            clust_fun = x$clust_fun,
            xdim = x$grid$xdim,
            ydim = x$grid$ydim,
            topo = x$grid$topo,
            neighbourhood.fct = as.character(x$grid$neighbourhood.fct),
            toroidal = x$grid$toroidal,
            eps = x$eps,
            minPts = x$minPts)

    cmm_args <- c(cmm_args, x$dots)

    ## function calls -------

    if(x$clust_fun %in% c('hclust', 'som')) {

      cmm_args$clust_fun <- NULL

    } else if(x$clust_fun == 'dbscan') {

      cmm_args$k <- NULL

      cmm_args$clust_fun <- NULL

    }

    cv_call <-
      call2('cv_cluster',
            !!!compact(cmm_args),
            clustering_fun = switch(x$clust_fun,
                                    hclust = hcluster,
                                    kmeans = kcluster,
                                    pam = kcluster,
                                    dbscan = dbscan_cluster,
                                    som = som_cluster))

    eval(cv_call)

  }

#' @rdname cv
#' @export cv.combi_analysis
#' @export

  cv.combi_analysis <- function(x,
                                nfolds = 5,
                                kNN = 5,
                                simple_vote = TRUE,
                                resolve_ties = FALSE,
                                kernel_fun = function(x) 1/x,
                                seed = 1234,
                                .parallel = FALSE, ...) {

    ## entry control -------

    stopifnot(is_combi_analysis(x))

    nfolds <- as.integer(nfolds)

    kNN <- as.integer(kNN)

    stopifnot(is.logical(simple_vote))
    stopifnot(is.logical(resolve_ties))
    stopifnot(is.logical(.parallel))

    stopifnot(is.function(kernel_fun))

    ## cross-validation ---------

    args <-
      list2(data = eval_tidy(x$clust_analyses$observation$data),
            nfolds = 10,
            kNN = 5,
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
            kernel_fun = kernel_fun,
            simple_vote = simple_vote,
            resolve_ties = resolve_ties,
            seed = seed,
            .parallel = .parallel)

    args <- c(args, x$dots)

    cv_call <- call2('cv_cluster',
                     !!!args)

    eval(cv_call)

  }

# END -----
