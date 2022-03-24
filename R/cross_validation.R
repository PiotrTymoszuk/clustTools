# Stability of the clustering structure by cross-validation.

#' Cross-validate a clustering algorithm.
#'
#' @description Checks the stability of a clustering solution by
#' cross-validation (CV) and the classification error as a measure of the cluster
#' stability.
#' @details By principle similar to cross-validation of any machine learning
#' multi-level classifier. The training portion of a CV split is used to develop
#' of a cluster structure and the projection on the test portion is accomplished
#' by k-nearest neighbor (kNN) label propagation algorithm. For its
#' implementation details, see: \code{\link{propagate}}.
#' The fold are generated with \code{\link[caret]{createFolds}}.
#' @param data a numeric data frame, matrix or a red_analysis object. If a
#' red_analysis object is provided as the data argument, the observation
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
#' clust_analysis object.
#' @param seed initial setting of the random number generator.
#' @param .parallel logical, should the CV be run in parallel? Experimental.
#' @param ... extra arguments passed to the clustering_fun.
#' @return a list containing the global clust_analysis object, projection
#' (prediction) results and prediction summary for each fold and a prediction
#' summary for the whole CV.
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

    ## entry control

    if(all(class(data) == 'red_analysis')) {

      data <- tibble::column_to_rownames(data$component_tbl,
                                         'observation')

    }

    clustTools:::check_numeric(data)

    nfolds <- as.integer(nfolds)

    kNN <- as.integer(kNN)

    stopifnot(is.logical(simple_vote))
    stopifnot(is.logical(resolve_ties))
    stopifnot(is.logical(.parallel))

    stopifnot(is.function(kernel_fun))

    ## tests stability of the clustering algorithm by cross validation

    start_time <- Sys.time()
    message(paste('CV: ', nfolds, 'folds'))
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))

    set.seed(seed = seed)

    ## data split generation

    set.seed(seed)

    fold_ids <- caret::createFolds(1:nrow(data), k = nfolds)

    fold_set <- purrr::map(fold_ids, ~list(train = data[-.x, ],
                                           test = data[.x, ]))

    fold_set <- purrr::transpose(fold_set)

    ## defining the global classifier

    glob_classif <- clustering_fun(data = data, ...)

    glob_assign <- glob_classif$clust_assignment[c('observation', 'clust_id')]

    glob_assign <- rlang::set_names(glob_assign,
                                    c('observation', 'global_clust'))

    ## creating the training set classifiers

    if(.parallel) {

      future::plan('multisession')

      train_classif <- furrr::future_map(fold_set$train ,
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

      train_classif <- purrr::map(fold_set$train,
                                  function(x) clustering_fun(data = x, ...))

    }

    ## obtaining the predictions and comparing with the global classifier

    test_preds <- purrr::pmap(list(object = train_classif,
                                   newdata = fold_set$test),
                              purrr::safely(predict),
                              type = 'propagation',
                              kNN = kNN,
                              simple_vote = simple_vote,
                              resolve_ties = resolve_ties,
                              kernel_fun = kernel_fun)

    test_preds <- purrr::map(test_preds, ~.x$result)

    test_preds <- purrr::compact(test_preds)

    test_preds <- purrr::map(test_preds,
                             ~rlang::set_names(.x$clust_assignment[c('observation', 'clust_id')],
                                               c('observation', 'fold_clust')))

    test_preds <- purrr::map(test_preds,
                             ~dplyr::left_join(.x,
                                               glob_assign,
                                               by = 'observation'))

    test_preds <- purrr::map(test_preds,
                             dplyr::mutate,
                             correct = as.character(global_clust) == as.character(fold_clust))

    ## summary stats

    test_stats <- purrr::map(test_preds,
                             dplyr::summarise,
                             corr_rate = mean(as.numeric(correct),
                                              na.rm = TRUE))

    test_stats <- purrr::map2_dfr(test_stats,
                                  names(test_stats),
                                  ~dplyr::mutate(.x,
                                                 err_rate = 1 - corr_rate,
                                                 fold = .y))

    test_stats <- test_stats[c('fold', 'corr_rate', 'err_rate')]

    bca_err <- coxed::bca(test_stats$err_rate,
                          conf.level = 0.95)

    test_summary <- tibble::tibble(mean_error = mean(test_stats$err_rate,
                                                     na.rm = TRUE),
                                   lower_ci = bca_err[1],
                                   upper_ci = bca_err[2])

    test_preds <- purrr::map2_dfr(test_preds,
                                  names(test_preds),
                                  ~dplyr::mutate(.x, fold = .y))

    ## output

    list(clust_analysis_object = glob_classif,
         predictions = test_preds,
         fold_stats = test_stats,
         summary = test_summary)

  }

# END -----
