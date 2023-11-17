# Stability of the clustering structure by cross-validation.

# Out-of fold predictions -----------

#' Cross-validate a clustering algorithm.
#'
#' @description
#' Checks the quality of a clustering solution by
#' cross-validation (CV) with k-nearest neighbors (kNN) out-of-fold predictions
#' or predictions made by a self-organizing map (SOM).
#' Stability of the clustering structure is measured by cluster assignment
#' classification error in the out-of-fold predictions as compared with the
#' genuine clustering structure. Explanatory value and cluster separation are
#' determined by clustering variance and silhouette statistics.
#'
#' @details
#' By principle, similar to cross-validation of any machine learning
#' multi-class classifier. The training portion of a CV split is used to develop
#' of a cluster structure and the projection on the test portion is accomplished
#' by k-nearest neighbor (kNN) label propagation algorithm or derives the
#' cluster assignment from a trained SOM.
#' For implementation details, see: \code{\link{propagate}}
#' and \code{\link{map_som}}.
#' The folds are generated with \code{\link[caret]{createFolds}}.
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
#' @references
#' Rousseeuw PJ. Silhouettes: A graphical aid to the interpretation and
#' validation of cluster analysis. J Comput Appl Math (1987) 20:53–65.
#' doi:10.1016/0377-0427(87)90125-7
#' @references
#' Kohonen T. Self-Organizing Maps. Berlin, Heidelberg: Springer Berlin
#' Heidelberg (1995). doi:10.1007/978-3-642-97610-0
#' @references
#' Wehrens R, Kruisselbrink J. Flexible self-organizing maps in kohonen 3.0.
#' J Stat Softw (2018) 87:1–18. doi:10.18637/jss.v087.i07
#' @references
#' Venna J, Kaski S. Neighborhood preservation in nonlinear projection
#' methods: An experimental study. Lect Notes Comput Sci (including Subser Lect
#' Notes Artif Intell Lect Notes Bioinformatics) (2001) 2130:485–491.
#' doi:10.1007/3-540-44668-0_68
#'
#' @param data a numeric data frame, matrix or a `red_analysis` object. If a
#' `red_analysis` object is provided as the data argument, the observation
#' component/score table is subjected to clustering. For multi-layer SOM
#' analysis, the `data` argument is a list of such objects.
#' @param clust_assignment an optional data frame with two variables:
#' `observation` and `clust_id`, which defines the global cluster assignment.
#' If `NULL` (default), the global cluster assignment will be derived by
#' fitting the `clustering_fun` to the `data` - which can be a bit slower.
#' @param nfolds number of CV folds.
#' @param type type of the prediction algorithm: k-nearest neighbors
#' (propagation) or via the self-organizing map ('som', available only
#' for SOM and combined SOM clustering).
#' @param kNN number of the nearest neighbors.
#' @param simple_vote logical, should classical unweighted k-NN classification
#' be applied? If FALSE, distance-weighted k-NN is used with the provided kernel
#' function.
#' @param resolve_ties logical, should the ties be resolved at random? Applies
#' only to the simple unweighted voting algorithm.
#' @param kernel_fun kernel function transforming the distance into weight.
#' @param clustering_fun clustering function. Should return a
#' `clust_analysis` object.
#' @param kNN_data number of the nearest neighbors in the genuine data set
#' used for calculation of neighborhood preservation. See \code{\link{np}}
#' for details.
#' @param kNN_cluster number of the nearest neighbors of the given cluster used
#' for calculation of neighborhood preservation. See \code{\link{np}} for
#' details.
#' @param seed initial setting of the random number generator.
#' @param .parallel logical, should the CV be run in parallel?
#' @param ... extra arguments passed to the clustering_fun.
#'
#' @return a list containing the following elements:
#'
#' * the global \code{\link{clust_analysis}} object (`clust_analysis_object`)
#'
#' * kNN projection (prediction) results (`predictions`)
#'
#' * a data frame with the classification error, accuracy, fraction of
#' explained clustering variance and silhouette for the out-of-fold
#' predictions (`fold_stats`)
#'
#' * means and BCA's 95% confidence intervals for the classification error,
#' accuracy, fraction of explained variance and silhouette (`summary`)

  cv_cluster <- function(data,
                         clust_assignment = NULL,
                         nfolds = 5,
                         type = c('propagation', 'som'),
                         kNN = 5,
                         simple_vote = TRUE,
                         resolve_ties = FALSE,
                         kernel_fun = function(x) 1/x,
                         clustering_fun = clustTools::kcluster,
                         kNN_data = 5,
                         kNN_cluster = NULL,
                         seed = 1234,
                         .parallel = FALSE, ...) {

    ## entry control -------

    type <- match.arg(type[1], c('propagation', 'som'))

    nfolds <- as.integer(nfolds)

    kNN <- as.integer(kNN)

    stopifnot(is.logical(simple_vote))
    stopifnot(is.logical(resolve_ties))
    stopifnot(is.logical(.parallel))

    stopifnot(is.function(kernel_fun))

    global_clust <- NULL
    fold_clust <- NULL
    correct <- NULL

    start_time <- Sys.time()
    message(paste('CV: ', nfolds, 'folds'))
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))

    set.seed(seed = seed)

    ## data control --------

    if(is_red_analysis(data)) {

      data <- column_to_rownames(data$component_tbl,
                                 'observation')

    }

    if(is.data.frame(data) | is.matrix(data)) {

      check_numeric(data)

    } else if(is.list(data)) {

      purrr::walk(data, check_numeric)

    } else {

      stop('Unsupported data format', call. = FALSE)

    }

    ## defining the global classifier --------

    if(is.null(clust_assignment)) {

      glob_classif <- clustering_fun(data = data, ...)

      glob_assign <- glob_classif$clust_assignment[c('observation', 'clust_id')]

      glob_assign <- set_names(glob_assign,
                               c('observation', 'global_clust'))

    } else {

      stopifnot(is.data.frame(clust_assignment))

      if(any(!c('observation', 'clust_id') %in% names(clust_assignment))) {

        stop(paste("The 'clust_assignment' data frame has to contain",
                   "'observation' and 'clust_id' columns."),
             call. = FALSE)

      }

      glob_assign <- set_names(clust_assignment[c('observation', 'clust_id')],
                               c('observation', 'global_clust'))

      glob_classif <- NULL

    }

    ## creating the training set classifiers and out-of-fold predictions ------

    fold_set <- create_clust_folds(data, k = nfolds, seed = seed)

    if(.parallel) {

      future::plan('multisession')

      exports <- c('tidyverse',
                   'rlang',
                   'cluster',
                   'kohonen',
                   'dbscan',
                   'Rcpp',
                   'somKernels',
                   'clustTools',
                   'stringi',
                   'nomclust')

      train_classif <-
        furrr::future_map(fold_set$train ,
                          function(x) clustering_fun(data = x, ...),
                          .options = furrr::furrr_options(seed = TRUE,
                                                          packages = exports))

      test_obj <-
        furrr::future_pmap(list(object = train_classif,
                                newdata = fold_set$test),
                           purrr::safely(predict),
                           type = type,
                           kNN = kNN,
                           simple_vote = simple_vote,
                           resolve_ties = resolve_ties,
                           kernel_fun = kernel_fun,
                           .options = furrr::furrr_options(seed = TRUE,
                                                           packages = exports))

     future:: plan('sequential')

    } else {

      train_classif <- map(fold_set$train,
                           function(x) clustering_fun(data = x, ...))

      test_obj <- pmap(list(object = train_classif,
                            newdata = fold_set$test),
                       purrr::safely(predict),
                       type = type,
                       kNN = kNN,
                       simple_vote = simple_vote,
                       resolve_ties = resolve_ties,
                       kernel_fun = kernel_fun)

    }

    test_obj <- compact(map(test_obj, ~.x$result))

    ## classification errors ------

    test_preds <-
      map(test_obj,
          ~set_names(.x$clust_assignment[c('observation', 'clust_id')],
                     c('observation', 'fold_clust')))

    test_preds <- map(test_preds,
                      ~left_join(.x,
                                 glob_assign,
                                 by = 'observation'))

    correct <- NULL
    accuracy <- NULL
    error <- NULL
    fold <- NULL

    test_preds <-
      map2_dfr(test_preds,
               names(test_preds),
               ~mutate(.x,
                       correct = as.character(global_clust) == as.character(fold_clust),
                       fold = .y))

    test_error <- dplyr::summarise(test_preds,
                                   accuracy = mean(as.numeric(correct),
                                                   na.rm = TRUE),
                                   error = 1 - accuracy,
                                   .by = 'fold')

    ## variances, silhouettes and neighborhood preservation, out-of-fold  -------

    ## neighborhood preservation seems to be a bit tricky statistic to calculate
    ## with some distance metrics for a large dimension number, hence,
    ## I'm working with a safe solution

    test_stats <- map(test_obj,
                      purrr::safely(summary),
                      kNN_data = kNN_data,
                      kNN_cluster = kNN_cluster)

    test_stats <- map(test_stats, ~.x$result)

    test_stats <- map2_dfr(test_stats, names(test_stats),
                           ~mutate(.x, fold = .y))

    ## summary stats --------

    test_stats <- left_join(test_error,
                            test_stats,
                            by = 'fold')

    stat_names <-
      c('accuracy', 'error',
        'sil_width', 'frac_misclassified', 'frac_var', 'frac_np')

    test_stats <-
      select(test_stats,
             fold,
             dplyr::all_of(stat_names))

    oof_means <-
      map(select(test_stats,
                 dplyr::all_of(stat_names)),
          mean,
          na.rm = TRUE)

    bca_err <-
      map(select(test_stats,
                 dplyr::all_of(stat_names)),
          ~.x[!is.na(.x)])

    bca_err <- map(bca_err, coxed::bca, conf.level = 0.95)

    test_summary <-
      pmap(list(stat_names[stat_names %in% names(test_stats)],
                oof_means,
                bca_err),
           function(x, y, z) tibble(!!paste0(x, '_mean') := y,
                                    !!paste0(x, '_lower_ci') := z[1],
                                    !!paste0(x, '_upper_ci') := z[2]))

    ## output

    list(clust_analysis_object = glob_classif,
         predictions = test_preds,
         fold_stats = test_stats,
         summary = as_tibble(reduce(test_summary, cbind)))

  }

# S3 generics and methods --------

#' Cross-validate the clustering analysis object.
#'
#' @description
#' Checks the quality of a clustering solution by
#' cross-validation (CV) with k-nearest neighbors (kNN) out-of-fold predictions
#' or predictions made by a self-organizing map (SOM).
#' Stability of the clustering structure is measured by cluster assignment
#' classification error in the out-of-fold predictions as compared with the
#' genuine clustering structure. Explanatory value and cluster separation are
#' determined by clustering variance and silhouette statistics.
#'
#' @details
#' `cv()` is a S3 generic function.
#' By principle, cross-validation of a clustering structure is similar to
#' cross-validation of any machine learning multi-class classifier.
#' The training portion of a CV split is used to develop
#' of a cluster structure and the projection on the test portion is accomplished
#' by k-nearest neighbor (kNN) label propagation algorithm or derives the
#' cluster assignment from a trained SOM.
#' For implementation details, see: \code{\link{propagate}},
#' \code{\link{map_som}} and \code{\link{map_supersom}}.
#' The folds are generated with \code{\link[caret]{createFolds}}.
#' For `combi_analysis` objects, assignment of the observations to the CV folds
#' is done with the kNN algorithm for the 'top' assignment of the observations
#' to the clusters: nodes are ignored!
#' For `clust_analysis` and `combi_analysis` objects with multi-layered data
#' and clustering of U matrix, the SOM prediction method is the sole option.
#' Currently, it is not possible to cross-validate clustering analysis objects
#' generated with an user-provided dissimilarity matrices
#' (subclass `min_analysis` of `clust_analysis`).
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
#' @references
#' Rousseeuw PJ. Silhouettes: A graphical aid to the interpretation and
#' validation of cluster analysis. J Comput Appl Math (1987) 20:53–65.
#' doi:10.1016/0377-0427(87)90125-7
#' @references
#' Kohonen T. Self-Organizing Maps. Berlin, Heidelberg: Springer Berlin
#' Heidelberg (1995). doi:10.1007/978-3-642-97610-0
#' @references
#' Wehrens R, Kruisselbrink J. Flexible self-organizing maps in kohonen 3.0.
#' J Stat Softw (2018) 87:1–18. doi:10.18637/jss.v087.i07
#' @references
#' Venna J, Kaski S. Neighborhood preservation in nonlinear projection
#' methods: An experimental study. Lect Notes Comput Sci (including Subser Lect
#' Notes Artif Intell Lect Notes Bioinformatics) (2001) 2130:485–491.
#' doi:10.1007/3-540-44668-0_68
#'
#' @param x an object.
#' @param nfolds number of CV folds.
#' @param type type of the prediction algorithm: k-nearest neighbors
#' (propagation) or via the self-organizing map ('som', available only
#' for SOM and combined SOM clustering).
#' @param kNN number of the nearest neighbors.
#' @param simple_vote logical, should classical unweighted k-NN classification
#' be applied? If FALSE, distance-weighted k-NN is used with the provided kernel
#' function.
#' @param resolve_ties logical, should the ties be resolved at random? Applies
#' only to the simple unweighted voting algorithm.
#' @param kernel_fun kernel function transforming the distance into weight.
#' @param kNN_data number of the nearest neighbors in the genuine data set
#' used for calculation of neighborhood preservation. See \code{\link{np}}
#' for details.
#' @param kNN_cluster number of the nearest neighbors of the given cluster used
#' for calculation of neighborhood preservation. See \code{\link{np}} for
#' details.
#' @param seed initial setting of the random number generator.
#' @param .parallel logical, should the CV be run in parallel?
#' @param ... extra arguments, currently none.
#'
#' @return a list of class `cluster_cv` containing the following elements:
#'
#' * the global \code{\link{clust_analysis}} object (`clust_analysis_object`)
#'
#' * kNN projection (prediction) results (`predictions`)
#'
#' * a data frame with the classification error, accuracy, fraction of
#' explained clustering variance, silhouette and neighbor preservation for
#' the out-of-fold predictions (`fold_stats`)
#'
#' * means and BCA's 95% confidence intervals for the classification error,
#' accuracy, fraction of explained variance, silhouette and neighborhood
#' preservation (`summary`)
#'
#' Note the \code{\link{summary.cluster_cv}} and
#' \code{\link{extract.cluster_cv}} methods.
#'
#' @export

  cv <- function(x, ...) UseMethod('cv')

#' @rdname cv
#' @export

  cv.clust_analysis <- function(x,
                                nfolds = 5,
                                type = c('propagation', 'som'),
                                kNN = 5,
                                simple_vote = TRUE,
                                resolve_ties = FALSE,
                                kernel_fun = function(x) 1/x,
                                kNN_data = 5,
                                kNN_cluster = NULL,
                                seed = 1234,
                                .parallel = FALSE, ...) {

    ## entry control ---------

    stopifnot(is_clust_analysis(x))

    type <- match.arg(type[1], c('propagation', 'som'))

    nfolds <- as.integer(nfolds)

    kNN <- as.integer(kNN)
    kNN_data <- as.integer(kNN_data)

    if(is.null(kNN_cluster)) {

      if(x$clust_fun %in% c('som', 'supersom')) {

        kNN_cluster <- kNN_data

      } else {

        kNN_cluster <- 1

      }

    }

    stopifnot(is.logical(simple_vote))
    stopifnot(is.logical(resolve_ties))
    stopifnot(is.logical(.parallel))

    stopifnot(is.function(kernel_fun))

    if(type == 'som' & !x$clust_fun %in% c('som', 'supersom')) {

      warning('SOM predictions are available only for SOM clustering analyses.',
              call. = FALSE)

      return(NULL)

    }

    if(type == 'propagation' & x$clust_fun == 'supersom') {

      warning("kNN label propagation is not available for multi-layer SOM.",
              call. = FALSE)

      return(NULL)

    }

    if(x$clust_fun %in% c('prediction', 'supersom_prediction')) {

      warning("Cross-validation is not possible for predictions.",
              call. = FALSE)

      return(NULL)

    }

    ## common parameters --------

    if(x$clust_fun != 'supersom') {

      distance_names <- x$dist_method

    } else {

      distance_names <- x$clust_obj$dist.fcts

    }

    cmm_args <-
      list2(data = model.frame(x),
            nfolds = nfolds,
            type = type,
            kNN = kNN,
            simple_vote = simple_vote,
            resolve_ties = resolve_ties,
            kernel_fun = kernel_fun,
            distance_method = distance_names,
            kNN_data = kNN_data,
            kNN_cluster = kNN_cluster,
            seed = seed,
            .parallel = .parallel,
            k = nrow(ngroups(x)),
            hc_method = x$hc_method,
            clust_fun = x$clust_fun,
            lambdas = x$lambdas,
            xdim = x$grid$xdim,
            ydim = x$grid$ydim,
            topo = x$grid$topo,
            neighbourhood.fct = as.character(x$grid$neighbourhood.fct),
            toroidal = x$grid$toroidal,
            eps = x$eps,
            minPts = x$minPts)

    cmm_args <- c(cmm_args, x$dots)

    ## function calls -------

    if(x$clust_fun %in% c('hclust', 'som', 'supersom', 'htk')) {

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
                                    htk = htk_cluster,
                                    pam = kcluster,
                                    dbscan = dbscan_cluster,
                                    som = som_cluster,
                                    supersom = som_cluster))

    cluster_cv(eval(cv_call))

  }

#' @rdname cv
#' @export

  cv.min_analysis <- function(x, ...) {

    warning(paste('Cross-validation of clustering structures created',
                  'with user-provided dissimilarity matrices is not supported',
                  'at the moment.'),
            call. = FALSE)


    return(NULL)

  }

#' @rdname cv
#' @export

  cv.combi_analysis <- function(x,
                                nfolds = 5,
                                type = c('propagation', 'som'),
                                kNN = 5,
                                simple_vote = TRUE,
                                resolve_ties = FALSE,
                                kernel_fun = function(x) 1/x,
                                kNN_data = 5,
                                kNN_cluster = NULL,
                                seed = 1234,
                                .parallel = FALSE, ...) {

    ## entry control -------

    stopifnot(is_combi_analysis(x))

    type <- match.arg(type[1], c('propagation', 'som'))

    nfolds <- as.integer(nfolds)

    kNN <- as.integer(kNN)
    kNN_data <- as.integer(kNN_data)

    if(is.null(kNN_cluster)) kNN_cluster <- 1

    stopifnot(is.logical(simple_vote))
    stopifnot(is.logical(resolve_ties))
    stopifnot(is.logical(.parallel))

    stopifnot(is.function(kernel_fun))

    if(is_umatrix_analysis(x) & type == 'propagation') {

      warning(paste("Cross-validation of multi-layer SOM analyses is not",
                    "possible with the kNN label propagation algorithm.",
                    "Consider type = 'som' instead."),
              call. = FALSE)

      return(NULL)

    }

    ## cross-validation ---------

    if(is_umatrix_analysis(x)) {

      clustering_fun <- multi_cluster

      distance_som <- NULL

      distance_method <- x$clust_analyses$observation$clust_obj$dist.fcts

    } else {

      clustering_fun <- combi_cluster

      distance_som <- x$clust_analyses$observation$dist_method

      distance_method <- NULL

    }

    args <-
      list2(data = model.frame(x)$observation,
            nfolds = nfolds,
            type = type,
            kNN = 5,
            clustering_fun = clustering_fun,
            distance_som = distance_som,
            distance_method = distance_method,
            xdim = x$clust_analyses$observation$grid$xdim,
            ydim = x$clust_analyses$observation$grid$ydim,
            topo = x$clust_analyses$observation$grid$topo,
            neighbourhood.fct = as.character(x$clust_analyses$observation$grid$neighbourhood.fct),
            toroidal = x$clust_analyses$observation$grid$toroidal,
            rlen = nrow(x$clust_analyses$observation$clust_obj$changes),
            node_clust_fun = switch(x$clust_analyses$node$clust_fun,
                                    hclust = hcluster,
                                    kmeans = kcluster,
                                    htk = htk_cluster,
                                    pam = kcluster,
                                    dbscan = dbscan_cluster,
                                    som = som_cluster),
            kernel_fun = kernel_fun,
            simple_vote = simple_vote,
            resolve_ties = resolve_ties,
            kNN_data = kNN_data,
            kNN_cluster = kNN_cluster,
            seed = seed,
            .parallel = .parallel)

    args <- c(args, x$dots)

    args <- compact(args)

    cv_call <- call2('cv_cluster',
                     !!!args)

    cluster_cv(eval(cv_call))

  }

# END -----
