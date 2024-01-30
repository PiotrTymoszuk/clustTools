# Functions for regularized clustering

# Regularized K-means clustering -----

#' Regularized hard threshold KMEANS clustering.
#'
#' @description
#' Implements the hard threshold KMEANS algorithm proposed by Raymaekers
#' and Zamar and provided by the `clusterHD` package.
#'
#' @details
#' The algorithm offers an interesting approach to clustering of
#'  multi-dimensional data, especially those containing variables
#' of little relevance for the clustering analysis (e.g. as investigated by
#' permutation importance with \code{\link{impact}}). For details, please refer
#' to the genuine R function \code{\link[clusterHD]{HTKmeans}} and the seminal
#' paper. The function works with the squared Euclidean distance metric and
#' accepts only a numeric data frame as the input data.
#' There are two crucial parameters to be provided by the user, the number of
#' centers/clusters `k`, which can be determined by methods such as peak mean
#' silhouette width, and the regularization argument `lambdas`, whose value can
#' be found by tuning (e.g. comparing silhouette widths or clustering variances
#' for various lambda values). If `lambdas` is set to `NULL` or provided as
#' a numeric vector, the best lambda value is found with the
#' \code{\link[clusterHD]{getLambda}}. If a single value is provided, it
#' will be used for clustering.
#' Tuning of the lambda parameter using explained clustering variance,
#' silhouette widths and neighbor preservation statistic is facilitated by the
#' `tune_htk()` function.
#'
#' @references
#' Raymaekers J, Zamar RH. Regularized K-means Through Hard-Thresholding.
#' J Mach Learn Res (2022) 23:1–48.
#' Available at: http://jmlr.org/papers/v23/21-0052.html
#'
#' @param data a numeric data frame with observations in the rows and
#' variables in the columns.
#' @param k number of centers (clusters).
#' @param lambdas a numeric vector of the regularization parameter. See Details.
#' @param select_stat statistic used for selection of the optimal lambda value.
#' Ignored if `lambdas` is a single numeric value. For `htk_cluster()` they are
#' 'AIC' (Akaike Information Criterion) or 'BIC' (Bayesian Information
#' Criterion). For `tune_htk()` they are silhouette width ('silhouette'),
#' fraction of observations with negative silhouette values
#' ('misclassification'), fraction of explained clustering variance ('variance'),
#' or neighborhood preservation ('np').
#' @param seed root of the random number generator.
#' @param type type of the tuning procedure. When set to 'train' (default),
#' cluster structure quality statistics are computed for the entire data set.
#' When set to 'cv', cross-validated statistics for subsequent lambda values
#' are calculated.
#' @param nfolds number of CV folds.
#' @param kNN number of the nearest neighbors used by the cluster assignment
#' classifier in the cross-validation folds. Ignored if `type = 'train'`.
#' @param simple_vote logical, should classical unweighted k-NN classification
#' be applied? If FALSE, distance-weighted k-NN is used with the provided kernel
#' function. Ignored if `type = 'train'`.
#' @param resolve_ties logical, should the ties be resolved at random? Applies
#' only to the simple unweighted voting algorithm. Ignored if `type = 'train'`.
#' @param kernel_fun kernel function transforming the distance into weight.
#' Ignored if `type = 'train'`.
#' @param kNN_data number of the nearest neighbors in the genuine data set
#' used for calculation of neighborhood preservation. See \code{\link{np}}
#' for details.
#' @param kNN_cluster number of the nearest neighbors of the given cluster used
#' for calculation of neighborhood preservation. See \code{\link{np}} for
#' details.
#' @param .parallel logical, shoudl the analysis be run in parallel?
#' @param ... extra arguments provided to \code{\link[clusterHD]{HTKmeans}}.
#'
#' @return `htk_cluster()` returns an object of the
#' class \code{\link{clust_analysis}}.
#'
#' @export

  htk_cluster <- function(data,
                          k = 2,
                          lambdas = NULL,
                          select_stat = c('AIC', 'BIC'),
                          seed = 1234, ...) {

    ## entry control -------

    check_numeric(data)

    stopifnot(is.numeric(k))

    k <- as.integer(k[1])

    if(!is.null(lambdas)) {

      if(!is.numeric(lambdas)) {

        stop(paste("'lambdas' has to be a numeric vector, NULL or a",
                   "single numeric value."),
             call. = FALSE)

      }

    }

    select_stat <- match.arg(select_stat[1], c('AIC', 'BIC'))

    ## check for the data ------

    check_numeric(data)

    mat_data <- as.matrix(data)

    if(is.null(rownames(mat_data))) {

      obs_id <- 1:nrow(mat_data)

    } else {

      obs_id <- rownames(mat_data)

    }

    model_frame <- enexpr(data)

    ## clustering ------

    ## handling distance method requests, to make the function compatible
    ## with combined SOM procedures

    extra_args <- list2(...)

    extra_args$distance_method <- NULL

    extra_args <- compact(extra_args)

    fun <- clusterHD::HTKmeans

    res_call <- call2('fun',
                      X = mat_data,
                      k = k,
                      lambdas = lambdas,
                      !!!extra_args)

    res <- eval(res_call)

    ## lambdas as a single number --------

    if(length(lambdas) == 1) {

      set.seed(seed)

      clust_ass <-
        tibble(observation = obs_id,
               clust_id = factor(res$HTKmeans.out[[1]]$cluster))

      ## determining the distance of centers from variable means
      ## and active variables

      clust_vars <- colnames(mat_data)

      center_mat <- res$HTKmeans.out[[1]]$centers

      colnames(center_mat) <- clust_vars

      center_distances <- colSums(center_mat^2)

      active_variables <- clust_vars[center_distances != 0]

      return(clust_analysis(list(data = quo(!!model_frame),
                                 dist_mtx = calculate_dist(data, 'squared_euclidean'),
                                 dist_method = 'squared_euclidean',
                                 lambdas = lambdas,
                                 center_distances = center_distances,
                                 active_variables = active_variables,
                                 clust_fun = 'htk',
                                 clust_obj = res,
                                 clust_assignment = clust_ass,
                                 dots = list2(...))))

    }

    ## lambda as NULL or a numeric vector ------

    opt_lambda <- clusterHD::getLambda(res, type = select_stat)

    return(htk_cluster(data = data,
                       k = k,
                       lambdas = opt_lambda,
                       seed = seed, ...))

  }

#' @rdname htk_cluster
#' @export

  tune_htk <- function(data,
                       k = 2,
                       lambdas = NULL,
                       select_stat = c('silhouette',
                                       'misclassification',
                                       'variance',
                                       'np'),
                       type = c('train', 'cv'),
                       nfolds = 5,
                       kNN = 5,
                       simple_vote = TRUE,
                       resolve_ties = FALSE,
                       kernel_fun = function(x) 1/x,
                       kNN_data = 5,
                       kNN_cluster = NULL,
                       seed = 1234,
                       .parallel = FALSE, ...) {

    ## input check ------

    start_time <- Sys.time()

    check_numeric(data)

    clust_features <- names(data)

    stopifnot(is.numeric(k))

    k <- as.integer(k[1])

    type <- match.arg(type[1], c('train', 'cv'))

    stopifnot(is.logical(.parallel))

    select_stat <- match.arg(select_stat[1],
                             c('silhouette',
                               'misclassification',
                               'variance',
                               'np'))

    ## the CV arguments are checked by a downstream function

    ## the lambda vector is stripped from the output of the HTKmeans()
    ## function

    if(is.null(lambdas)) {

      res <- clusterHD::HTKmeans(X = as.matrix(data),
                                 k = k,
                                 lambdas = lambdas, ...)

      lambdas <- res$lambdas

    }

    stopifnot(is.numeric(lambdas))

    ## benchmarking and parallel backend ---------

    message(paste('Tuning for', length(lambdas), 'lambda values'))

    if(.parallel) future::plan('multisession')

    on.exit(message(paste('Elapsed:', Sys.time() - start_time)),
            add = TRUE,
            after = FALSE)

    on.exit(future::plan('sequential'))

    ## tuning objects --------

    extra_args <- list2(...)

    tune_calls <-
      map(lambdas,
          ~call2('htk_cluster',
                 data = data,
                 k = k,
                 lambdas = .x,
                 select_stat = 'AIC',
                 seed = seed,
                 !!!compact(extra_args)))

    tune_objects <-
      furrr::future_map(tune_calls, eval,
                        .options = furrr::furrr_options(seed = seed))

    tune_objects <- set_names(tune_objects,
                              paste0('obj_', 1:length(tune_objects)))

    lambda <- NULL
    object_id <- NULL

    lambda_assign <- tibble(lambda = lambdas,
                            object_id = names(tune_objects))

    ## training data stats ------

    if(type == 'train') {

      lambda <- NULL

      ## working with safely(), because clustering for some lambda values
      ## can fail

      stats <-
        furrr::future_map2(tune_objects,
                           names(tune_objects),
                           ~purrr::safely(tibble)(mutate(summary(.x,
                                                                 kNN_data = kNN_data,
                                                                 kNN_cluster = kNN_cluster)),
                                                  object_id = .y),
                           .options = furrr::furrr_options(seed = TRUE))

      stats <- map_dfr(stats, ~.x$result)

      stats <- left_join(stats, lambda_assign, by = 'object_id')

    }

    ## CV stats -------

    if(type == 'cv') {

      if(.parallel) {

        if(nfolds > length(lambdas)) {

          cv_objects <-
            map(tune_objects,
                ~purrr::safely(cv)(.x,
                                   type = 'propagation',
                                   nfolds = nfolds,
                                   kNN = kNN,
                                   simple_vote = simple_vote,
                                   resolve_ties = resolve_ties,
                                   kernel_fun = kernel_fun,
                                   kNN_data = kNN_data,
                                   kNN_cluster = kNN_cluster,
                                   seed = seed,
                                   .parallel = TRUE))

        } else {

          cv_objects <-
            furrr::future_map(tune_objects,
                              ~purrr::safely(cv)(.x,
                                                 type = 'propagation',
                                                 nfolds = nfolds,
                                                 kNN = kNN,
                                                 simple_vote = simple_vote,
                                                 resolve_ties = resolve_ties,
                                                 kernel_fun = kernel_fun,
                                                 kNN_data = kNN_data,
                                                 kNN_cluster = kNN_cluster,
                                                 seed = seed,
                                                 .parallel = FALSE),
                              .options = furrr::furrr_options(seed = TRUE))

        }

      } else {

        cv_objects <-
          map(tune_objects,
              ~purrr::safely(cv)(.x,
                                 type = 'propagation',
                                 nfolds = nfolds,
                                 kNN = kNN,
                                 simple_vote = simple_vote,
                                 resolve_ties = resolve_ties,
                                 kernel_fun = kernel_fun,
                                 kNN_data = kNN_data,
                                 kNN_cluster = kNN_cluster,
                                 seed = seed,
                                 .parallel = FALSE))


      }


      correct <- map_lgl(cv_objects, ~is.null(.x$error))

      cv_objects <- compact(map(cv_objects, ~.x$result))

      stats <- map(cv_objects, summary)

      stats <- map2_dfr(stats, names(stats),
                        ~mutate(.x, object_id = .y))

      stats <- left_join(stats, lambda_assign, by = 'object_id')

      frac_var_mean <- NULL
      sil_width_mean <- NULL
      frac_np_mean <- NULL
      frac_misclassified_mean <- NULL

      sil_width <- NULL
      frac_misclassified <- NULL
      frac_var <- NULL
      frac_np <- NULL

      stats <- dplyr::transmute(stats,
                                object_id = object_id,
                                lambda = lambda,
                                sil_width = sil_width_mean,
                                frac_misclassified = frac_misclassified_mean,
                                frac_var = frac_var_mean,
                                frac_np = frac_np_mean)

      stats <- mutate(stats, lambda = lambdas[correct])

    }

    stats <- dplyr::relocate(stats, lambda)

    ## centers and identification of active variables -------

    ## active variables are those for which all center coordinates
    ## are non zero.
    ## I'm extracting the squared Euclidean distance to the zero point (mean)

    center_distances <- map(tune_objects, ~.x$center_distances)

    n_active_vars <- map_dbl(tune_objects, ~length(.x$active_variables))

    center_distances <- as_tibble(do.call('rbind', center_distances))

    center_distances <- mutate(center_distances,
                               lambda = lambdas,
                               n_active_vars = n_active_vars)

    stats <- left_join(stats,
                       center_distances[c('lambda', 'n_active_vars')],
                       by = 'lambda')

    ## best tune and the analysis object -------

    select_var <- switch(select_stat,
                         silhouette = 'sil_width',
                         misclassification = 'frac_misclassified',
                         variance = 'frac_var',
                         np = 'frac_np')

    select_fun <-
      switch(select_stat,
             silhouette = function(x) max(x, na.rm = TRUE),
             misclassification = function(x) min(x, na.rm = TRUE),
             variance = function(x) max(x, na.rm = TRUE),
             np = function(x) max(x, na.rm = TRUE))

    select_lab <- switch(select_stat,
                         silhouette = 'max',
                         misclassification = 'min',
                         variance = 'max',
                         np = 'max')

    ## in case there are more than one optimal tunes selected by the numeric
    ## statistic criterion, the scenario with the lowest number of active
    ## variables is selected

    best_stats <-
      filter(stats,
             .data[[select_var]] == select_fun(.data[[select_var]]))

    best_stats <-
      filter(best_stats,
             n_active_vars == min(n_active_vars))

    best_stats <- best_stats[1, ]

    ## the output object ---------

    tuner(list(analysis = tune_objects[[best_stats$object_id[[1]]]],
               stats = stats,
               center_distances = center_distances,
               fun = 'tune_htk',
               dataset = type,
               type = 'development',
               clust_vars = clust_features,
               tune_params = 'lambda',
               tune_criteria = tibble(!!select_var := select_lab,
                                      n_active_vars = 'min'),
               best_tune = tibble(lambda = best_stats$lambda[[1]])))

  }

# END -------
