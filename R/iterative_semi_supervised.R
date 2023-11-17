# Iterative kNN label propagation algorithm used for semi-supervised clustering

# S3 interface -------

#' Iterative k-nearest neighbor label propagation algorithm.
#'
#' @description
#' Prediction of cluster assignment by the k-nearest neighbor label propagation
#' algorithm with an automated choice of the k value based on a loss function.
#'
#' @details
#' The function finds the optimal value of k for the k-nearest neighbor
#' classifier by iteratively checking quality of the cluster assignment. The
#' quality check is accomplished by one of the loss functions:
#' silhouette width (`select_stat = 'silhouette'`, default),
#' percentage of observations with negative silhouette widths
#' ('misclassification'),
#' fraction of explained clustering variance (i.e. ratio of the between-cluster
#' sum of squares to the total sum of squares, `select_stat = 'variance'`),
#' or neighbor preservation (`select_stat = 'np'`).
#' The `prediter()` function is a S3 generic function.
#' The function works only for clustering analyses and combined SOM - clustering
#' analyses with the data provided as data frames but not as distance matrices.
#'
#' @references
#' Leng M, Wang J, Cheng J, Zhou H, Chen X. Adaptive semi-supervised
#' clustering algorithm with label propagation. J Softw Eng (2014) 8:14–22.
#' doi:10.3923/jse.2014.14.22
#'
#' @return an object of the \code{\link{tuner}} class with the `plot()` and `summary()`
#' methods.
#'
#' @param x a `clust_analysis` or `combi_analysis` object.
#' @param newdata a data frame or a matrix with the new data. If `NULL`,
#' the training data is used for fitting the clustering structure.
#' @param select_stat a name of the loss function defining the quality measure
#' of the prediction. For details, see: Details.
#' @param max_k the maximal number of the nearest neighbors to be tested.
#' @param kNN_data number of the nearest neighbors in the dataset, used for
#' determination of neighborhood preservation statistic. See: \code{\link{np}}
#' for details.
#' @param kNN_cluster number of the nearest neighbors in the cluster, used for
#' determination of neighborhood preservation statistic. See: \code{\link{np}}
#' for details.
#' @param .parallel logical, should the analysis be run in parallel?
#' @param ... extra arguments passed to \code{\link{predict.clust_analysis}} or
#' \code{\link{predict.combi_analysis}} such as kernel weighting. Note that you
#' cannot specify the `kNN` and `type` arguments.
#'
#' @export

  prediter <- function(x, ...) UseMethod('prediter')

#' @rdname prediter
#' @export

  prediter.clust_analysis <- function(x,
                                      newdata = NULL,
                                      select_stat = c('silhouette',
                                                      'misclassification',
                                                      'variance',
                                                      'np'),
                                      max_k = 20,
                                      kNN_data = 5,
                                      kNN_cluster = NULL,
                                      .parallel = FALSE, ...) {

    ## entry control -------

    start_time <- Sys.time()

    stopifnot(is_clust_analysis(x) | is_combi_analysis(x))

    select_stat <- match.arg(select_stat[1],
                             c('silhouette',
                               'misclassification',
                               'variance',
                               'np'))

    stopifnot(is.numeric(max_k))
    stopifnot(is.logical(.parallel))

    n_observations <- nobs(x)$observations

    if(max_k > n_observations) {

      stop("'max_k' has to be lower than the observation number.",
           call. = FALSE)

    }

    max_k <- as.integer(max_k)

    k_vec <- 1:max_k

    ## Benchmarking ---------

    message(paste('Finding optimal number of nearest neighbors for',
                  length(k_vec), 'k values'))

    on.exit(message(paste('Elapsed:', Sys.time() - start_time)),
            add = TRUE,
            after = FALSE)

    on.exit(future::plan('sequential'))

    ## data check: essentially done by `predict()` ------

    if(is.null(newdata)) newdata <- model.frame(x)

    check_numeric(newdata)

    ## predictions and quality stats ---------

    if(.parallel) {

      future::plan('multisession')

      preds <-
        furrr::future_pmap(list(kNN = k_vec),
                           propagate,
                           object = x,
                           newdata = newdata,
                           ...,
                           .options = furrr::furrr_options(seed = TRUE))

      stats <-
        furrr::future_map_dfr(preds,
                              summary,
                              kNN_data = kNN_data,
                              kNN_cluster = kNN_cluster,
                              .options = furrr::furrr_options(seed = TRUE))

    } else {

      preds <- pmap(list(kNN = k_vec),
                    propagate,
                    object = x,
                    newdata = newdata, ...)

      stats <- map_dfr(preds,
                       summary,
                       kNN_data = kNN_data,
                       kNN_cluster = kNN_cluster)

    }

    kNN <- NULL

    stats <- mutate(stats, kNN = k_vec)

    stats <- dplyr::relocate(stats, kNN)

    ## choice of the optimal k and output -------

    select_var <- switch(select_stat,
                         silhouette = 'sil_width',
                         misclassification = 'frac_misclassified',
                         variance = 'frac_var',
                         np = 'frac_np')

    select_fun <- switch(select_stat,
                         silhouette = function(x) max(x, na.rm = TRUE),
                         misclassification = function(x) min(x, na.rm = TRUE),
                         variance = function(x) max(x, na.rm = TRUE),
                         np = function(x) max(x, na.rm = TRUE))

    select_lab <- switch(select_stat,
                         silhouette = 'max',
                         misclassification = 'min',
                         variance = 'max',
                         np = 'max')

    opt_k <-
      filter(stats,
             .data[[select_var]] == select_fun(.data[[select_var]]))

    opt_k <- opt_k$kNN[1]

    tuner(list(analysis = preds[[opt_k]],
               stats = stats,
               fun = 'prediter',
               dataset = 'train',
               type = 'prediction',
               clust_vars = names(model.frame(x)),
               tune_params = 'kNN',
               tune_criteria = tibble(!!select_var := select_lab),
               best_tune = tibble(kNN = opt_k)))

  }

#' @rdname prediter
#' @export

  prediter.min_analysis <- function(x, ...) {

    warning(paste('Iterative predictions are not available for clustering',
                  'analyses done with distance matrices.'),
            call. = FALSE)

    return(NULL)

  }

#' @rdname prediter
#' @export

  prediter.umatrix_analysis <- function(x, ...) {

    warning(paste('Iterative predictions are not available for clustering',
                  'analyses done with multiple data layers.'),
            call. = FALSE)

    return(NULL)

  }

#' @rdname prediter
#' @export

  prediter.combi_analysis <- function(x,
                                      newdata = NULL,
                                      select_stat = c('silhouette',
                                                      'misclassification',
                                                      'variance',
                                                      'np'),
                                      max_k = 20,
                                      kNN_data = 5,
                                      kNN_cluster = NULL,
                                      .parallel = FALSE, ...) {

    ## entry control -------

    prediter.clust_analysis(x,
                            newdata = newdata,
                            select_stat = select_stat,
                            max_k = max_k,
                            .parallel = .parallel, ...)

  }

# END -----
