# Contains the class constructor functions.

# Dimensionality reduction analysis class ------

#' Reduction analysis class object.
#'
#' @description
#' Constructs a `red_analysis` class object given a list storing
#' results of a dimensionality reduction analysis.
#'
#' @details
#' A named list with the following elements is required as the `x` argument:
#'
#' * `red_obj` with the analysis output,
#'
#' * `red_fun` name of the reduction function,
#'
#' * `dist_method` name of the distance metric,
#'
#' * `component_tbl` a data frame with component/score values for the observations,
#'
#' * `loadings` a data frame with variable loadings, relevant e.g. for PCA or
#' factor analysis,
#'
#' `data` a quosure calling the original data set.
#'
#' If the `component_tbl` data frame contains the cluster assignment information
#' an object of the `clust_red` sub-class is created. This sub-class inherits
#' almost all methods from the superclass `red_analysis`. The only difference
#' concerns the `plot()` method, which, by default generates layout (score)
#' scatter plots with the cluster assignment coded by the point color.
#'
#' @param x a named list, see Details.
#'
#' @return a `red_analysis` object with the elements listed in Details.
#'
#' @export

  red_analysis <- function(x) {

    ## entry control

    if(!is.list(x)) {

      stop('Only named lists can be converted to a clust_analysis objects',
           call. = FALSE)

    }

    if(is.null(names(x))) {

      stop('Only named lists can be converted to a clust_analysis objects',
           call. = FALSE)

    }

    if(!all(c('red_obj',
              'red_fun',
              'dist_method',
              'component_tbl',
              'loadings',
              'data') %in% names(x))) {

      stop(paste('The input list needs to provide red_obj,',
                 'red_fun, component_tbl, loadings and data elements'),
           call. = FALSE)

    }

    stopifnot(rlang::is_quosure(x$data))

    stopifnot(x$dist.method %in% c(get_kernel_info(), 'weighted_som', 'custom'))

    ## output

    x <- structure(x, class = 'red_analysis')

    if(is.data.frame(x$component_tbl)) {

      if('clust_id' %in% names(x$component_tbl)) {

        x <- structure(x, class = c('clust_red', class(x)))

      }

    }

    x

  }

# Simple clustering class -----

#' Clustering analysis object.
#'
#' @description
#' Constructs a `clust_analysis` class object given a list with
#' results of a clustering analysis.
#'
#' @details A named list with the following elements is required
#' as the `x` argument:
#'
#' * `data`: a quosure calling the original data set,
#'
#' * `dist_mtx`: a numeric matrix with the distances between the observations,
#'
#' * `dist_method`: name of the distance statistic,
#'
#' * `clust_obj`: the output object of the clustering analysis,
#'
#' * `clust_fun`: the name of the clustering function or prediction,
#'
#' * `clust_assignment`: a data frame with the cluster assignment of the
#' observations. It has to contain the variables `observation` and `clust_id`,
#'
#' * `dots`: additional arguments passed to the clustering function.
#'
#' The `clust_analysis` object can be created for clustering solutions based on
#' data frames or matrices or user-provided distance matrices (clustering
#' functions working with dissimilarity objects of the `dist` class).
#' In the later case an instance of the subclass `min_analysis` is returned.
#' As such, the `min_analysis` class inherits most of the methods for
#' `clust_analysis` objects. However, some methods requiring source tabular data
#' like reduction analysis with the genuine data frame, variable importance
#' or heterologous cross-distances will not be available.
#' Semi-supervised clustering and cross-validation will be implemented for
#' the `min_analysis` class in the future.
#'
#' @param x a named list, see Details.
#'
#' @return a `clust_analysis` object with the elements listed in Details.
#'
#' @export

  clust_analysis <- function(x) {

    ## entry control -------

    if(!is.list(x)) {

      stop('Only named lists can be converted to a clust_analysis objects',
           call. = FALSE)

    }

    if(is.null(names(x))) {

      stop('Only named lists can be converted to a clust_analysis objects',
           call. = FALSE)

    }

    if(!all(c('data',
              'dist_mtx',
              'dist_method',
              'clust_obj',
              'clust_fun',
              'clust_assignment',
              'dots') %in% names(x))) {

      stop(paste('The input list needs to provide data, dist_mtx,',
                 'dist_method, clust_fun, clust_obj, clust_assignment',
                 'and dots elements'),
           call. = FALSE)

    }

    stopifnot(rlang::is_quosure(x$data))

    stopifnot(is.matrix(x$dist_mtx))

    stopifnot(x$dist.method %in% c(get_kernel_info(), 'weighted_som', 'custom'))

    stopifnot(x$clust_fun %in% c('hclust',
                                 'kmeans',
                                 'htk',
                                 'dbscan',
                                 'som',
                                 'pam',
                                 'prediction',
                                 'supersom',
                                 'supersom_prediction'))

    stopifnot(is.data.frame(x$clust_assignment))

    stopifnot(all(c('observation', 'clust_id') %in% names(x$clust_assignment)))

    ## output ------

    dist_check <- inherits(eval_tidy(x$data), 'dist')

    if(!dist_check) {

      x <- structure(x, class = 'clust_analysis')

    } else {

      x$dist_method <- 'custom'

      x <- structure(x, class = c('min_analysis', 'clust_analysis'))

    }

    if(!is.null(x$active_variables)) {

      x <- structure(x, class = c('hard_threshold', class(x)))

    }

    x

  }

# Combined reduction analysis clustering class -----

#' Combined SOM - clustering analysis.
#'
#' @description
#' Constructs a `combi_anlysis()` class object given a list with
#' results of reduction analysis or self-organizing map and clustering analysis.
#'
#' @details
#' A named list with the following elements is required as the `x` argument:
#'
#' * `clust_analyses`: a list of `red_analysis` or `clust_analysis` objects,
#'
#' * `clust_assignment`: a data frame with the cluster assignment with the
#' `observation` and `clust_id` variables.
#'
#' For combined solutions involving unsupervised clustering of the SOM U matrix,
#' the function returns an object of subclass
#' `umatrix_analysis`, which inherits virtually all methods from the
#' `combi_analysis` class.
#'
#' @param x a named list, see Details.
#'
#' @return a `combi_analysis` object with the elements specified in Details.

  combi_analysis <- function(x) {

    ## entry control

    if(!is.list(x)) stop('Only named lists can be converted to a clust_analysis objects',
                         call. = FALSE)

    if(is.null(names(x))) stop('Only named lists can be converted to a clust_analysis objects',
                               call. = FALSE)

    if(!all(c('clust_analyses', 'clust_assignment') %in% names(x))) {

      stop('The input list needs to provide clust_analyses and clust_assignment elements',
           call. = FALSE)

    }

    stopifnot(is.data.frame(x$clust_assignment))

    stopifnot(all(c('observation',
                    'clust_id') %in% names(x$clust_assignment)))

    ## output

    min_check <- map_lgl(x$clust_analyses, is_min_analysis)

    if(any(min_check)) {

      return(structure(x, class = c('umatrix_analysis', 'combi_analysis')))

    } else {

      return(structure(x, class = 'combi_analysis'))

    }

  }

# Importance table class -----

#' Construct an importance object.
#'
#' @description
#' Constructs an object of class `importance` on the top of
#' a tibble with the clustering variable importance testing results.
#'
#' @param x a data frame with the following columns:
#'
#' * `variable` with the names of variables
#'
#' * `total_wss` with total within-cluster sum of squares
#'
#' * `total_ss` with total sum of squares
#'
#' * `between_ss` with between-cluster sum of squares
#'
#' * `frac_var` with fraction of explained clustering variance
#'
#' * `frac_diff` with difference in fraction of explained clustering variance
#' between the genuine clustering analysis and the clustering analysis done with
#' the variable of interest modified. This variable stores the actual variable
#' importance metric.
#'
#' @return a tibble of the `importance` class with the variables listed
#' in Details.

  importance <- function(x) {

    if(!is.data.frame(x)) {

      stop('Please provide a data frame as an input.',
           call. = FALSE)

    }

    if(!all(c('total_wss',
              'total_ss',
              'between_ss',
              'frac_var',
              'variable',
              'frac_diff') %in% names(x))) {

      stop(paste('The input data frame needs to have total_wss, total_ss,',
                 'between_ss, frac_var, variable and frac_diff columns.'),
           call. = FALSE)

    }

    ## output

    x <- as_tibble(x)

    structure(x, class = c('importance', class(x)))

  }

# Cross-distance class -------

#' Construct a cross_dist object.
#'
#' @description
#' The `cross_dist` class objects a created on the top of a named list
#' of matrices of cross-distances between clusters.
#' The list elements have to be named after the compared clusters following
#' the 'clust1 vs clust2' scheme.
#'
#' @param x a named list of cross-distance matrices.#'
#' @param type type of the cross-distances:
#' `homologous` for comparison of within the same clustering structure or
#' `heterologous` for comparison of two clustering structures.
#' @param method name of the distance metric.
#' @param x_levels order of the clusters of the
#' first cluster/combi analysis object.
#' @param y_levels order of the clusters of the
#' second cluster/combi analysis object.
#' @param ... extra arguments, currently none defined.
#'
#' @return an object of class `cross-distance` being a list of
#' cross-distance matrices.
#' Information on the comparison type and distance
#' metric are stored as the `type` and `dist_method` attributes.
#' Information on the cluster order is stored as the
#' `x_levels` and `y_levels` attributes.

  cross_dist <- function(x,
                         type = c('homologous', 'heterologous'),
                         method = 'euclidean',
                         x_levels = NULL,
                         y_levels = NULL, ...) {

    ## entry control -------

    err_txt <- "The 'x' argument has to be a named list of numeric matrices."

    if(!is.list(x)) stop(err_txt, call. = FALSE)

    mtx_check <- map_lgl(x, is.matrix)

    if(!all(mtx_check)) stop(err_txt, call. = FALSE)

    num_check <- map_lgl(x, is.numeric)

    if(!all(num_check)) stop(err_txt, call. = FALSE)

    type <- match.arg(type[1], c('homologous', 'heterologous'))

    if(!method %in% c(get_kernel_info(), 'weighted_som', 'custom')) {

      stop(paste("Unsupported distance measure.",
                 "Please consult 'get_kernel_info()."),
           call. = FALSE)

    }

    ## object -------

    x <- structure(x, class = c('cross_dist', class(x)))

    attr(x, 'type') <- type
    attr(x, 'dist_method') <- method
    attr(x, 'x_levels') <- x_levels
    attr(x, 'y_levels') <- y_levels

    x

  }

# Silhouette ------

#' Generate a sil_extra object.
#'
#' @description
#' Extends the \code{\link[cluster]{silhouette}} object by
#' cluster order and names.
#'
#' @details
#' The `sil_extra` class has \code{\link{summary.sil_extra}}
#' and \code{\link{plot.sil_extra}}
#' methods compatible with the tidyverse environment.
#'
#' @param x an object of the \code{\link[cluster]{silhouette}} class.
#' @param assignment an data frame with the `clust_id` and `observation`
#' columns defining the cluster assignment, e.g. obtained by the
#' \code{\link{extract}} function applied to a `clust_analysis` or
#' `combi_analysis` object.
#'
#' @return an object of the `sil_extra` class. Technically, a tibble with
#' the observation ID, cluster name, neighbor cluster name and silhouette width.

  sil_extra <- function(x, assignment) {

    ## entry control -----

    if(!inherits(x, 'silhouette')) {

      stop("'x' must be a 'silhouette' class object.",
           call. = FALSE)

    }

    err_txt <- paste("'assignment' needs to be a data frame with the
                     columns 'observation' and 'clust_id'.")

    if(!is.data.frame(assignment)) stop(err_txt, call. = FALSE)

    if(!all(c('observation', 'clust_id') %in% names(assignment))) {

      stop(err_txt, call. = FALSE)

    }

    if(nrow(x) != nrow(assignment)) {

      stop("The number of rows in 'x' and 'assignment' must be equal.",
           call. = FALSE)

    }

    ## the output tibble ------

    clust_levs <- levels(assignment$clust_id)

    observation <- NULL
    clust_id <- NULL
    neighbor_id <- NULL
    sil_width <- NULL

    sil_extra <-
      tibble(observation = assignment$observation,
             clust_id = clust_levs[x[, 'cluster']],
             neighbor_id = clust_levs[x[, 'neighbor']],
             sil_width = x[, 'sil_width'])

    sil_extra <- mutate(sil_extra,
                        clust_id = factor(clust_id, clust_levs),
                        neighbor_id = factor(neighbor_id, clust_levs))

    structure(sil_extra,
              class = c('sil_extra', class(sil_extra)))

  }

# Cross-validation -------

#' Cross validation results.
#'
#' @description
#' Creates an object of the `cluster_cv` class on a top of a list.
#'
#' @details
#' The `x` argument has to be a list with the following elements:
#'
#' * `clust_analysis_object` which stores a `clust_analysis` or
#' `combi_analysis` object with the global cluster assignment
#'
#' * `predictions` with a data frame storing out-of-fold predictions with the
#' variables `observation`, `fold_clust` (out-of-fold cluster assignment),
#' `global_clust` (cluster assignment in the global clustering structure),
#' `correct` (a logical indicating the out-of-fold - global assignment
#' concordance) and `fold`(fold ID).
#'
#' * `fold_stats`, a data frame which stores fold means of
#' accuracy (`corr_rate`), classification error (`err_rate`), fraction of
#' explained clustering variance (`frac_var`) and silhouette width (`sil_width`)
#'
#' * `summary` storing the global means and BCA 95% confidence intervals listed
#' in `fold_stats`
#'
#' @param x a named list with elements listed in Details.
#' @param ... extra arguments, currently none
#'
#' @return an instance of the `cluster_cv` class with the elements specified
#' in Details.

  cluster_cv <- function(x, ...) {

    ## input control ------

    err_txt <-
      paste("'x' has to be a named list with 'clust_analysis_object'",
            "'predictions', 'fold_stats' and 'summary' elements.")

    if(!is.list(x)) stop(err_txt, call. = FALSE)
    if(is.null(names(x))) stop(err_txt, call. = FALSE)

    if(any(!c("clust_analysis_object", "predictions", "fold_stats", "summary") %in% names(x))) {

      stop(err_txt, call. = FALSE)

    }

    ## checking the list elements -------

    if(!is_clust_analysis(x$clust_analysis_object)) {

      if(!is_combi_analysis(x$clust_analysis_object)) {

        stop(paste("'x' has to be an instance of 'clust_analysis'",
                   "or 'combi_analysis' class."),
             call. = FALSE)

      }

    }

    err_txt <-
      paste("The 'predictions' element has to be a data frame with the",
            "'observation', 'fold_clust', 'global_clust', 'correct' and",
            "'fold' variables.")

    if(!is.data.frame(x$predictions)) stop(err_txt, call. = FALSE)

    if(any(!c("observation", "fold_clust", "global_clust", "correct", "fold") %in% names(x$predictions))) {

      stop(err_txt, call. = FALSE)

    }

    err_txt <-
      paste("The 'fold_stats' element has to be a data frame with the",
            "'fold', 'accuracy', 'error', 'frac_var' and 'sil_width'",
            "variables")

    if(!is.data.frame(x$fold_stats)) stop(err_txt, call. = FALSE)

    if(any(!c('fold', 'accuracy', 'error', 'frac_var', 'sil_width') %in% names(x$fold_stats))) {

      stop(err_txt, call. = FALSE)

    }

    if(!is.data.frame(x$summary)) {

      stop("The 'summary' element has to be data frame.",
           call. = FALSE)

    }

    ## the object

    structure(x, class = 'cluster_cv')

  }

# Neighborhood -------

#' Create a `knb` class object.
#'
#' @description
#' Creates an instance of the `knb` class at the top of a data frame.
#'
#' @details
#' `knb` objects store results of testing neighborhood preservation by
#' self-organizing maps and clustering analyses and inherit most of their
#' methods from data frames.
#' The input data frame has to contain the following columns:
#'
#' * `observation` with observation identifiers
#'
#' * `clust_id` with assignment of the observations to the clusters
#'
#' * `kNN_data` with the number of nearest data point neighbors
#'
#' * `kNN_cluster` with the number of nearest cluster neighbors
#'
#' @param x a data frame as specified in Details.
#' @param ... extra arguments, currently none.
#'
#' @return an instance of the `knb` class as described in Details.

  knb <- function(x, ...) {

    err_txt <-
      paste("'x' has to be a data frame with the `observation`,",
            "'clust_id', 'kNN_data' and 'kNN_cluster'.")

    if(!is.data.frame(x)) stop(err_txt, call. = FALSE)

    if(any(!c('observation', 'clust_id', 'kNN_data', 'kNN_cluster') %in% names(x))) {

      stop(err_txt, call. = FALSE)

    }

    if(!all(is.na(x$clust_id))) {

      if(!is.factor(x$clust_id)) {

        stop("The 'clust_id' variable has to be a factor or NA.", call. = FALSE)

      }

    }

    if(!all(is.na(x$kNN_data))) {

      if(!is.numeric(x$kNN_data)) {

        stop("The 'kNN_data' variable has to be numeric.", call. = FALSE)

      }

    }

    if(!all(is.na(x$kNN_cluster))) {

      if(!is.numeric(x$kNN_cluster)) {

        stop("The 'kNN_clust' variable has to be numeric or NA.", call. = FALSE)

      }


    }

    structure(x, class = c('knb', class(x)))

  }

# Tuning results ------

#' Create an `tuner` class object.
#'
#' @description
#' Creates and object of the `tuner` class on the top of a list with tuning
#' of parameters of cluster analysis or prediction.
#'
#' @details
#' The input list has to have three elements:
#'
#' * `analysis`: a `clust_analysis` or `combi_analysis` object created with the
#' best set of the tuning parameters
#'
#' * `stats`: a data frame with values of quality stats (silhouette width,
#' fraction of potentially misclassified observations/negative silhouette width,
#' fraction of explained clustering variance, and fraction of preserved nearest
#' neighbors)
#'
#' * `fun`: name of the tuning function
#'
#' * `dataset`: a string specified which data was used during the tuning: the
#' training data set or cross-validation
#'
#' * `type`: type of analysis, development or prediction
#'
#' * `clust_vars`: a vector of names of clustering variables
#'
#' * `tune_params`: a vector of names of the tuning parameters
#'
#' * `tune_criteria`: a data frame that specifies which criteria were applied
#' to select the best combination of the tuning parameters
#'
#' * `best_tune`: a data frame storing the best values of the tuning parameters
#'
#' @param x a list with elements specified in Details.
#'
#' @return an instance of the `tuner` class as described in Details.

  tuner <- function(x) {

    ## entry control -------

    elements <- c('analysis',
                  'stats',
                  'fun',
                  'dataset',
                  'type',
                  'clust_vars',
                  'tune_params',
                  'tune_criteria',
                  'best_tune')

    err_txt <-
      paste("'x' has to be a list with the",
            paste(elements, collapse = ', '),
            "elements.")

    if(!is.list(x)) stop(err_txt, call. = FALSE)

    if(any(!elements %in% names(x))) {

      stop(err_txt, call. = FALSE)

    }

    stat_vars <- c('sil_width',
                   'frac_misclassified',
                   'frac_var',
                   'frac_np')

    err_txt <-
      paste("The `stats` element of 'x' has to be a data frame with the",
            paste(stat_vars, collapse = ', '),
            "columns")

    if(!is.data.frame(x$stats)) stop(err_txt, call. = FALSE)

    if(any(!stat_vars %in% names(x$stats))) {

      stop(err_txt, call. = FALSE)

    }

    if(!is_clust_analysis(x$analysis) & !is_combi_analysis(x$analysis)) {

      stop(paste("The `analysis` element of 'x' has to be an instance of",
                 "of the 'clust_analysis' or 'combi_analysis' class."),
           call. = FALSE)

    }

    ## class setting -------

    structure(x, class = 'tuner')

  }

# Spectralization result container -------

#' Create a `spectre` object.
#'
#' @description
#' Creates an instance of `spectre` class.
#'
#' @details
#' Technically, a `spectre` object is a list with the following components:
#'
#' * `degrees`: a numeric vector with degrees of the graph nodes stored in the
#' diagonal of the degree matrix
#'
#' * `eigen_values`: a numeric vector of eigenvalues sorted from the largest to
#' the smallest one
#'
#' * `eigen_vectors`: a numeric matrix whose rows represent the observations in
#' he initial data set and columns representing the eigenvectors.
#' The eigenvectors are sorted by their einegvalues: the largest come first
#'
#' @param x a list with elements specified in Details.
#'
#' @return an instance of the `spectre` class as described in Details.

  spectre <- function(x) {

    ## entry control -------

    error_txt <-
      paste("'x' has to be a list with the following elements: 'degree',",
            "eigen_values', and 'eigen_vectors'.")

    if(!is.list(x)) stop(error_txt, call. = FALSE)

    if(is.null(names(x))) stop(error_txt, call = FALSE)

    if(any(!c('degrees', 'eigen_values', 'eigen_vectors') %in% names(x))) {

      stop(error_txt, call. = FALSE)

    }

    if(!is.numeric(x$degrees) | !is.numeric(x$eigen_values)) {

      stop("'degree' and 'eigen_values' have to be numeric vectors.",
           call. = FALSE)

    }

    error_txt <- "'eigen_vectors' has to be a numeric matrix."

    if(!is.matrix(x$eigen_vector)) stop(error_txt, call. = FALSE)

    if(!is.numeric(x$eigen_vector)) stop(error_txt, call. = FALSE)

    ## construction --------

    structure(x, class = 'spectre')

  }

# END ------
