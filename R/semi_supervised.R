# Semi supervised clustering

# k-NN label propagation ------

#' Project the cluster assignment with k-NN label propagation.
#'
#' @description
#' Projects the cluster assignment with a k-nearest neighbor
#' classifier onto a new data set.
#'
#' @details
#' If a `red_analysis` object is provided as new data, the cluster
#' assignment is projected onto the component/score table. The newdata input has
#' to have the same variables as those used for development of the input
#' cluster_analysis object.
#'
#' @references
#' Leng M, Wang J, Cheng J, Zhou H, Chen X. Adaptive
#' semi-supervised clustering algorithm with label propagation.
#' J Softw Eng (2014) 8:14–22. doi:10.3923/jse.2014.14.22
#'
#' @param object a `clust_analysis` or a `combi_analysis` object
#' @param newdata a numeric data frame, matrix or a `red_analysis` object.
#' If NULL (default), the bare cluster assignment table is returned.
#' @param distance_method a distance metric, by default it is retrieved from
#' the input `clust_analysis` or `combi_analysis` object. For the later, the
#' distance used for observation clustering is used in the projection.
#' @param kNN number of the nearest neighbors.
#' @param simple_vote logical, should classical unweighted k-NN classification
#' be applied? If FALSE, distance-weighted k-NN is used with the provided kernel
#' function.
#' @param resolve_ties logical, should the ties be resolved at random? Applies
#' only to the simple unweighted voting algorithm.
#' @param kernel_fun kernel function transforming the distance into weight.
#' @param detailed logical, should a detailed output including the kNN table and
#' voting scheme be returned. If FALSE, the bare `clust_analysis` object with the
#' predictions is returned.
#'
#' @return a \code{\link{clust_analysis}} object or,
#' if `detailed` = TRUE, a list with the kNN
#' table and the voting results.
#'
#' @export

  propagate <- function(object,
                        newdata = NULL,
                        distance_method = NULL,
                        kNN = 5,
                        simple_vote = TRUE,
                        resolve_ties = FALSE,
                        kernel_fun = function(x) 1/x,
                        detailed = FALSE) {

    ## entry control -------

    if(!is_clust_analysis(object)) {

      if(!is_combi_analysis(object)) {

        stop('A valid clust_analysis or combi_analysis object required.',
             call. = FALSE)

      }

    }

    if(is.null(newdata)) {

      return(object$clust_assignment)

    }

    if(is.null(distance_method)) {

      if(is_clust_analysis(object)) {

        distance_method <- object$dist_method

      } else {

        distance_method <- object$clust_analyses$observation$dist_method

      }

    }

    if(!distance_method %in% get_kernel_info()) {

      stop('Invalid distance method.', call. = FALSE)

    }

    if(is_red_analysis(newdata)) {

      newdata <- column_to_rownames(newdata$component_tbl,
                                    'observation')

    }

    check_numeric(newdata)

    kNN <- as.integer(kNN)

    stopifnot(is.logical(simple_vote))
    stopifnot(is.logical(detailed))
    stopifnot(is.function(kernel_fun))

    clust_id <- NULL

    ## extracting the training data set -------
    ## constrained to the cases with the cluster assignment (non-NA)

    if(is_clust_analysis(object)) {

      train_set <- as.matrix(eval_tidy(object$data))

    } else {

      train_set <- as.matrix(eval_tidy(object$clust_analyses$observation$data))

    }

    if(ncol(train_set) != ncol(newdata)) {

      stop(paste('The numbers of columns in new data and the table used',
                 'for cluster development must be equal.'),
           call. = FALSE)

    }

    if(any(!names(train_set) %in% names(newdata))) {

      stop('Clustering features missing from the newdata.', call. = FALSE)

    }

    if(is.null(rownames(train_set))) {

      rownames(train_set) <- paste0('train_', 1:nrow(train_set))

    }

    train_assign <- mutate(object$clust_assignment,
                           .rowname = rownames(train_set))

    train_assign <- filter(train_assign,
                           !is.na(clust_id))

    label_vec <- set_names(train_assign$clust_id,
                           train_assign$.rowname)

    train_set <- train_set[train_assign$.rowname, ]

    ## constructing the mixed test-train matrix and calculating a distance object ---------

    if(is.null(rownames(newdata))) {

      rownames(newdata) <- paste0('test_', 1:nrow(train_set))

    }

    mix_set <- rbind(train_set,
                     as.matrix(newdata[, colnames(train_set)]))

    mix_diss <- calculate_dist(data = mix_set,
                               method = distance_method)

    ## kNN calculation and label assignment ---------

    knn_dists <- dbscan::kNN(as.dist(mix_diss),
                             k = nrow(mix_diss) - 1)

    knn_test <- list(dist = knn_dists$dist[rownames(newdata), ],
                     id = knn_dists$id[rownames(newdata), ])

    knn_test$annot_id <- matrix(rownames(mix_diss)[knn_test$id],
                                ncol = nrow(mix_diss) - 1)

    rownames(knn_test$annot_id) <- rownames(knn_test$id)

    knn_test$labels <- matrix(label_vec[knn_test$annot_id],
                              ncol = nrow(mix_diss) - 1)

    rownames(knn_test$labels) <- rownames(knn_test$id)

    ## voting, constrained to the k-nearest neighbors from the train dataset --------

    if(simple_vote) {

      clust_assignment <- map(rownames(knn_test$labels),
                              ~knn_test$labels[.x, ])

      clust_assignment <- map(clust_assignment,
                              ~vote_simple(.x[!is.na(.x)][1:kNN],
                                           resolve_ties = resolve_ties))

      clust_assignment <- unlist(clust_assignment)

    } else {

      rows <- map(rownames(knn_test$labels),
                  ~knn_test$labels[.x, ])

      dists <- map(rownames(knn_test$labels),
                   ~knn_test$dist[.x, ])

      non_na <- map(rows,
                    ~!is.na(.x))

      rows <- map2(rows, non_na, ~.x[.y][1:kNN])

      dists <- map2(dists, non_na, ~.x[.y][1:kNN])

      clust_assignment <-
        pmap(list(vector = rows,
                  dist_vec = dists),
             vote_kernel,
             kernel_fun = kernel_fun)

      clust_assignment <- unlist(clust_assignment)

    }

    model_frame <- enexpr(newdata)

    clust_out <-
      clust_analysis(list(data = quo(!!model_frame),
                          dist_mtx = as.matrix(mix_diss)[rownames(newdata), rownames(newdata)],
                          dist_method = distance_method,
                          clust_fun = 'prediction',
                          clust_obj = NULL,
                          clust_assignment = tibble(observation = rownames(knn_test$labels),
                                                    clust_id = factor(clust_assignment)),
                          dots = list2()))

    if(!detailed) {

      return(clust_out)

    } else {

      return(list(mix_data = mix_set,
                  mix_diss = mix_diss,
                  knn = knn_test[c('dist', 'id', 'labels')],
                  clust_analysis_object = clust_out))

    }

  }


# S3 methods --------

#' Semi-supervised clustering.
#'
#' @description Projects the cluster assignment onto new data using simple
#' observation matching or a k-nearest neighbor (kNN) label propagation
#' algorithm.
#'
#' @details For the implementation details, see: \code{\link{propagate}}.
#' The default distance metric is extracted from the `clust_analysis` object.
#' For `combi_analysis` objects, the default distance metric is the distance
#' between observations (not nodes!).
#' The cluster projection is done on the top level, i.e. takes into account the
#' final assignment of the observations to the clusters and ignoring
#' the SOM nodes.
#'
#' @references
#' Leng M, Wang J, Cheng J, Zhou H, Chen X. Adaptive
#' semi-supervised clustering algorithm with label propagation.
#' J Softw Eng (2014) 8:14–22. doi:10.3923/jse.2014.14.22
#'
#' @param object an object.
#' @param newdata a numeric data frame, matrix or a red_analysis object.
#' If NULL (default), the bare cluster assignment table is returned.
#'
#' @param type type of the projection: simple observation matching
#' ('class', default) or kNN label propagation ('propagation').
#'
#' @param ... extra arguments passed to \code{\link{propagate}}.
#'
#' @return a \code{\link{clust_analysis}} object.
#'
#' @export predict.clust_analysis
#' @export

  predict.clust_analysis <- function(object,
                                     newdata = NULL,
                                     type = c('class', 'propagation'), ...) {

    ## entry control ---------

    stopifnot(is_clust_analysis(object))

    type <- match.arg(type[1],
                      c('class', 'propagation'))

    if(is.null(newdata)) {

      return(object$clust_assignment)

    }

    if(is_red_analysis(newdata)) {

      newdata <- column_to_rownames(newdata$component_tbl,
                                    'observation')

    }

    check_numeric(newdata)

    ## projections -----------

    if(type == 'class') {

      train_assignment <- column_to_rownames(object$clust_assignment,
                                             'observation')

      if(nrow(newdata) != nrow(train_assignment)) {

        stop('The numbers of rows in new data and the table used for cluster development must be equal',
             call. = FALSE)

      }

      newdata <- as.data.frame(newdata)

      if(!is.null(rownames(newdata))) {

        test_assignment <-
          tibble(observation = rownames(newdata),
                 clust_id = train_assignment[rownames(newdata), 'clust_id'])

      } else {

        warning('Unnamed observations in new data', call. = FALSE)

        test_assignment <- rownames_to_column(train_assignment,
                                              'observation')

        test_assignment <- as_tibble(test_assignment)

      }

      ## output --------

      model_frame <- enexpr(newdata)

      clust_obj <-
        clust_analysis(list(data = quo(!!model_frame),
                            dist_mtx = calculate_dist(newdata,
                                                      method = object$dist_method),
                            dist_method = object$dist_method,
                            clust_fun = 'prediction',
                            clust_obj = NULL,
                            clust_assignment = test_assignment,
                            dots = list2()))

      return(clust_obj)

    }

    return(propagate(object = object,
                     newdata = newdata, ...))

  }

#' @rdname predict.clust_analysis
#' @export predict.combi_analysis
#' @export

  predict.combi_analysis <- function(object,
                                     newdata = NULL,
                                     type = c('class', 'propagation'), ...) {

    ## entry control --------

    stopifnot(is_combi_analysis(object))

    type <- match.arg(type[1],
                      c('class', 'propagation'))

    if(is.null(newdata)) {

      return(object$clust_assignment)

    }

    if(is_red_analysis(newdata)) {

      newdata <- column_to_rownames(newdata$component_tbl,
                                    'observation')

    }

    check_numeric(newdata)

    ## prediction ---------

    if(type == 'class') {

      train_assignment <- column_to_rownames(object$clust_assignment,
                                             'observation')

      if(nrow(newdata) != nrow(train_assignment)) {

        stop(paste('The numbers of rows in new data and the table used',
                   'for cluster development must be equal'),
             call. = FALSE)

      }

      newdata <- as.data.frame(newdata)

      if(!is.null(rownames(newdata))) {

        test_assignment <-
          tibble(observation = rownames(newdata),
                 clust_id = train_assignment[rownames(newdata),
                                             'clust_id'])

      } else {

        warning('Unnamed observations in new data', call. = FALSE)

        test_assignment <- rownames_to_column(train_assignment,
                                              'observation')

        train_assignment <- as_tibble(train_assignment)

      }

      ## output --------

      model_frame <- enexpr(newdata)

      clust_obj <-
        clust_analysis(list(data = quo(!!model_frame),
                            dist_mtx = calculate_dist(newdata,
                                                      method = object$clust_analyses$observation$dist_method),
                            dist_method = object$clust_analyses$observation$dist_method,
                            clust_fun = 'prediction',
                            clust_obj = NULL,
                            clust_assignment = test_assignment))

      return(clust_obj)

    }

    return(propagate(object = object,
                     newdata = newdata, ...))

  }

# END ------
