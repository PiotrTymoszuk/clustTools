# Functions for semi supervised clustering.

# k-NN label propagation ------

#' Project the cluster assignment with k-NN label propagation.
#'
#' @description Projects the cluster assignment with a k-nearest neighbor
#' classifier onto a new data set.
#' @details If a red_analysis object is provided as new data, the cluster
#' assignment is projected onto the component/score table. The newdata input has
#' to have the same variables as those used for development of the input
#' cluster_analysis object.
#' @param object a clust_analysis or a combi_analysis object
#' @param newdata a numeric data frame, matrix or a red_analysis object. If NULL
#' (default), the bare cluster assignment table is returned.
#' @param distance_method a distance metric, by default it is retrieved from
#' the input clust_analysis or combi_analysis object. For the later, the
#' distance used for observation clustering is used in the projection.
#' @param kNN number of the nearest neighbors.
#' @param simple_vote logical, should classical unweighted k-NN classification
#' be applied? If FALSE, distance-weighted k-NN is used with the provided kernel
#' function.
#' @param resolve_ties logical, should the ties be resolved at random? Applies
#' only to the simple unweighted voting algorithm.
#' @param kernel_fun kernel function transforming the distance into weight.
#' @param detailed logical, should a detailed output including the kNN table and
#' voting scheme be returned. If FALSE, the bare clust_analysis object with the
#' predictions is returned.
#' @return a clust_analysis object or, if detailed = TRUE, a list with the kNN
#' table and the voting results.
#' @export

  propagate <- function(object,
                        newdata = NULL,
                        distance_method = NULL,
                        kNN = 5,
                        simple_vote = TRUE,
                        resolve_ties = FALSE,
                        kernel_fun = function(x) 1/x,
                        detailed = FALSE) {

    ## entry control

    if(!all(class(object) == 'clust_analysis')) {

      if(!all(class(object) == 'combi_analysis')) {

        stop('A valid clust_analysis or combi_analysis object required')

      }

    }

    if(is.null(newdata)) {

      return(object$clust_assignment)

    }

    if(is.null(distance_method)) {

      if(all(class(object) == 'clust_analysis')) {

        distance_method <- object$dist_method

      } else {

        distance_method <- object$clust_analyses$observation$dist_method

      }

    }

    if(!distance_method %in% clustTools::get_kernel_info()) {

      stop('Invalid distance method.', call. = FALSE)

    }

    if(all(class(newdata) == 'red_analysis')) {

      newdata <- tibble::column_to_rownames(newdata$component_tbl,
                                            'observation')

    }

    clustTools:::check_numeric(newdata)

    kNN <- as.integer(kNN)

    stopifnot(is.logical(simple_vote))
    stopifnot(is.logical(detailed))
    stopifnot(is.function(kernel_fun))

    ## extracting the training data set
    ## constrained to the cases with the cluster assignment (non-NA)

    if(all(class(object) == 'clust_analysis')) {

      train_set <- as.matrix(rlang::eval_tidy(object$data))

    } else {

      train_set <- as.matrix(rlang::eval_tidy(object$clust_analyses$observation$data))

    }

    if(ncol(train_set) != ncol(newdata)) {

      stop('The numbers of columns in new data and the table used for cluster development must be equal.',
           call. = FALSE)

    }

    if(any(!names(train_set) %in% names(newdata))) {

      stop('Clustering features missing from the newdata.', call. = FALSE)

    }

    if(is.null(rownames(train_set))) {

      rownames(train_set) <- paste0('train_', 1:nrow(train_set))

    }

    train_assign <- dplyr::mutate(object$clust_assignment,
                                  .rowname = rownames(train_set))

    train_assign <- dplyr::filter(train_assign,
                                  !is.na(clust_id))

    label_vec <- rlang::set_names(train_assign$clust_id,
                                  train_assign$.rowname)

    train_set <- train_set[train_assign$.rowname, ]

    ## constructing the mixed test-train matrix and calculating a distance object

    if(is.null(rownames(newdata))) {

      rownames(newdata) <- paste0('test_', 1:nrow(train_set))

    }

    mix_set <- rbind(train_set,
                     as.matrix(newdata[, colnames(train_set)]))

    mix_diss <- clustTools::calculate_dist(data = mix_set,
                                           method = distance_method)

    ## kNN calculation and label assignment

    knn_dists <- dbscan::kNN(as.dist(mix_diss), k = nrow(mix_diss) - 1)

    knn_test <- list(dist = knn_dists$dist[rownames(newdata), ],
                     id = knn_dists$id[rownames(newdata), ])

    knn_test$annot_id <- matrix(rownames(mix_diss)[knn_test$id],
                                ncol = nrow(mix_diss) - 1)

    rownames(knn_test$annot_id) <- rownames(knn_test$id)

    knn_test$labels <- matrix(label_vec[knn_test$annot_id],
                              ncol = nrow(mix_diss) - 1)

    rownames(knn_test$labels) <- rownames(knn_test$id)

    ## voting, constrained to the k-nearest neighbors from the train dataset

    if(simple_vote) {

      clust_assignment <- purrr::map(rownames(knn_test$labels),
                                     ~knn_test$labels[.x, ])

      clust_assignment <- purrr::map(clust_assignment,
                                     ~clustTools:::vote_simple(.x[!is.na(.x)][1:kNN],
                                                               resolve_ties = resolve_ties))

      clust_assignment <- unlist(clust_assignment)

    } else {

      rows <- purrr::map(rownames(knn_test$labels),
                         ~knn_test$labels[.x, ])

      dists <- purrr::map(rownames(knn_test$labels),
                          ~knn_test$dist[.x, ])

      non_na <- purrr::map(rows,
                           ~!is.na(.x))

      rows <- purrr::map2(rows, non_na, ~.x[.y][1:kNN])

      dists <- purrr::map2(dists, non_na, ~.x[.y][1:kNN])

      clust_assignment <- purrr::pmap(list(vector = rows,
                                           dist_vec = dists),
                                      clustTools:::vote_kernel,
                                      kernel_fun = kernel_fun)

      clust_assignment <- unlist(clust_assignment)

    }

    model_frame <- rlang::enexpr(newdata)

    clust_out <- clustTools::clust_analysis(list(data = rlang::quo(!!model_frame),
                                                 dist_mtx = as.matrix(mix_diss)[rownames(newdata), rownames(newdata)],
                                                 dist_method = distance_method,
                                                 clust_fun = 'prediction',
                                                 clust_obj = NULL,
                                                 clust_assignment = tibble::tibble(observation = rownames(knn_test$labels),
                                                                                   clust_id = factor(clust_assignment)),
                                                 dots = rlang::list2()))

    if(!detailed) {

      return(clust_out)

    } else {

      return(list(mix_data = mix_set,
                  mix_diss = mix_diss,
                  knn = knn_test[c('dist', 'id', 'labels')],
                  clust_analysis_object = clust_out))

    }

  }

# END -----
