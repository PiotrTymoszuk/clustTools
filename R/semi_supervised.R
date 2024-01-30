# Semi supervised clustering

# k-NN label propagation ------

#' Project the cluster assignment with k-NN label propagation.
#'
#' @description
#' Projects the cluster assignment with a k-nearest neighbor
#' classifier onto a new data set.
#'
#' @details
#' If a `red_analysis` object is provided as `newdata`, the cluster
#' assignment is projected onto the component/score table. The `newdata` input
#' has to have the same variables as those used for development of the input
#' cluster_analysis object. This algorithm is not available for multi-layer
#' SOM.
#'
#' @references
#' Leng M, Wang J, Cheng J, Zhou H, Chen X. Adaptive
#' semi-supervised clustering algorithm with label propagation.
#' J Softw Eng (2014) 8:14–22. doi:10.3923/jse.2014.14.22
#'
#' @param object a `clust_analysis` or a `combi_analysis` object
#' @param newdata a numeric data frame, matrix or a `red_analysis` object.
#' If NULL (default), the bare cluster assignment table is returned.
#' @param variables an optional vector with names of variables to be used for
#' the cluster assignment prediction. If `NULL` (default), all variables will
#' be used.
#' @param active_variables logical, should the prediction be done with active
#' variables only? refers only to objects created with hard threshold
#' regularized algorithms.
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

  propagate <- function(object,
                        newdata = NULL,
                        variables = NULL,
                        active_variables = FALSE,
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

    } else {

      if(object$clust_fun == 'supersom') {

        stop('kNN propagarion not implemented for multi-layer SOM.',
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
    stopifnot(is.logical(active_variables))

    clust_id <- NULL

    ## extracting the training data set -------
    ## constrained to the cases with the cluster assignment (non-NA)

    if(is_clust_analysis(object)) {

      train_set <- as.matrix(model.frame(object))

    } else {

      train_set <- as.matrix(model.frame(object)$observation)

    }

    if(active_variables) {

      if(is_clust_analysis(object)) {

        variables <- object$active_variables

      } else {

        variables <- object$clust_analyses$node$active_variables

      }

    }

    if(!is.null(variables)) {

      stopifnot(is.character(variables))

      if(any(!variables %in% colnames(train_set)) |
         any(!variables %in% colnames(newdata))) {

        stop(paste("Some features specified by the 'variables' vector",
                   "are missing from the data."),
             call. = FALSE)

      }

      train_set <- train_set[, variables, drop = FALSE]

      test_set <- newdata[, variables, drop = FALSE]

    } else {

      test_set <- newdata

    }

    if(ncol(train_set) != ncol(test_set)) {

      stop(paste('The numbers of columns in newdata and the data used',
                 'for cluster development must be equal.'),
           call. = FALSE)

    }

    if(any(!colnames(train_set) %in% colnames(test_set))) {

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

    if(is.null(rownames(test_set))) {

      rownames(test_set) <- paste0('test_', 1:nrow(train_set))

    }

    mix_set <- rbind(train_set,
                     as.matrix(test_set[, colnames(train_set)]))

    mix_diss <- calculate_dist(data = mix_set,
                               method = distance_method)

    ## kNN calculation and label assignment ---------

    knn_dists <- dbscan::kNN(as.dist(mix_diss),
                             k = nrow(mix_diss) - 1)

    knn_test <- list(dist = knn_dists$dist[rownames(test_set), ],
                     id = knn_dists$id[rownames(test_set), ])

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

    clust_ass <-
      tibble(observation = rownames(knn_test$labels),
             clust_id = factor(clust_assignment,
                               levels(object$clust_assignment$clust_id)))

    clust_out <-
      clust_analysis(list(data = quo(!!model_frame),
                          dist_mtx = calculate_dist(newdata, distance_method),
                          dist_method = distance_method,
                          clust_fun = 'prediction',
                          clust_obj = NULL,
                          clust_assignment = clust_ass,
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


# SOM prediction ---------

#' Predict node assignment for self-organizing map.
#'
#' @description
#' Predictions of assignment of the observations to the self-organizing map
#' SOM nodes are made with the trained SOM neuronal network. The distances,
#' weights and SOM architecture are extracted from the provided `clust_analysis`
#' or `combi_analysis` object.
#'
#' @details
#' Uses \code{\link[kohonen]{map.kohonen}} for mapping of the observations onto
#' the SOM nodes.
#' If a `red_analysis` object is provided as `newdata`, the cluster
#' assignment is predicted for the component/score table. The `newdata` input
#' has to have the same variables as those used for development of the input
#' cluster_analysis object.
#'
#' @references
#' Kohonen T. Self-Organizing Maps. Berlin, Heidelberg: Springer Berlin
#' Heidelberg (1995). doi:10.1007/978-3-642-97610-0
#'
#' @references
#' Wehrens R, Kruisselbrink J. Flexible self-organizing maps in kohonen 3.0.
#' J Stat Softw (2018) 87:1–18. doi:10.18637/jss.v087.i07
#'
#' @return a \code{\link{clust_analysis}} object.
#'
#' @param object a `clust_analysis` or a `combi_analysis` object.
#' @param newdata a numeric data frame, matrix or a `red_analysis` object or a
#' list containing such objects.
#' If NULL (default), the bare cluster assignment table is returned.
#' @param ... extra arguments passed to \code{\link[kohonen]{map.kohonen}}.

  map_som <- function(object, newdata = NULL, ...) {

    ## entry control: the object -------

    if(!is_clust_analysis(object)) {

      if(!is_combi_analysis(object)) {

        stop('A valid clust_analysis or combi_analysis object required.',
             call. = FALSE)

      }

    }

    if(is_clust_analysis(object)) {

      if(object$clust_fun != 'som') {

        warning(paste("SOM predictions are only possible for clustering",
                      "objects generated with the SOM algorithm."),
                call. = FALSE)

        return(NULL)

      }

    }

    if(is.null(newdata)) {

      return(object$clust_assignment)

    }

    if(is_clust_analysis(object)) {

      distance_method <- object$dist_method

    } else {

      distance_method <- object$clust_analyses$observation$dist_method

    }

    ## entry control: newdata -------

    if(is_red_analysis(newdata)) {

      newdata <- column_to_rownames(newdata$component_tbl,
                                    'observation')

    }

    check_numeric(newdata)

    ## checking for compatibility of the training and test data -----

    if(is_clust_analysis(object)) {

      train_set <- model.frame(object)

    } else {

      train_set <- model.frame(object)$observation

    }

    train_set <- as.matrix(train_set)

    newdata <- as.matrix(newdata)

    if(ncol(train_set) != ncol(newdata)) {

      stop(paste('The numbers of columns in newdata and the data used',
                 'for cluster development must be equal.'),
           call. = FALSE)

    }

    if(any(!colnames(train_set) %in% colnames(newdata))) {

      stop('Clustering features missing from the newdata.', call. = FALSE)

    }

    if(is.null(rownames(train_set))) {

      rownames(train_set) <- paste0('train_', 1:nrow(train_set))

    }

    if(is.null(rownames(newdata))) {

      rownames(newdata) <- paste0('train_', 1:nrow(newdata))

    }

    ## node mapping, clust_analysis ------

    if(is_clust_analysis(object)) {

      preds <- kohonen::map(object$clust_obj, newdata = newdata, ...)

      observation <- NULL
      clust_id <- NULL
      node <- NULL
      neuro_dist <- NULL

      clust_ass <-
        tibble(observation = rownames(newdata),
               clust_id = factor(preds$unit.classif),
               node = preds$unit.classif,
               neuro_dist = preds$distances)

    }

    ## node and cluster mapping, combi_analysis -------

    if(is_combi_analysis(object)) {

      preds <- kohonen::map(object$clust_analyses$observation$clust_obj,
                            newdata = newdata, ...)

      observation <- NULL
      node <- NULL

      clust_ass <-
        tibble(observation = rownames(newdata),
               node = as.character(preds$unit.classif))

      node_ass <- filter(object$clust_assignment,
                         !duplicated(node))

      clust_ass <- left_join(clust_ass,
                             node_ass[c('node', 'clust_id')],
                             by = 'node')

    }

    ## cluster analysis object -----

    model_frame <- enexpr(newdata)

    dist_matrix <- calculate_dist(newdata, distance_method)

    clust_analysis(list(data = quo(!!model_frame),
                        dist_mtx = dist_matrix,
                        dist_method = distance_method,
                        clust_fun = 'prediction',
                        clust_obj = NULL,
                        clust_assignment = clust_ass,
                        dots = list2()))

  }

#' @rdname map_som

  map_supersom <- function(object, newdata = NULL, ...) {

    ## entry control: the object -------

    if(!is_clust_analysis(object) & !is_umatrix_analysis(object)) {

      stop('A valid clust_analysis or combi_analysis object required.',
           call. = FALSE)

    }

    if(is_clust_analysis(object)) {

      if(object$clust_fun != 'supersom') {

        warning(paste("SOM predictions are only possible for clustering",
                      "objects generated with the SOM algorithm."),
                call. = FALSE)

        return(NULL)

      }

    }

    if(is.null(newdata)) {

      return(object$clust_assignment)

    }

    # entry control: newdata --------

    if(is.matrix(newdata) | is.data.frame(newdata) | is_red_analysis(newdata)) {

      stop(paste("'newdata' has to be a list of numeric data frames, matrices",
                 "or 'red_analysis' objects."),
           call. = FALSE)

    }

    newdata <-
      map(newdata,
          function(x) if(is_red_analysis(x)) column_to_rownames(newdata$component_tbl,
                                                                'observation') else x)

    purrr::walk(newdata, check_numeric)

    newdata <- map(newdata, as.matrix)

    ## the compatibility check is done by map.kohonen()

    ## predictions ------

    preds <- kohonen::map(object$clust_obj, newdata = newdata, ...)

    observation <- NULL
    clust_id <- NULL
    node <- NULL
    neuro_dist <- NULL

    clust_ass <-
      tibble(observation = rownames(newdata[[1]]),
             clust_id = factor(preds$unit.classif),
             node = preds$unit.classif,
             neuro_dist = preds$distances)

    ## weighted distance matrix ------

    distance_names <- object$clust_obj$dist.fcts

    distance_weights <-
      object$clust_obj$user.weights * object$clust_obj$distance.weights

    dist_matrix <-
      calculate_weighted_dist(newdata,
                              method = distance_names,
                              weights = distance_weights,
                              FUN = `+`)

    ## the output clust_analysis object -------

    model_frame <- enexpr(newdata)

    clust_analysis(list(data = quo(!!model_frame),
                        dist_mtx = dist_matrix,
                        dist_method = 'weighted_som',
                        clust_fun = 'supersom_prediction',
                        clust_obj = NULL,
                        clust_assignment = clust_ass,
                        dots = list2(dist.fcts = distance_names,
                                     user.weights = object$clust_obj$user.weights,
                                     distance.weights = object$clust_obj$distance.weights)))

  }

# S3 methods --------

#' Semi-supervised clustering.
#'
#' @description Projects the cluster assignment onto new data using simple
#' observation matching, a k-nearest neighbor (kNN) label propagation
#' algorithm or, specifically for self-organizing map (SOM), predicts the node
#' assignment based on the trained SOM neuronal network.
#'
#' @details For the implementation details of the kNN label propagation
#' algorithm, see: \code{\link{propagate}}.
#'
#' The default distance metric is extracted from the `clust_analysis` object.
#' For `combi_analysis` objects, the default distance metric is the distance
#' between observations (not nodes!).
#' In case of clustering analyses performed with hard threshold regularization
#' algorithms (currently only \code{\link{htk_cluster}}), the prediction by
#' the kNN classifier is done by default by using all available variables.
#' However, by setting `active_variables = TRUE`, the user may switch to
#' prediction of the cluster assignment only with variables contributing to
#' development of the clustering structure. See the paper by Raymaekers and
#' Zamar for rationale of the hard thresholding regularization and active
#' variable selection.
#'
#' Currently, it is not possible to perform semi-supervised clustering for
#' clustering analysis objects generated with user-provided dissimilarity
#' objects (subclass `min_analysis` of `clust_analysis`). In such cases, `NULL`
#' is returned with a warning.
#'
#' For the kNN propagation, the cluster projection is done on the top level,
#' i.e. takes into account the final assignment of the observations to the
#' clusters and ignoring the SOM nodes.
#' Predictions via the trained SOM neuronal network are accomplished with
#' \code{\link{map_som}} (single-layer SOM) or \code{\link{map_supersom}}
#' (multi-layer SOM), which internally use
#' the \code{\link[kohonen]{map.kohonen}} function. In this case, the distances,
#' weights and SOM architecture are extracted from the `clust_analysis` or
#' `combi_analysis` object.
#' If the SOM prediction method is applied to a `combi_analysis` object, the
#' cluster assignment is done in a bottom - top direction: the observations
#' are mapped onto the SOM nodes and the nodes assigned to the clusters as
#' specified by the assignment data frame (component `clust_assignment` of
#' the `combi_analysis` object).
#' If the user tries to apply the SOM method with a `clust_analysis` method
#' that was not generated with a non-SOM algorithm, `NULL` is returned
#' with a warning.
#' The SOM method is also the only method applicable to analyses employing
#' multi-layer SOM.
#'
#'
#' @references
#' Leng M, Wang J, Cheng J, Zhou H, Chen X. Adaptive
#' semi-supervised clustering algorithm with label propagation.
#' J Softw Eng (2014) 8:14–22. doi:10.3923/jse.2014.14.22
#' @references
#' Kohonen T. Self-Organizing Maps. Berlin, Heidelberg: Springer Berlin
#' Heidelberg (1995). doi:10.1007/978-3-642-97610-0
#' @references
#' Wehrens R, Kruisselbrink J. Flexible self-organizing maps in kohonen 3.0.
#' J Stat Softw (2018) 87:1–18. doi:10.18637/jss.v087.i07
#' @references
#' Raymaekers J, Zamar RH. Regularized K-means Through Hard-Thresholding.
#' J Mach Learn Res (2022) 23:1–48. Available at:
#' http://jmlr.org/papers/v23/21-0052.html
#'
#' @param object an object.
#' @param newdata a numeric data frame, matrix or a red_analysis object.
#' If NULL (default), the bare cluster assignment table is returned.
#' @param type type of the projection: simple observation matching
#' ('class', default), kNN label propagation ('propagation') or prediction
#' via SOM neuronal network ('som'). The SOM prediction method is the sole
#' prediction algorithm implemented for multi-layer SOM.
#' @param active_variables logical, should only active variables be used for the
#' cluster assignment prediction? Relevant only for objects created with hard threshold regularization algorithms and ignored otherwise  See Details.
#'
#' @param ... extra arguments passed to \code{\link{propagate}}.
#'
#' @return a \code{\link{clust_analysis}} object.
#'
#' @export predict.clust_analysis
#' @export

  predict.clust_analysis <- function(object,
                                     newdata = NULL,
                                     type = c('class', 'propagation', 'som'),
                                     active_variables = FALSE,
                                     ...) {

    ## entry control ---------

    stopifnot(is_clust_analysis(object))

    type <- match.arg(type[1],
                      c('class', 'propagation', 'som'))

    if(is.null(newdata)) {

      return(object$clust_assignment)

    }

    if(object$clust_fun != 'supersom') {

      if(is_red_analysis(newdata)) {

        newdata <- column_to_rownames(newdata$component_tbl,
                                      'observation')

      }

      check_numeric(newdata)

    } else {

      if(is.data.frame(newdata) | is.matrix(newdata) | is_red_analysis(newdata)) {

        stop(paste("'newdata' has to be a list of numeric data frames,",
                   "matrices or 'red_analysis' objects."),
             call. = FALSE)

      }

      newdata <-
        map(newdata,
            function(x) if(is_red_analysis(x)) column_to_rownames(newdata$component_tbl,
                                                                  'observation') else x)


      purrr::walk(newdata, check_numeric)

    }

    stopifnot(is.logical(active_variables))

    ## projections -----------

    if(type %in% c('class', 'propagation') & object$clust_fun == 'supersom') {

      warning(paste("Simple class assignment and kNN propagation are",
                    "not implemented for multi-layer SOM."),
              call. = FALSE)

      return(NULL)

    }

    if(type == 'class') {

      train_assignment <- column_to_rownames(object$clust_assignment,
                                             'observation')

      if(nrow(newdata) != nrow(train_assignment)) {

        stop('The numbers of rows in newdata and the table used for cluster development must be equal',
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

    ## cluster assignment via propagation and SOM ------

    if(type == 'propagation') {

      return(propagate(object = object,
                       newdata = newdata,
                       active_variables = active_variables, ...))

    } else if(type == 'som') {

      if(object$clust_fun == 'som') {

        return(map_som(object = object,
                       newdata = newdata, ...))

      } else {

        return(map_supersom(object = object,
                            newdata = newdata, ...))

      }

    }

  }

#' @rdname predict.clust_analysis
#' @export

  predict.min_analysis <- function(object, ...) {

    warning(paste('Currently, semi-supervised clustering for objects',
                  'generated with user-provided dissimilarity matrices',
                  'is not supported.'),
            call. = FALSE)

    return(NULL)

  }

#' @rdname predict.clust_analysis
#' @export predict.combi_analysis
#' @export

  predict.combi_analysis <- function(object,
                                     newdata = NULL,
                                     type = c('class', 'propagation', 'som'),
                                     active_variables = FALSE, ...) {

    ## entry control --------

    stopifnot(is_combi_analysis(object))

    type <- match.arg(type[1],
                      c('class', 'propagation', 'som'))

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

    switch(type,
           propagation = propagate(object = object,
                                   newdata = newdata,
                                   active_variables = active_variables, ...),
           som = map_som(object = object,
                         newdata = newdata, ...))

  }

#' @rdname predict.clust_analysis
#' @export

  predict.umatrix_analysis <- function(object,
                                       newdata = NULL, ...) {

    ## entry control --------

    stopifnot(is_combi_analysis(object))
    stopifnot(is_umatrix_analysis(object))

    ## input check for 'newdata' is done by the downstream function
    ## SOM node predictions --------

    som_preds <- map_supersom(object$clust_analyses$observation,
                              newdata, ...)

    ## node - cluster assignment ------

    node <- NULL

    assignments <- map(list(newdata = som_preds,
                            object = object),
                       extract, 'assignment')

    assignments <- map(assignments,
                       mutate,
                       node = as.character(node))

    assignments <-
      map2(assignments,
           list(c('observation', 'node'),
                c('node', 'clust_id')),
           ~.x[.y])

    assignments$object <- filter(assignments$object,
                                 !duplicated(node))

    clust_ass <-
      reduce(assignments,
             left_join,
             by = 'node')

    ## the output cluster analysis object --------

    som_preds$clust_assignment <- clust_ass

    som_preds

  }

# END ------
