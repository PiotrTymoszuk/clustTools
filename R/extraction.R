# Access of elements (components) of the S3 objects

#' Extract features of a clust_analysis object.
#'
#' @description
#' A general extractor method for accessing properties and features
#' of a `clust_analysis`, `combi_analysis` and `red_analysis` object,
#' and specific methods for accessing the modeling data frame
#' and distance matrix.
#'
#' @details
#' `extract()` is a S3 generic function.
#'
#' @param x an object.
#' @param formula an object.
#' @param type the feature name:
#'
#' * `distance` extracts the matrix with distances between the observations,
#'
#' * `data` the data set used for the analysis,
#'
#' * `assignment` assignment of the observations to the clusters,
#'
#' * `clust_object` or `object` returns the wrapped clustering object.
#'
#' * `clust_object` or `scores` return the component pr score tables for the
#' observations,
#'
#' * `loadings` retrieves the table of variable loadings,
#'
#' * `sdev` returns standard deviations, associated with the
#' components.
#'
#' * `umatrix` computes the U-matrix, i.e. weighted distance between
#' the self-organizing map (SOM) nodes. Available only for clustering analyses
#' done with SOM.
#' The U matrix is computed with \code{\link[kohonen]{object.distances}}.
#'
#' @param ... extra arguments, currently none.
#'
#' @return the requested feature/property.
#' @export

  extract <- function(x, ...) UseMethod('extract')

#' @rdname extract
#' @export extract.clust_analysis
#' @export

  extract.clust_analysis <- function(x,
                                     type = c('distance',
                                              'assignment',
                                              'clust_object',
                                              'data',
                                              'object',
                                              'umatrix'), ...) {

    ## entry control

    stopifnot(is_clust_analysis(x))

    type <- match.arg(type[1],
                      c('distance',
                        'assignment',
                        'clust_object',
                        'data',
                        'object',
                        'umatrix'))

    ## output

    if(type != 'umatrix') {

      return(switch(type,
                    distance = as.dist(x$dist_mtx),
                    assignment = x$clust_assignment,
                    clust_object = x$clust_obj,
                    data = eval_tidy(x$data),
                    object = x$clust_obj))

    }

    if(!x$clust_fun %in% c('som', 'supersom')) {

      warning("The 'umatrix' option is only available for SOM analyses.",
              call. = FALSE)

      return(NULL)

    }

    kohonen::object.distances(x$clust_obj, type = 'codes')

  }

#' @rdname extract
#' @export model.frame.clust_analysis
#' @export

  model.frame.clust_analysis <- function(formula, ...) {

    stopifnot(is_clust_analysis(formula))

    eval_tidy(formula$data)

  }

#' @rdname extract
#' @export dist.clust_analysis
#' @export

  dist.clust_analysis <- function(x, type = c('distance', 'umatrix'), ...) {

    stopifnot(is_clust_analysis(x))

    type <- match.arg(type[1], c('distance', 'umatrix'))

    extract(x, type = type)

  }

#' @rdname extract
#' @export

  dist.min_analysis <- function(x, ...) {

    NextMethod()

  }

#' @rdname extract
#' @export extract.combi_analysis
#' @export

  extract.combi_analysis <- function(x,
                                     type = c('distance',
                                              'assignment',
                                              'clust_object',
                                              'data',
                                              'object',
                                              'umatrix'), ...) {

    stopifnot(is_combi_analysis(x))

    type <- match.arg(type[1],
                      c('distance',
                        'assignment',
                        'clust_object',
                        'data',
                        'object',
                        'umatrix'))

    if(type == 'assignment') {

      return(x$clust_assignment)

    } else if(type == 'umatirx') {

      return(kohonen::object.distances(x$clust_analyses$observation$clust_obj,
                                       type = 'codes'))

    } else {

      map(x$clust_analyses,
          extract,
          type = type)

    }

  }

#' @rdname extract
#' @export model.frame.combi_analysis
#' @export

  model.frame.combi_analysis <- function(formula, ...) {

    stopifnot(is_combi_analysis(formula))

    map(formula$clust_analyses, model.frame)

  }

#' @rdname extract
#' @export dist.combi_analysis
#' @export

  dist.combi_analysis <- function(x, type = c('distance', 'umatrix'), ...) {

    stopifnot(is_combi_analysis(x))

    type <- match.arg(type[1], c('distance', 'umatrix'))

    compact(extract(x, type = type))

  }

#' @rdname extract
#' @export extract.red_analysis
#' @export

  extract.red_analysis <- function(x,
                                   type = c('component_tbl',
                                            'scores',
                                            'loadings',
                                            'data',
                                            'sdev',
                                            'object'), ...) {

    stopifnot(is_red_analysis(x))

    type <- match.arg(type[1],
                      c('component_tbl',
                        'scores',
                        'loadings',
                        'data',
                        'sdev',
                        'object'))

    if(type != 'sdev') {

      return(switch(type,
                    component_tbl = x$component_tbl,
                    scores = x$component_tbl,
                    loadings = x$loadings,
                    data = eval_tidy(x$data),
                    object = x$red_obj))

    }

    if(x$red_fun == 'pca') {

      return(tibble(component = 1:length(x$red_obj$sdev),
                    sdev = x$red_obj$sdev,
                    perc_sdev = x$red_obj$sdev/sum(x$red_obj$sdev) * 100,
                    var = x$red_obj$sdev^2,
                    perc_var = x$red_obj$sdev^2/sum(x$red_obj$sdev^2) * 100))

    }

    if(x$red_fun == 'fa') {

      p <- nrow(x$red_obj$loadings)

      sum_sq <- colSums(x$red_obj$loadings^2)

      return(tibble(component = 1:ncol(x$red_obj$loadings),
                    sdev = sqrt(sum_sq/p),
                    perc_sdev = sqrt(sum_sq/p) * 100,
                    var = sum_sq/p,
                    perc_var = sum_sq/p * 100))

    }

    score_tbl <- x$component_tbl

    dims <- stringi::stri_extract(names(score_tbl), regex = 'comp_\\d+')

    dims <- dims[!is.na(dims)]

    variances <- map_dbl(score_tbl[dims], var)

    std_devs <- map_dbl(score_tbl[dims], stats::sd)

    return(tibble(component = 1:length(dims),
                  sdev = std_devs,
                  perc_sdev = std_devs/sum(std_devs) * 100,
                  var = variances,
                  perc_var = variances/sum(variances) * 100))

  }

#' @rdname extract
#' @export model.frame.red_analysis
#' @export

  model.frame.red_analysis <- function(formula, ...) {

    stopifnot(is_red_analysis(formula))

    eval_tidy(formula$data)

  }

# END --------
