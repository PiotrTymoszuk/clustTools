# Appearance for the S3 class objects

#' Printing of objects.
#'
#' @description
#' Prints a `clust_analysis`, `combi_analysis`, `red_analysis`
#' or `cross_dist` object.
#'
#' @param x an object.
#' @param ... extra arguments, currently none.
#'
#' @return nothing, called for side effects.
#' @export

  print.clust_analysis <- function(x, ...) {

    stopifnot(is_clust_analysis(x))

    cat(paste0('Clustering analysis with ',
               toupper(x$clust_fun), ' and ',
               x$dist_method, ' distance method.'))

    cat('\nCluster assignment:\n')

    print(as_tibble(x$clust_assignment))

  }

#' @rdname print.clust_analysis
#' @export

  print.combi_analysis <- function(x, ...) {

    stopifnot(is_combi_analysis(x))

    cat(paste('Combined SOM/clusetring analysis with',
               x$dist_method, 'distance method.'))

    cat('\nCluster assignment:\n')

    print(as_tibble(x$clust_assignment))

  }

#' @rdname print.clust_analysis
#' @export

  print.cross_dist <- function(x, ...) {

    stopifnot(is_cross_dist(x))

    cat(paste0(attr(x, 'type'), ' ',
               attr(x, 'dist_method'),
               ' cross-distances between clusters.'))

    cat('\nClusters:\n')
    cat(paste(names(x), collapse = ', '))

  }

#' @rdname print.clust_analysis
#' @export

  print.red_analysis <- function(x, ...) {

    stopifnot(is_red_analysis(x))

    cat(paste0(toupper(x$red_fun), ' reduction analysis object.'))

    cat('\nComponents:\n')

    print(as_tibble(x$component_tbl))

  }

#' @rdname print.clust_analysis
#' @export

  print.tuner <- function(x, ...) {

    stopifnot(is_tuner(x))

    cat(paste('Tuning results for',
              paste(x$tune_params, collapse = ', ')))

    cat('\nBest parameter combination:\n')

    print(as_tibble(x$best_tune))

  }

# END ------
