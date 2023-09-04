# Calculation of silhouettes

# S3 methods --------

#' Silhouette statistic.
#'
#' @description
#' Computes silhouette statistics for a `clust_analysis` or `combi_analysis`
#' objects.
#'
#' @details
#' The function employs the default method of the
#' \code{\link[cluster]{silhouette}} generics and is hence agnostic to the
#' clustering method.
#' For SOM clustering, i.e. `combi_analysis` objects, the calculation is done
#' for the simple assignment of observations to the clusters
#' (SOM nodes are ignored).
#' For the extended output of the function, an object of the class
#' \code{\link{sil_extra}} is returned, which preserves cluster order
#' and names and offers tidyverse-friendly \code{\link{plot.sil_extra}}
#' and \code{\link{summary.sil_extra}} methods.
#'
#' @references
#' Rousseeuw PJ. Silhouettes: A graphical aid to the interpretation and
#' validation of cluster analysis. J Comput Appl Math (1987) 20:53–65.
#' doi:10.1016/0377-0427(87)90125-7
#' @references
#' Schubert E, Rousseeuw PJ. Faster k-Medoids Clustering: Improving the PAM,
#' CLARA, and CLARANS Algorithms. in Lecture Notes in Computer Science
#' (including subseries Lecture Notes in Artificial Intelligence and Lecture
#' Notes in Bioinformatics) (Springer), 171–187.
#' doi:10.1007/978-3-030-32047-8_16
#'
#' @param x an object of the \code{\link{clust_analysis}} or
#' \code{\link{combi_analysis}} class.
#' @param output the function output.
#' For `silhouette`, an object of the canonical cluster's class
#' \code{\link[cluster]{silhouette}} is returned.
#' For `extended`, an object of the class \code{\link{sil_extra}} is returned.
#' See the Details.
#' @param ... extra arguments passed to \code{\link[cluster]{silhouette}}.
#'
#' @return an object of the class \code{\link[cluster]{silhouette}} or
#' \code{\link{sil_extra}}.
#'
#' @export silhouette.clust_analysis
#' @export

  silhouette.clust_analysis <- function(x,
                                        output = c('extended', 'silhouette'),
                                        ...) {

    ## entry control -------

    stopifnot(is_clust_analysis(x))

    output <- match.arg(output[1], c('extended', 'silhouette'))

    ## computation --------

    sil_obj <-
      silhouette(x = as.numeric(extract(x, 'assignment')$clust_id),
                 dist = dist(x), ...)

    if(output == 'silhouette') return(sil_obj)

    sil_extra(sil_obj, extract(x, 'assignment'))

  }

#' @rdname silhouette.clust_analysis
#' @export silhouette.combi_analysis
#' @export

  silhouette.combi_analysis <- function(x,
                                        output = c('extended', 'silhouette'),
                                        ...) {

    ## entry control -------

    stopifnot(is_combi_analysis(x))

    output <- match.arg(output[1], c('extended', 'silhouette'))

    ## computation --------

    sil_obj <-
      silhouette(x = as.numeric(extract(x, 'assignment')$clust_id),
                 dist = dist(x)$observation, ...)

    if(output == 'silhouette') return(sil_obj)

    sil_extra(sil_obj, extract(x, 'assignment'))

  }

# END -----
