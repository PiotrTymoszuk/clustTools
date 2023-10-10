# Set custom cluster names

#' Set cluster names.
#'
#' @description
#' Sets custom cluster names.
#'
#' @details
#' The package's clustering functions of the clustTools package
#' name clusters with integer numbers by default. This method poses a handy
#' tool to set custom cluster names with a named character vector.
#' The cluster order (i.e. vector levels) is defined by the order of
#' the naming vector's elements.
#'
#' @param .data a `clust_analysis` or `combi_analysis` object.
#' @param nm a named character vector with the new names as elements and old
#' cluster names as names.
#' @param ... extra arguments, currently none.
#'
#' @return a \code{\link{clust_analysis}} or \code{\link{combi_analysis}} object.
#'
#' @export rename.clust_analysis
#' @export

  rename.clust_analysis <- function(.data, nm, ...) {

    ## entry control -------

    stopifnot(is_clust_analysis(.data) | is_combi_analysis(.data))

    clust_ass <- extract(.data, 'assignment')

    err_txt <- paste("Names of the 'nm' vector must match old",
                     "names of the clustering object.")

    if(!all(names(nm) %in% levels(clust_ass$clust_id))) {

      stop(err_txt, call. = FALSE)

    }

    if(length(nm) != length(levels(clust_ass$clust_id))) {

      stop(err_txt, call. = FALSE)

    }

    ## renaming --------

    clust_id <- NULL

    .data$clust_assignment <-
      mutate(.data$clust_assignment,
             clust_id = nm[as.character(clust_id)],
             clust_id = unname(clust_id),
             clust_id = factor(clust_id, unname(nm)))

    .data

  }

#' @rdname rename.clust_analysis
#' @export rename.combi_analysis
#' @export

  rename.combi_analysis <- function(.data, nm, ...) {

    rename.clust_analysis(.data, nm)

  }

# END ------
