# Functions used for data pre-processing

#' Set row names.
#'
#' @description Sets row names in a data frame.
#' @details a tibble is silently converted to a data frame.
#' @param data a data frame or a tibble.
#' @param row_names a character vector of the proper length.
#' @return a data frame.
#' @export

  set_rownames <- function(data, row_names = as.character(1:nrow(data))) {

    ## entry control

    stopifnot(is.data.frame(data))

    row_names <- as.character(row_names)

    stopifnot(length(row_names) == nrow(data))

    if(sum(duplicated(row_names)) > 0) {

      stop('Duplicated row names are not allowed.', call. = FALSE)

    }

    ## output

    data <- as.data.frame(data)

    rownames(data) <- row_names

    data

  }

#' Normalization of a data frame.
#'
#' @description
#' Normalization with median or mean centering of a data frame
#' or tibble (`center_data()`) or simple min/max normalization (`min_max()`).
#' Preserves the row names.
#'
#' @details
#' A wrapper around \code{\link[base]{scale}}. Mean scaling is equal
#' to canonical Z-score normalization.
#'
#' @param data a data frame or a tibble. All variables need to be numeric.
#' @param type type of the centering, mean (default) or median.
#' @param complete_cases logical, should the observations with the complete
#' variable record only be included as an output?
#'
#' @return a data frame or a tibble.
#'
#' @export

  center_data <- function(data,
                          type = c('mean', 'median'),
                          complete_cases = FALSE) {

    ## entry control

    stopifnot(is.data.frame(data))
    stopifnot(is.logical(complete_cases))

    check_numeric(data)

    type <- match.arg(type[1], c('mean', 'median'))

    ## centering

    if(complete_cases) {

      data <- filter(data, complete.cases(data))

    }

    if(!is_tibble(data)) {

      row_names <- rownames(data)

    }

    center_fun <- switch(type,
                         mean = function(x) mean(x, na.rm = TRUE),
                         median = function(x) stats::median(x, na.rm = TRUE))

    new_data <- map_dfc(data,
                        function(x) scale(x, center = center_fun(x))[, 1])

    if(!is_tibble(data)) {

      new_data <- set_rownames(new_data, row_names = row_names)

    }

    new_data

  }

#' @rdname center_data
#' @export

  min_max <- function(data,
                      complete_cases = FALSE) {

    ## entry control

    stopifnot(is.data.frame(data))
    stopifnot(is.logical(complete_cases))

    check_numeric(data)

    ## scaling

    if(complete_cases) {

      data <- filter(data, complete.cases(data))

    }

    if(!is_tibble(data)) {

      row_names <- rownames(data)

    }

    new_data <-
      map_dfc(data,
              function(x) (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))

    if(!is_tibble(data)) {

      new_data <- set_rownames(new_data, row_names = row_names)

    }

    new_data

  }

# END ------
