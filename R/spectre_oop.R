# S3 methods for the `spectre` class

# Appearance and number of observations ------

#' Print method for the `spectre` class.
#'
#' @description Prints a `spectre` object.
#'
#' @param x an object of class `spectre`.
#' @param ... extra arguments, currently none.
#'
#' @return none, called for its side effects.
#'
#' @export

  print.spectre <- function(x, ...) {

    stopifnot(is_spectre(x))

    cat(paste('SPECTRE object for n =',
              ncol(x$eigen_vectors),
              'eigenvectors'))

  }

#' Numbers of observations used for spectral decomposition.
#'
#' @description
#' Number of observations used for spectral decomposition.
#'
#' @param object an object of class `spectre`.
#' @param ... extra arguments, currently none.
#'
#' @return a numeric value.
#'
#' @export

  nobs.spectre <- function(object, ...) nrow(object$eigen_vectors)

# Type conversion -------

#' Convert an object to a `red_analysis` instance.
#'
#' @description
#' Convert an object to an instance of \code{\link{red_analysis}} class.
#'
#' @details
#' In case of \code{\link{spectre}} objects, the `kdim` trailing eigenvectors,
#' i.e. eigenvectors with the smallest eigenvalues will be used as components
#' of the `red_analysis` objects. They will be sorted in an ascending order,
#' i.e. the last eigenvector will be named 'comp_1', the one before last
#' 'comp_2' and so on.
#'
#'
#' @param x an object.
#' @param kdim number of the smallest eigenvectors included used as components
#' of the `red_analysis` object.
#' @param skip_last logical, should the last eigenvector be omitted before
#' drawing the k trailing eigenvectors?
#' @param ... extra arguments passed to the methods.
#'
#' @export

  as_red_analysis <- function(x, ...) UseMethod('as_red_analysis')

#' @rdname as_red_analysis
#' @export as_red_analysis.spectre
#' @export

  as_red_analysis.spectre <- function(x,
                                      kdim = 3,
                                      skip_last = TRUE, ...) {

    stopifnot(is_spectre(x))
    stopifnot(is.numeric(kdim))
    stopifnot(is.logical(skip_last))

    ## eigenvectors of interest ------

    eigen_vectors <- x$eigen_vectors

    observations <- rownames(eigen_vectors)

    if(skip_last) eigen_vectors <- eigen_vectors[, -ncol(eigen_vectors)]

    indexes <- seq(ncol(eigen_vectors), ncol(eigen_vectors) - kdim + 1)

    eigen_vectors <- eigen_vectors[, indexes]

    colnames(eigen_vectors) <- paste0('comp_', 1:kdim)

    observation <- NULL

    eigen_vectors <- mutate(as_tibble(eigen_vectors),
                            observation = observations)

    eigen_vectors <- relocate(eigen_vectors, observation)

    ## output --------

    model_frame <- enexpr(x)

    red_analysis(list(red_obj = NULL,
                      red_fun = 'spectralize',
                      dist_method = 'custom',
                      component_tbl = eigen_vectors,
                      loadings = NULL,
                      data = quo(!!model_frame)))

  }

# Plotting ------

#' Diagnostic plots for `spectre` objects.
#'
#' @description
#' Generates scatter plots of eigenvalues and selected eigenvectors.
#'
#' @details
#' Two types of plots are generated as specified by the `type` argument:
#'
#' * _eigenvalues_: the k smallest eigenvalues (by default the trailing 50 ones).
#' This plot type may be usefull at assessment of the number of eigenvectors
#' used in other analyses such as clustering.
#'
#' * _eigenvectors_: a scatter plot of two selected eigenvectors, be default,
#' the last and one before last eigenvector.
#'
#' @return a `ggplot` graphic object.
#'
#' @param x an object of `spectre` class.
#' @param type type of the plot as specified in the Details. `eigenvalues` by
#' default.
#' @param cust_theme custom `ggplot` `theme` object.
#' @param ... extra arguments passed to \code{\link{plot_eigenvalues}} and
#' \code{\link{plot_eigenvectors}} which control e.g. the number of the trailing
#' eigenvalues or the eigenvector pair to be plotted.
#'
#' @export plot.spectre
#' @export

  plot.spectre <- function(x,
                           type = c('eigenvalues', 'eigenvectors'),
                           cust_theme = ggplot2::theme_classic(), ...) {

    stopifnot(is_spectre(x))

    type <- match.arg(type[1], c('eigenvalues', 'eigenvectors'))

    if(!inherits(cust_theme, 'theme')) {

      stop("'cust_theme' hast to be a valid ggplot theme object.",
           call. = FALSE)

    }

    if(type == 'eigenvalues') {

      out_plot <- plot_eigenvalues(x, cust_theme = cust_theme, ...)

    } else {

      out_plot <- plot_eigenvectors(x, cust_theme = cust_theme, ...)

    }

    out_plot

  }

# END ------
