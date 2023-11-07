# SOM-specific stats

#' Quantization error and population-based convergence.
#'
#' @description
#' The functions compute two numeric statistics helpful at assessing convergence
#' of a self-organizing maps (SOM): quantization error (`qe()`) and
#' population-based convergence (`pbc()`).
#'
#' @details
#' Quantization error is computed as a sum of distances of data points to their
#' winning units (node).
#' The idea behind population-base convergence is to check whether data points
#' and SOM nodes are drawn from the same population as assessed by comparing
#' means and variances of the data points and nodes in a variable-wise manner.
#' The `qe()` and `pbc()` functions are S3 generics.
#'
#' @references
#' 1. Breard GT. Evaluating Self-Organizing Map Quality Measures as Convergence
#' Criteria Criteria. Open Access Master’s Theses. Paper 1033.
#' Available at: https://digitalcommons.uri.edu/theses/1033
#' @references
#' Kohonen T. Self-Organizing Maps. Berlin, Heidelberg: Springer Berlin
#' Heidelberg (1995). doi:10.1007/978-3-642-97610-0
#' @references
#' Hamel L, Ott BH. A Population Based Convergence Criterion for
#' Self-Organizing Maps. (2012)
#'
#' @return `qe()` returns a single numeric value. `pbc()` returns a list with
#' two components, `variable_stats` with a mean differences, variance ratios,
#' their lower and upper bounds of 95% confidence intervals and logical
#' variables indicating if the convergence was reached.
#'
#' @param x an object.
#' @param ... extra arguments passed to methods.
#'
#' @export

  qe <- function(x, ...) UseMethod('qe')

#' @rdname qe
#' @export

  qe.clust_analysis <- function(x, ...) {

    ## input control ------

    stopifnot(is_clust_analysis(x))

    if(!x$clust_fun %in% c('som', 'supersom')) {

      warning('Quantization errors are available only for SOM analyses.',
              call. = FALSE)

    }

    ## stat value ------

    sum(x$clust_obj$distances)

  }

#' @rdname qe
#' @export

  qe.combi_analysis <- function(x, ...) {

    stopifnot(is_combi_analysis(x))

    qe(x$clust_analyses$observation)

  }

#' @rdname qe
#' @export

  pbc <- function(x, ...) UseMethod('pbc')

#' @rdname qe
#' @export

  pbc.clust_analysis <- function(x, ...) {

    ## input control ------

    stopifnot(is_clust_analysis(x))

    if(!x$clust_fun %in% c('som', 'supersom')) {

      warning('Quantization errors are available only for SOM analyses.',
              call. = FALSE)

    }

    # observation and node data ------

    pbc_data <-
      list(data = model.frame(x),
           neurons = x$clust_obj$codes)

    if(x$clust_fun != 'som') {

      pbc_data <- map(pbc_data, reduce, cbind)

    }

    pbc_data <- map(pbc_data, as.data.frame)

    ## variances, variance ratios and their confidence intervals ------

    variances <- map(pbc_data, map_dbl, var)

    variance_ratios <- map2_dbl(variances[[1]], variances[[2]], `/`)

    variable <- NULL
    variance_ratio <- NULL

    f_stat <- NULL
    variance_lower <- NULL
    variance_upper <- NULL
    variance_converged <- NULL

    variance_ratios <-
      tibble(variable = names(variance_ratios),
             variance_ratio = variance_ratios)

    variance_ratios <-
      mutate(variance_ratios,
             f_stat = stats::qf(0.025,
                                nrow(pbc_data[[1]]) - 1,
                                nrow(pbc_data[[2]]) - 1,
                                lower.tail = FALSE),
             variance_lower = variance_ratio/f_stat,
             variance_upper = variance_ratio * f_stat,
             variance_converged = (1 > variance_lower) & (1 < variance_upper))

    ## means, mean difference and confidence interval of the mean delta ------

    means <- map(pbc_data, colMeans)

    mean_diff <- map2_dbl(means[[1]], means[[2]], `-`)

    mean_difference <- NULL

    mean_diff <- tibble(variable = names(mean_diff),
                        mean_difference = mean_diff)
    var_square <-
      map2_dbl(variances[[1]], variances[[2]],
               ~sqrt(.x/nrow(pbc_data[[1]]) + .y/nrow(pbc_data[[2]])))

    z_stat <- NULL
    mean_lower <- NULL
    mean_upper <- NULL
    mean_converged <- NULL

    mean_diff <-
      mutate(mean_diff,
             z_stat = stats::qnorm(0.025, lower.tail = FALSE),
             var_square = var_square,
             mean_lower = mean_difference - z_stat * var_square,
             mean_upper = mean_difference + z_stat * var_square,
             mean_converged =  (0 > mean_lower) & (0 < mean_upper))

    ## overall convergence -------

    converge <- NULL

    pbc_stats <-
      mutate(left_join(mean_diff,
                       variance_ratios,
                       by = 'variable'),
             converged = mean_converged & variance_converged)

    list(variable_stats = pbc_stats,
         convergence_stat = mean(pbc_stats$converged))

  }

#' @rdname qe
#' @export

  pbc.combi_analysis <- function(x, ...) {

    stopifnot(is_combi_analysis(x))

    pbc(x$clust_analyses$observation)

  }

# END ------
