# Utilities and front-end functions for matrix decomposition

# Computation utilities --------

#' Compute an affinity matrix.
#'
#' @description
#' The function computes an affinity matrix given a dissimilarity/distance
#' matrix.
#'
#' @details
#' The code is inspired by two sources:
#' https://rpubs.com/gargeejagtap/SpectralClustering
#' http://www.di.fc.ul.pt/~jpn/r/spectralclustering/spectralclustering.html.
#'
#' In brief, affinity matrices A in two flavors are computed:
#'
#' * _weighted_: if observations `i` and `j` are nearest neighbors as defined by
#' `kNN`, `A[i,j] = similarity(i,j)`. Otherwise, `A[i,j] = 0`
#'
#' * _unweighted_: if observations `i` and `j` are nearest neighbors as
#' defined by `kNN`, `A[i,j] = 1`. Otherwise, `A[i,j] = 0`
#'
#' Note, that the conversion of of pairwise distances to pairwise similarities
#' is done with the user-provided function `simil_fun`. The default `siml_fun`
#' is a popular transformation of Euclidean distances to Euclidean similarities.
#' It is, however, recommended to experiment with other functions such as
#' `1 - x` for binary distances, `2 - x` for cosine distances, or even with
#' some common kernel functions (e.g. Gaussian kernel).
#'
#' @param dist_mtx a square distance matrix.
#' @param kNN numeric, number of the nearest neighbors., has to be lower than
#' the dimension of `dist_mtx`.
#' @param weighted logical, should the affinity matrix be weighted by
#' similarity? Defaults to `TRUE`.
#' @param simil_fun a function used to convert pairwise distances
#' to pairwise similarities.
#'
#' @return a numeric matrix.
#'
#' @export

  dist2affi <- function(dist_mtx,
                        kNN = 5,
                        weighted = TRUE,
                        simil_fun = function(x) 1/(1 + x)) {

    ## converts a distance matrix to an affinity matrix:
    ## unweighted (neighborhood 0/1 coded)
    ## or similarity-weighted (neighbors get the similarity, non-neighbors 0)

    ## entry control -------

    error_txt <- "'dist_mtx' has to be a numeric matrix."

    if(!is.matrix(dist_mtx)) stop(error_txt, call. = FALSE)

    if(!is.numeric(dist_mtx)) stop(error_txt, call. = FALSE)

    stopifnot(is.numeric(kNN))

    kNN <- as.integer(kNN)

    stopifnot(is.logical(weighted))

    if(!is.function(simil_fun)) {

      stop("'simil_fun' has to be a function.", call. = FALSE)

    }

    ## affinity matrix calculation ------

    n <- nrow(dist_mtx)

    A <- matrix(0, nrow = n, ncol = n)

    stopifnot(kNN < n)

    if(!weighted) {

      ## A: stores 0/1 coded information if two nodes are
      ## in the neighborhood

      for (i in 1:n) {

        # the nearest neighbors get 1

        nn_index <- order(dist_mtx[i,])[2:kNN]

        A[i,][nn_index] <- 1

      }

      # find mutual neighbors

      A <- A + t(A)
      A[A == 2] <- 1

    } else {

      S <- simil_fun(dist_mtx)

      for(i in 1:n) {

        best.similarities <- sort(S[i,], decreasing = TRUE)[2:kNN]

        for (s in best.similarities) {

          # undirected graph -> symmetric matrix

          j <- which(S[i,] == s)

          A[i,j] <- S[i,j]
          A[j,i] <- S[i,j]

        }

      }

    }

    dimnames(A) <- dimnames(dist_mtx)

    return(A)

  }

# Plotting utils ----------

#' Plot trailing eigenvalues.
#'
#' @description
#' Plots k least eigenvalues.
#'
#' @return a `ggplot` graphic.
#'
#' @param spectre_object an object of `spectre` class
#' @param k numeric that specifies the number of the least eigenvalues
#' to be plotted.
#' @param plot_title plot title.
#' @param point_size size of the data points.
#' @param point_color point color.
#' @param point_alpha point alpha.
#' @param cust_theme custom `ggplot` `theme` object.

  plot_eigenvalues <- function(spectre_object,
                               k = 50,
                               plot_title = NULL,
                               point_size = 2,
                               point_color = 'steelblue',
                               point_alpha = 0.75,
                               cust_theme = ggplot2::theme_classic()) {

    ## entry control --------

    stopifnot(is_spectre(spectre_object))

    stopifnot(is.numeric(k))

    if(!inherits(cust_theme, 'theme')) {

      stop("'cust_theme' has to be a valid ggplot theme object.",
           call. = FALSE)

    }

    ## plot data and labels -------

    e_values <- spectre_object$eigen_values

    n <- length(e_values)

    indexes <- seq(n - k, n)

    if(is.null(plot_title)) plot_title <- 'Smallest eigenvalues'

    plot_subtitle <-
      paste0('eigenvalues, total: n = ', n,
             ', displayed: n = ', length(indexes))

    value_number <- NULL
    value <- NULL

    plot_data <-
      tibble(value_number = indexes,
             value = e_values[indexes])

    ## plotting -------

    ggplot(plot_data,
           aes(x = value_number,
               y = value)) +
      geom_point(size = point_size,
                 shape = 21,
                 alpha = point_alpha,
                 fill = point_color) +
      cust_theme +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = 'Eigenvalue number',
           y = 'Eigenvalue')

  }

#' Scatter plot of selected eigenvectors.
#'
#' @description
#' Generates a two-dimensional scatter plot for selected eigenvectors.
#'
#' @return a `ggplot` graphic object.
#'
#' @param spectre_object an object of `spectre` class
#' @param x_eigen a numeric. The number of eigenvector to be presented in
#' the x axis. If `NULL`, this is the last eigenvector available in the object.
#' @param y_eigen a numeric. The number of eigenvector to be presented in the
#' y axis. If `NULL`, this is the one before last eigenvectro available in
#' the object.
#' @param plot_title plot title.
#' @param point_size size of the data points.
#' @param point_color point color.
#' @param point_alpha point alpha.
#' @param point_wjitter jittering width for the data points.
#' @param point_hjitter jittering height for the data points.
#' @param cust_theme custom `ggplot` `theme` object.

  plot_eigenvectors <- function(spectre_object,
                                x_eigen = NULL,
                                y_eigen = NULL,
                                plot_title = NULL,
                                point_size = 2,
                                point_color = 'steelblue',
                                point_alpha = 0.75,
                                point_wjitter = 0,
                                point_hjitter = 0,
                                cust_theme = ggplot2::theme_classic()) {


    ## entry control ------

    stopifnot(is_spectre(spectre_object))

    if(!inherits(cust_theme, 'theme')) {

      stop("'cust_theme' has to be a valid ggplot theme object.",
           call. = FALSE)

    }

    n <- ncol(spectre_object$eigen_vectors)

    if(is.null(x_eigen)) x_eigen <- n

    if(is.null(y_eigen)) y_eigen <- n - 1

    stopifnot(is.numeric(x_eigen))
    stopifnot(is.numeric(y_eigen))

    if(x_eigen > n | y_eigen > n) {

      error_txt <- paste("Maximal allowed number of eigenvectors is", n)

      stop(error_txt, call. = FALSE)

    }

    ## plotting data and labels -------

    if(is.null(plot_title)) plot_title <- 'Selected eigenvectors'

    plot_subtitle <- paste('observations: n =', n)

    x_lab <- paste('Eigenvector', x_eigen)

    y_lab  <- paste('Eigenvector', y_eigen)

    plot_data <-
      as.data.frame(spectre_object$eigen_vectors[, c(x_eigen, y_eigen)])

    plot_data <- set_names(plot_data, c('comp_1', 'comp_2'))

    ## plotting --------

    comp_1 <- NULL
    comp_2 <- NULL

    ggplot(plot_data,
           aes(x = comp_1,
               y = comp_2)) +
      geom_point(shape = 21,
                 size = point_size,
                 fill = point_color,
                 alpha = point_alpha,
                 position = position_jitter(width = point_wjitter,
                                            height = point_hjitter)) +
      cust_theme +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab)

  }

# Laplacian Spectralization ---------

#' Laplacian decomposition.
#'
#' @description
#' Performs a spectral decomposition of a pairwise distances between
#' observations of a data set.
#'
#' @details
#' The code is inspired by two sources:
#' https://rpubs.com/gargeejagtap/SpectralClustering
#' http://www.di.fc.ul.pt/~jpn/r/spectralclustering/spectralclustering.html.
#'
#' The function combines the following computation steps:
#'
#' * _distance matrix calculation_. The distance matrix is calculated via
#' \code{\link{calculate_dist}}, you may check available distance types by
#' calling \code{\link{get_kernel_info}}.
#'
#' * _calculation of affinity matrix_. This step is accomplished by
#' \code{\link{dist2affi}} and results in a numeric representation of the
#' nearest neighborhood, i.e. affinity matrix A. Such matrix may be unweighted
#' (nearest neighbors get `1`, all other pairs 0) or weighted by similarity
#' (nearest neighbors: similarity statistic, all other pairs 0). It is important
#' to note, that conversion of the distance matrix to similarity matrix is done
#' by applying a user-provided function `simil_fun`. Its default value
#' corresponds to a common transformation of Euclidean distance to Euclidean
#' similarity, which must not be optimal for all distance measures.
#'
#' * _calculation of degree matrix and Laplacian_. The degree matrix D is a
#' diagonal matrix, whose diagonal are sums of columns of the affinity matrix A.
#' As such, D stores the numbers of neighbors of a given observation or, from
#' a graph perspective, degree of the graph nodes. The Laplacian matrix U is
#' computed as a simple difference: \eqn{U = D - A}. It is, optionally,
#' normalized with the following formula: \eqn{U_{norm} = D^{-1/2} U D^{-1/2}}.
#'
#' * _decomposition of the Laplacian_. This step is done via base R's
#' \code{\link[base]{eigen}} and generates eigenvectors and eigenvalues of the
#' Laplacian matrix U. The trailing smallest eigenvectors with eigenvalues
#' not departing substantially from 0 represent the base dimensions of the data
#' set and may be used for further analysis steps such as clustering.
#' The `plot()` method called for the `spectralize()` function output may be
#' helpful ad finding such trailing eigenvectors.
#'
#' @return a list of class \code{\link{spectre}} with the following elements:
#'
#' * `degrees`: degrees of the graph nodes stored in the diagonal of the degree
#' matrix
#'
#' * `eigen_values`: a numeric vector of eigenvalues sorted from the largest to
#' the smallest one
#'
#' * `eigen_vectors`: a numeric matrix whose rows represent the observations in
#' he initial data set and columns representing the eigenvectors.
#' The eigenvectors are sorted by their einegvalues: the largest come first.
#'
#' Optionally, if `return_laplacian == TRUE`, the Laplacian matrix
#' (with or without normalization) is returned as well.
#'
#' @param data a numeric data frame or matrix.
#' @param distance_method distance method. Call `get_kernel_info` for
#' available distance metrics.
#' @param kNN numeric, number of the nearest neighbors., has to be lower than
#' the dimension of `dist_mtx`.
#' @param weighted logical, should the affinity matrix be weighted by
#' similarity? Defaults to `TRUE`.
#' @param simil_fun a function used to convert pairwise distances
#' to pairwise similarities.
#' @param norm_laplacian logical, should the Laplacian matrix be normalized
#' prior to decomposition?
#' @param return_laplacian logical, should the Laplacian matrix be included in
#' the function output?
#'
#' @export

  spectralize <- function(data,
                          distance_method = 'euclidean',
                          kNN = 5,
                          weighted = FALSE,
                          simil_fun = function(x) 1/(1 + x),
                          norm_laplacian = TRUE,
                          return_laplacian = FALSE) {

    ## inspired by
    ## https://rpubs.com/gargeejagtap/SpectralClustering
    ## http://www.di.fc.ul.pt/~jpn/r/spectralclustering/spectralclustering.html

    ## entry control: in part done by downstream functions -------

    check_numeric(data)

    stopifnot(is.logical(weighted))
    stopifnot(is.logical(norm_laplacian))

    ## computation of the distance, affinity, and degree matrix --------

    dist_mtx <- calculate_dist(data = data,
                               method = distance_method)

    n <- nrow(data)

    ## affinity matrices for weighted and unweighted cases

    A <- dist2affi(dist_mtx,
                   kNN = kNN,
                   weighted = weighted,
                   simil_fun = simil_fun)

    ## degrees of the nodes, i.e. the numbers of neighbors
    ## or weight sums per node

    degrees <- colSums(A)

    D <- diag(degrees)

    dimnames(D) <- dimnames(dist_mtx)

    ## Laplacian and its decomposition ------

    U <- D - A

    if(norm_laplacian) {

      U <- (D %^% -0.5) %*% U %*% (D %^% -0.5)

    }

    eigen_object <- eigen(U, symmetric = TRUE)

    ## output ---------

    comp_names <- paste0('eigen_', 1:n)

    vectors <- eigen_object$vectors

    dimnames(vectors) <- list(rownames(data), comp_names)

    values <- set_names(eigen_object$values, comp_names)

    out <- list(degrees = degrees,
                eigen_values = values,
                eigen_vectors = vectors)

    if(return_laplacian) {

      out <- c(out,
               list(laplacian = U))

    }

    spectre(out)

  }

# END -----
