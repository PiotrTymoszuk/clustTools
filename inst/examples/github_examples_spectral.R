# Example of spectral clustering

# tools and analysis data -------

  library(kernlab)
  library(tidyverse)
  library(clustTools)

  extract <- clustTools::extract

  ## the spiral data set: amplified

  data("spirals")

  colnames(spirals) <- c('dim1', 'dim2')

  spirals <- as_tibble(spirals)

# traditional clustering solutions: hierarchical and K-MEANS clustering -------

  kclust_object <- kcluster(spirals,
                            distance_method = 'euclidean',
                            clust_fun = 'kmeans',
                            k = 2)

  hclust_object <- hcluster(spirals,
                            distance_method = 'euclidean',
                            k = 2,
                            hc_method = 'single')

  ## you may find the optimal `eps` value by calling `plot()` for the output
  ## of the dbscan_cluster

  dbscan_object <- dbscan_cluster(spirals,
                                  distance_method = 'euclidean',
                                  eps = 0.2)

# spectral clustering: matrix decomponsition -------

  ## decomposition and checking the data set base (the trailing eigenvectors)

  spectre_object <-
    spectralize(data = spirals,
                distance_method = 'euclidean',
                kNN = 4,
                weighted = FALSE,
                #simil_fun = function(x) exp(-3 * x),
                norm_laplacian = FALSE)

  plot(spectre_object,
       type = 'eigenvalues')

  plot(spectre_object,
       type = 'eigenvectors',
       point_hjitter = 0.001,
       point_wjitter = 0.001)

  ## `red_analysis` object, which will be used later in clustering

  red_analysis_object <- as_red_analysis(spectre_object,
                                         kdim = 2,
                                         skip_last = FALSE)

# Spectral clustering: K-MEANS -------

  kspectre_object <- red_analysis_object %>%
    kcluster(k = 2)

  ## diagnostic plots and numeric stats of cluster performance:
  ## k = 2 clusters is obviously the best solution!

  kspectre_object %>%
    plot

  kspectre_object %>%
    summary

  kspectre_object %>%
    ngroups

# comparison of the clustering solutions --------

  ## appending the `spirals` data set with the cluster labels

  spirals <- spirals %>%
    mutate(kclust_labs = extract(kclust_object, 'assignment')$clust_id,
           hclust_labs = extract(hclust_object, 'assignment')$clust_id,
           dbscan_labs = extract(dbscan_object, 'assignment')$clust_id,
           kspectre_labs = extract(kspectre_object, 'assignment')$clust_id)

  ## plots

  spiral_plots <-
    list(x = c('kclust_labs',
               'hclust_labs',
               'dbscan_labs',
               'kspectre_labs'),
         y = c('K-MEANS',
               'Single-linkage HCL',
               'DBSCAN',
               'Spectral K-MEANS')) %>%
    pmap(function(x, y) spirals %>%
           ggplot(aes(x = dim1,
                      y = dim2,
                      fill = .data[[x]])) +
           geom_point(shape = 21,
                      size = 2) +
           theme_light() +
           labs(title = y,
                fill = 'Cluster')) %>%
    set_names(c('kclust', 'hclust', 'dbscan', 'kspectre'))

# END ------
