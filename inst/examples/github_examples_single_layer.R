# Git Hub examples of a clustering analysis of single-layer data.
# Working with the biopsy data set from the MASS package

# libraries and data sources -------

  ## data source, the loading order matters
  ## load MASS as the first package

  library(MASS)

  library(tidyverse)
  library(clustTools)
  library(somKernels)

  library(patchwork)

  ## the MASS biopsy data set.
  ## Including only complete cases in the analysis
  ## The ID is not unique, since multiple samples
  ## are available for the same patient, so I'm redefining it
  ## no scaling of the clustering features V1 - V9, they are
  ## on the same scale

  my_biopsy <- biopsy %>%
    filter(complete.cases(.)) %>%
    mutate(observation = paste0('obs_', 1:nrow(.)),
           ID = paste(observation, ID, sep = '_')) %>%
    as_tibble

  ## training portion: a random subset of 383 observations
  ## the rest will be used for validation ('hold-out')

  set.seed(1234)

  train_ids <- sample(1:nrow(my_biopsy),
                      size = 383,
                      replace = FALSE)

  biopsy_lst <-
    list(train = my_biopsy[train_ids, ],
         test = my_biopsy[-train_ids, ]) %>%
    map(column_to_rownames, 'ID') %>%
    map(select, starts_with('V'))

# Clustering tendency of the training and test data subsets ------

  biopsy_clust_tendency <- biopsy_lst %>%
    map(get_clust_tendency,
        n = 200)

# Clustering analysis objects -------

  clust_objects <- list()

  ## Ward's hierarchical clustering

  clust_objects$hcl_euclidean <-
    hcluster(data = biopsy_lst$train,
             distance_method = 'euclidean',
             k = 2,
             hc_method = 'ward.D2')

  clust_objects$hcl_cosine <-
    hcluster(data = biopsy_lst$train,
             distance_method = 'cosine',
             k = 2,
             hc_method = 'ward.D2')

  ## kmeans

  clust_objects$kmeans_euclidean <-
    kcluster(data = biopsy_lst$train,
             distance_method = 'euclidean',
             k = 2,
             clust_fun = 'kmeans')

  clust_objects$kmeans_cosine <-
    kcluster(data = biopsy_lst$train,
             distance_method = 'cosine',
             k = 2,
             clust_fun = 'kmeans')

  ## PAM

  clust_objects$pam_euclidean <-
    kcluster(data = biopsy_lst$train,
             distance_method = 'euclidean',
             k = 2,
             clust_fun = 'pam')

  clust_objects$pam_cosine <-
    kcluster(data = biopsy_lst$train,
             distance_method = 'cosine',
             k = 2,
             clust_fun = 'pam')

# Cluster number choice --------

  ## plot(type = 'diagnostic') returns
  ## a list of diagnostic plots for determination
  ## of the optimal cluster number

  cl_number_plots <- clust_objects %>%
    map(plot,
        type = 'diagnostic') %>%
    map(map,
        ~.x +
          theme(plot.tag.position = 'bottom'))

# Silhouettes Clustering variance --------

  ## mean silhouettes

  cl_silhouettes <- clust_objects %>%
    map(silhouette) %>%
    map_dfr(summary) %>%
    filter(clust_id == 'global') %>%
    mutate(clustering_algorithm = names(clust_objects)) %>%
    relocate(clustering_algorithm)

  ## explained clustering variances

  cl_variances <- clust_objects %>%
    map(var) %>%
    map_dbl(~.x$frac_var)

# Heat maps of distances ------

  ## distance heat maps

  pam_heat_maps <- clust_objects[c("pam_euclidean", "pam_cosine")] %>%
    map(plot, type = 'heat_map') %>%
    map(~.x +
          theme(plot.tag = element_blank(),
                axis.text = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank(),
                legend.position = 'bottom'))

  ## heat maps of mean distances between the clusters

  pam_mean_heat_maps <- clust_objects[c("pam_euclidean", "pam_cosine")] %>%
    map(cross_distance) %>%
    map(plot, type = 'mean') %>%
    map(~.x + theme(legend.position = 'bottom'))

# Cross-validation ------

  euclidean_cv_stats <-
    clust_objects[c("hcl_euclidean", "kmeans_euclidean", "pam_euclidean")] %>%
    map(cv) %>%
    map(summary)

# Semi-supervised clustering ------

  ## predictions for KMEANS/Euclidean

  kmeans_clusters <- list()

  kmeans_clusters$training <- clust_objects$kmeans_euclidean

  kmeans_clusters$test <-
    predict(clust_objects$kmeans_euclidean,
            newdata = biopsy_lst$test,
            type = 'propagation')

  ## explained variance and silhouette width

  kmeans_variance <- kmeans_clusters %>%
    map(var) %>%
    map_dbl(~.x$frac_var)

  kmeans_silhouette <- kmeans_clusters %>%
    map(silhouette) %>%
    map(summary)

  ## heat-maps

  kmeans_heat_maps <- kmeans_clusters %>%
    map(plot, type = 'heat_map') %>%
    map(~.x +
          theme(legend.position = 'bottom',
                plot.tag = element_blank(),
                axis.text = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank())) %>%
    map2(., c('Biopsy: training', 'Biopsy: test'),
         ~.x + labs(title = .y))

  ## cross-distances

  kmeans_cross_dists <-
    cross_distance(kmeans_clusters$training,
                   kmeans_clusters$test)

# Benign and malignant samples in the training and test clusters -----

  ## retrieving the cluster assignments

  kmeans_clust_assign <- kmeans_clusters %>%
    map(extract, 'assignment') %>%
    map(mutate, ID = observation)

  ## joining with the biopsy data

  kmeans_clust_assign <- kmeans_clust_assign %>%
    map(left_join,
        my_biopsy[c('ID', 'class')],
        by = 'ID')

  ## counting benign and malignant samples in the clusters

  kmeans_clust_counts <- kmeans_clust_assign %>%
    map(count, clust_id, class)

  ## kappa and ROC: I'm using `multiClassSummary()`
  ## which takes a data frame with the `obs` and `pred`
  ## variables storing the observed and predicted outcome, respectively

  library(caret)

  kmeans_clust_assign <- kmeans_clust_assign %>%
    map(mutate,
        obs = class,
        pred = car::recode(clust_id,
                           "'1' = 'benign'; '2' = 'malignant'"),
        pred = factor(pred, c('benign', 'malignant')))

  kmeans_roc_stats <- kmeans_clust_assign %>%
    map(select, obs, pred) %>%
    map(as.data.frame) %>%
    map(multiClassSummary,
        lev = c('benign', 'malignant'))

# Advanced visualization options ------

  ## plots of MDS of the distance matrix, as an alternative
  ## of distance heat map

  kmeans_mds_dist <- kmeans_clusters %>%
    map(plot,
        type = 'components',
        red_fun = 'mds') %>%
    map(~.x +
          theme(plot.tag = element_blank(),
                legend.position = 'bottom'))

  ## UMAP of the data set

  kmeans_umap_data <- kmeans_clusters %>%
    map(plot,
        type = 'components',
        red_fun = 'umap',
        with = 'data') %>%
    map(~.x +
          theme(plot.tag = element_blank(),
                legend.position = 'bottom'))

  ## heat maps of clustering features

  kmeans_hm_variables <- kmeans_clusters %>%
    map(plot_clust_hm) %>%
    map(~.x +
          scale_fill_gradient2(low = 'steelblue',
                               mid = 'black',
                               high = 'firebrick',
                               midpoint = 5.5,
                               limits = c(1, 10)) +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.line.x = element_blank(),
                legend.position = 'bottom')) %>%
    map2(., c('Biopsy: training', 'Biopsy: test'),
         ~.x + labs(title = .y))

# Clustering variable importance -------

  kmeans_variable_importance <-
    impact(kmeans_clusters$training,
           n_iter = 50,
           .parallel = TRUE)

# Density clustering ------

  ## generating the UMAP layouts for the training and test
  ## portion of the data set

  biopsy_density <- biopsy_lst %>%
    map(reduce_data,
        distance_method = 'manhattan',
        kdim = 2,
        red_fun = 'umap',
        random_state = 1234)

  ## clustering the training subset
  ## eps is a pure guess to begin with

  biopsy_dbscan <- list()

  biopsy_dbscan$training <-
    dbscan_cluster(data = biopsy_density$train,
                   distance_method = 'manhattan',
                   eps = 1,
                   minPts = 7)


  biopsy_dbscan$training <-
    dbscan_cluster(data = biopsy_density$train,
                   distance_method = 'manhattan',
                   eps = 0.7,
                   minPts = 7)

  ## predictions

  biopsy_dbscan$test <-
    predict(biopsy_dbscan$training,
            newdata = biopsy_density$test,
            type = 'propagation')

  ## clustering variance and silhouette width

  biopsy_dbscan_variance <- biopsy_dbscan %>%
    map(var) %>%
    map_dbl(~.x$frac_var)

  biopsy_dbscan_silhouette <- biopsy_dbscan %>%
    map(silhouette) %>%
    map(summary)

  ## UMAP plots

  biopsy_dbscan_umap_layout <- biopsy_dbscan %>%
    map(plot,
        type = 'data') %>%
    map(~.x +
          theme(plot.tag = element_blank(),
                legend.position = 'bottom'))

# END ------
