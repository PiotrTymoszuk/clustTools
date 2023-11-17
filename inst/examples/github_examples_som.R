# Single layer SOM clustering with the wine dataset

# tools and data -------

  ## kohonen will be used for diagnostic plots
  ## and as the data source
  ##
  ## the library order matters!

  library(kohonen)

  library(tidyverse)
  library(clustTools)
  library(somKernels)
  library(caret)

  ## patchwork for stitching the plots

  library(patchwork)

  ## the wines dataset pre-processing: elimination of invariant features
  ## (low variance to mean ratio) and median centering

  data("wines")

  my_wines <- wines %>%
    as.data.frame %>%
    mutate(ID = paste0('wine_', 1:nrow(.))) %>%
    column_to_rownames('ID')

  distr_stats <- my_wines %>%
    map_dfr(~tibble(variance = var(.x),
                    mean = mean(.x),
                    var_mean_ratio = var(.x)/mean(.x))) %>%
    mutate(variable = names(my_wines)) %>%
    relocate(variable)

  top_variant_features <- distr_stats %>%
    filter(var_mean_ratio > 0.1) %>%
    .$variable

  my_wines <- my_wines %>%
    select(all_of(top_variant_features)) %>%
    center_data('median')

  ## appending the parental data frame with wine classification
  ## it will be used for the final validation of the results

  my_wines <- my_wines %>%
    mutate(vintage = vintages)

  ## training: 100 randomly chosen observations
  ## the rest used for validation#
  ## created with caret's createDataPartition() to keep
  ## the vintage distribution

  set.seed(12345)

  train_ids <- createDataPartition(my_wines$vintage, p = 100/177)[[1]]

  wine_lst <-
    list(train = my_wines[train_ids, ],
         test = my_wines[-train_ids, ]) %>%
    map(select, -vintage)

# Clustering tendency ------

  clust_tendency <- wine_lst %>%
    map(get_clust_tendency,
        n = 60)

# Tuning of the clustering structure --------

  algos <- list()

  algos$som_pam_euclidean <-
    combi_cluster(data = wine_lst$train,
                  distance_som = 'euclidean',
                  xdim = 5,
                  ydim = 4,
                  topo = 'hexagonal',
                  neighbourhood.fct = 'gaussian',
                  toroidal = FALSE,
                  rlen = 3000,
                  node_clust_fun = kcluster,
                  k = 3,
                  clust_fun = 'pam')

  algos$som_pam_manhattan <-
    combi_cluster(data = wine_lst$train,
                  distance_som = 'manhattan',
                  xdim = 5,
                  ydim = 4,
                  topo = 'hexagonal',
                  neighbourhood.fct = 'gaussian',
                  toroidal = FALSE,
                  rlen = 3000,
                  node_clust_fun = kcluster,
                  k = 3,
                  clust_fun = 'pam')

  algos$som_pam_cosine <-
    combi_cluster(data = wine_lst$train,
                  distance_som = 'cosine',
                  xdim = 5,
                  ydim = 4,
                  topo = 'hexagonal',
                  neighbourhood.fct = 'gaussian',
                  toroidal = FALSE,
                  rlen = 3000,
                  node_clust_fun = kcluster,
                  k = 3,
                  clust_fun = 'pam')

# Checking the SOM convergence ------

  som_convergence_plots <- algos %>%
    map(plot, type = 'training') %>%
    map(~.x$observation +
          theme(plot.tag = element_blank(),
                legend.position = 'none')) %>%
    map2(., c('SOM/HCL, Euclidean distance',
              'SOM/HCL, Manhattan distance',
              'SOM/HCL, cosine distance'),
         ~.x +
           labs(title = .y,
                subtitle = 'SOM training convergence'))

  som_convergence_plots$som_pam_euclidean +
    som_convergence_plots$som_pam_manhattan +
    som_convergence_plots$som_pam_cosine +
    plot_layout(ncol = 2)

  algos %>%
    map_dbl(qe)

  algos %>%
    map(pbc) %>%
    map_dbl(~.x$convergence_stat)

  ## topology errors within SOM nodes and the final clusters

  algos %>%
    map(te, type = 'node') %>%
    map(summary) %>%
    map(filter,  clust_id == 'global')

  algos %>%
    map(te, type = 'final') %>%
    map(summary)

  ## neighborhood preservation:
  ## data -> SOM nodes
  ## SOM nodes -> clusters
  ## data -> clusters

  algos %>%
    map(np, type = 'data') %>%
    map(summary) %>%
    map(filter, clust_id == 'global')

  algos %>%
    map(np, type = 'node') %>%
    map(summary)

  algos %>%
    map(np, type = 'final') %>%
    map(summary)

# Clustering of the U matrix, cluster number ------

  node_cluster_number_plots <- algos %>%
    map(plot, type = 'diagnostic') %>%
    map(~.x$node$silhouette +
          theme(plot.tag = element_blank())) %>%
    map2(., c('SOM/HCL, Euclidean distance',
              'SOM/HCL, Manhattan distance',
              'SOM/HCL, cosine distance'),
         ~.x +
           scale_y_continuous(limits = c(0, 0.65)) +
           labs(title = .y,
                subtitle = 'SOM node clustering'))

  node_cluster_number_plots$som_pam_euclidean +
    node_cluster_number_plots$som_pam_manhattan +
    node_cluster_number_plots$som_pam_cosine +
    plot_layout(ncol = 2)

# Numeric stats of the algorithms --------

  algos_variance <- algos %>%
    map(var) %>%
    map_dbl(~.x$frac_var)

  algos_silhouette <- algos %>%
    map(silhouette) %>%
    map(summary)

# Cross-validation -------

  ## SOM prediction method is recommended

  algos_cv <- algos %>%
    map(cv, type = 'som') %>%
    map(summary) %>%
    map(select, ends_with('mean'))

# Semi-supervised clustering ------

  cosine_clusters <- list()

  cosine_clusters$train <- algos$som_pam_cosine

  cosine_clusters$test <- predict(cosine_clusters$train,
                                  newdata = wine_lst$test,
                                  type = 'som')

  ## comparison of variances and silhouette widths
  ## in the training and test data portions

  cosine_variance <- cosine_clusters %>%
    map(var) %>%
    map_dbl(~.x$frac_var)

  cosine_silhouettes <- cosine_clusters %>%
    map(silhouette) %>%
    map(summary)

  ## neighborhood preservation

  cosine_neighborhood <- list()

  cosine_neighborhood$train <- np(cosine_clusters$train,
                                  kNN_data = 5,
                                  type = 'final')

  cosine_neighborhood$test <- np(cosine_clusters$test,
                                 kNN_data = 5)

  cosine_neighborhood_plots <- cosine_neighborhood %>%
    map(plot) %>%
    map2(., c('Wines: training', 'Wines: test'),
         ~.x + labs(title = .y))

  cosine_neighborhood_plots$train +
    cosine_neighborhood_plots$test

# Visualization: heat maps -------

  ## U-matrix plots are available only for the training data
  ## and are stored in the 'node' element of the output plot list

  cosine_umatrix_hm <- plot(cosine_clusters$train,
                            type = 'heat_map')$node +
    labs(title = 'SOM/PAM, cosine distance',
         subtitle = 'U-matrix clustering, wines training subset') +
    theme(axis.text = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.tag = element_blank(),
          legend.position = 'bottom')

  ## test subset: heat map of pairwise distances between observations

  cosine_test_hm <- plot(cosine_clusters$test,
                         type = 'heat_map') +
    labs(title = 'Predictions',
         subtitle = 'Distances between observations, wines test subset') +
    theme(axis.text = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.tag = element_blank(),
          legend.position = 'bottom')

  cosine_umatrix_hm +
    cosine_test_hm

# Visualization: UMAP -------

  cosine_umap_plots <- cosine_clusters %>%
    map(plot,
        type = 'components',
        kdim = 2,
        with = 'data',
        red_fun = 'umap')

  ## for the training data subset:
  ## a list of UMAP layout plots is returned, we need only
  ## only the UMAP layout for the observations with color coding
  ## of the cluster assignment

  cosine_umap_plots$train <- cosine_umap_plots$train$final

  cosine_umap_plots <-
    map2(cosine_umap_plots,
         c('SOM/PAM, cosine distance, training',
           'SOM/PAM, cosine distance, test'),
         ~.x +
           labs(title = .y) +
           theme(plot.subtitle = element_blank(),
                 plot.tag = element_blank(),
                 legend.position = 'bottom'))

  cosine_umap_plots$train +
    cosine_umap_plots$test

# Visualization: heat map of the clustering features ------

  cosine_feature_hm <- cosine_clusters %>%
    map(plot_clust_hm) %>%
    map2(., c('SOM/PAM, cosine distance, training',
              'SOM/PAM, cosine distance, test'),
         ~.x +
           labs(title = .y) +
           scale_fill_gradient2(low = 'steelblue',
                                mid = 'black',
                                high = 'firebrick',
                                midpoint = 0,
                                limits = c(-3, 3),
                                oob = scales::squish) +
           theme(axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.line.x = element_blank(),
                 legend.position = 'bottom'))

  cosine_feature_hm$train +
    cosine_feature_hm$test

# Cross distances -----

  cross_distance(cosine_clusters$train,
                 cosine_clusters$test) %>%
    plot('mean')

# Feature importance -------

  cosine_importance <- impact(cosine_clusters$train,
                              n_iter = 50,
                              .parallel = TRUE)

  plot(cosine_importance)

# Vintage classes in the clusters --------

  ## assignment extraction works the same
  ## way as for non-SOM analyses

  vintage_assignment <- my_wines %>%
    rownames_to_column('observation') %>%
    select(observation, vintage)

  cosine_assignment <- cosine_clusters %>%
    map(extract, 'assignment') %>%
    map(left_join, vintage_assignment, by = 'observation')

  cosine_counts <- cosine_assignment %>%
    map(count, clust_id, vintage)

  ## kappa and ROC analysis
  ## renaming of the clusters after
  ## the predominant wine type

  cosine_assignment <- cosine_assignment %>%
    map(mutate,
        obs = vintage,
        pred = car::recode(clust_id,
                           "'1' = 'Barbera';
                           '2' = 'Barolo';
                           '3' = 'Grignolino'"),
        pred = factor(pred,
                      levels = c('Barbera', 'Barolo', 'Grignolino')))

  cosine_roc <- cosine_assignment %>%
    map(as.data.frame) %>%
    map(multiClassSummary,
        lev = c('Barbera', 'Barolo', 'Grignolino'))

# END -----
