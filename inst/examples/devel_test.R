# Tests the functionality during development

  library(clustTools)
  library(somKernels)
  library(MASS)
  library(tidyverse)

  select <- dplyr::select
  extract <- clustTools::extract

# globals -----

  test_data <- iris[1:4]

  test_data <- set_rownames(test_data,
                            paste0('obs_', 1:nrow(test_data)))

  test_data <- center_data(test_data,
                           type = 'median',
                           complete_cases = TRUE)

  test_data <- min_max(test_data)

  clustTools:::check_numeric(test_data)

  set.seed(1234)

  new_ids <- sample(1:nrow(test_data),
                    size = 50,
                    replace = FALSE)

  new_data <- test_data[new_ids, ]

  test_data <- test_data[-new_ids, ]

# testing the utils ------

  test_distances <- purrr::map(clustTools::get_kernel_info(),
                               calculate_dist,
                               data = test_data)

  test_distances <- rlang::set_names(test_distances,
                                     clustTools::get_kernel_info())

  clustTools:::plot_knn_distance(diss_obj = as.dist(test_distances$dice),
                                 k = 5,
                                 eps = 0.006)

  test_hcl <- hclust(d = as.dist(test_distances$dice))

  test_kohonen <- kohonen::som(as.matrix(test_data),
                               grid = kohonen::somgrid(xdim = 5, ydim = 4),
                               rlen = 1000,
                               dist.fcts = 'jaccard')

  #plot_dendro(test_hcl$clust_obj,
   #           k = 3,
    #          labels = FALSE,
     #         cluster_colors = c('firebrick',
      #                           'steelblue',
       #                          'darkolivegreen',
        #                         'black'))

  plot_nbclust(test_distances$sumofsquares,
               k = 3,
               FUNcluster = cluster::pam,
               method = 'silhouette')

  plot_som(test_kohonen)

  plot_train_som(test_kohonen)

  plot_point(data = iris,
             x_var = 'Sepal.Length',
             y_var = 'Sepal.Width',
             fill_var = 'Species')

  clustTools:::get_data_dim(test_distances$cosine)

  test_vec <- c(rep('a', 10), rep('b', 5), rep('c', 10))

  test_dists <- c(sample(5:15, 10), sample(1:10, 5), sample(1:10, 10))

  clustTools:::vote_simple(test_vec, resolve_ties = FALSE)

  clustTools:::vote_kernel(test_vec,
                           dist_vec = test_dists)


# testing the reduction methods -----

  ## reduction analysis objects

  test_red <- purrr::map(c('mds', 'pca', 'umap'),
                         reduce_data,
                         data = test_data,
                         distance_method = 'cosine',
                         kdim = 3)

  test_red <- rlang::set_names(test_red,
                               c('mds', 'pca', 'umap'))

  test_red$fa <-
    reduce_data(dplyr::mutate(MASS::painters[1:4],
                              dummy1 = sample(1:20,
                                              nrow(MASS::painters),
                                              replace = TRUE),
                              dummy2 = sample(1:20,
                                              nrow(MASS::painters),
                                              replace = TRUE)),
                red_fun = 'fa',
                kdim = 2,
                scores = 'Bartlett')

  ## their methods

  extract(test_red$pca, 'data')

  model.frame(test_red$umap)

  nobs(test_red$mds)

  var(test_red$umap)

  var(test_red$fa)
  nobs(test_red$fa)
  model.frame(test_red$fa)

  ## plotting

  plot(test_red$mds, type = 'scores')

  plot(test_red$umap, type = 'scores')

  plot(test_red$umap, type = 'scree')

  plot(test_red$pca,
       type = 'loadings',
       segment_color = 'gray60',
       point_color = 'coral2')

  plot(test_red$fa,
       type = 'loadings',
       segment_color = 'gray60',
       point_color = 'coral2')

  plot(test_red$fa,
       type = 'scores',
       segment_color = 'gray60',
       point_color = 'coral2')

  plot(test_red$fa,
       type = 'scree')

# Clustering tendency -----

  get_clust_tendency(test_data, n = 50)

# simple clustering ------

  test_hcl <- hcluster(data = test_red$umap,
                       distance_method = 'euclidean',
                       hc_method = 'complete',
                       k = 3)

  test_kmeans <- kcluster(data = test_red$pca,
                          distance_method = 'euclidean',
                          clust_fun = 'kmeans',
                          k = 3)

  test_pam <- kcluster(data = test_red$umap,
                       distance_method = 'manhattan',
                       clust_fun = 'pam',
                       k = 3)

  test_dbscan <- dbscan_cluster(data = test_red$umap,
                                distance_method = 'euclidean',
                                eps = 1,
                                minPts = 10)

  test_som <- som_cluster(data = test_data,
                          distance_method = 'cosine',
                          xdim = 4,
                          ydim = 4,
                          topo = 'hexagonal',
                          rlen = 1000)

# Clust_analysis OOP -----

  extract(test_som, 'assignment')

  model.frame(test_pam)

  dist(test_hcl)

  nobs(test_kmeans)

  ngroups(test_dbscan)

  components(test_som,
             kdim = 2,
             red_fun = 'pca',
             with = 'data')

# Clust_analysis plotting ------

  plot(test_pam, type = 'diagnostic')

  plot(test_dbscan, type = 'diagnostic')

  plot(test_som, type = 'diagnostic')

  plot(test_som,
       type = 'components',
       with = 'data',
       red_fun = 'umap',
       kdim = 2)

  plot(test_pam,
       type = 'components',
       with = 'data',
       red_fun = 'pca',
       kdim = 2)

  plot(test_dbscan,
       type = 'components',
       with = 'data',
       red_fun = 'pca',
       kdim = 2)

  plot(test_dbscan, type = 'heat_map')

  plot(test_som, type = 'training')

  plot(test_dbscan, type = 'data')

# Semi-supervised clustering ----

  ## predictions for the SOM clusters

  test_pred <- predict(object = test_som,
                       newdata = new_data,
                       type = 'propagation',
                       simple_vote = FALSE)

  test_neuro_pred <- predict(object = test_som,
                             newdata = new_data,
                             type = 'som')

  test_pred %>%
    plot('components',
         with = 'data',
         red_fun = 'umap',
         kdim = 2)

  test_neuro_pred %>%
    plot('components',
         with = 'data',
         red_fun = 'umap',
         kdim = 2)

  test_som %>%
    plot('components',
         with = 'data',
         red_fun = 'umap',
         kdim = 2)

  ## predictions for reduction analysis objects

  test_pred_red <- predict(object = test_pam,
                           newdata = reduce_data(data = new_data,
                                                 distance_method = 'cosine',
                                                 kdim = 3,
                                                 red_fun = 'umap'),
                           type = 'propagation')

  plot(test_pred_red, type = 'data')

  plot(test_pam, type = 'data')

  test_pred_pca <- predict(object = test_dbscan,
                           newdata = reduce_data(data = new_data,
                                                 distance_method = 'cosine',
                                                 kdim = 3,
                                                 red_fun = 'pca'),
                           type = 'propagation')

  plot(test_pred_pca, type = 'components', kdim = 2)

  plot(test_dbscan, type = 'components')

# Cross-validation ------

  cv(test_hcl,
     nfolds = 10,
     kNN = 5,
     type = 'propagation')

  cv(test_pam,
     nfolds = 10,
     kNN = 5)

  cv(test_som,
     nfolds = 10,
     kNN = 5)$summary

  cv(test_som,
     nfolds = 10,
     type = 'som')$summary

  cv(test_dbscan,
     nfolds = 10,
     kNN = 5)

# Variable importance and its OOP -----

  impact(test_dbscan)

  impact(test_som)

  plot(impact(test_hcl), 'scatter')

# SOM combi clustering and its OOP -----

  test_combi <- combi_cluster(data = test_data,
                              distance_som = 'manhattan',
                              xdim = 5,
                              ydim = 4,
                              topo = 'hexagonal',
                              neighbourhood.fct = 'gaussian',
                              toroidal = FALSE,
                              rlen = 1000,
                              node_clust_fun = hcluster,
                              k = 3,
                              hc_method = 'complete')

  extract(test_combi, 'clust_object')

  is_combi_analysis(test_combi)

  dist(test_combi)

  components(test_combi, kdim = 3, red_fun = 'umap', with = 'data')

  dist(test_combi, type = 'umatrix')

  components(test_combi, kdim = 2, red_fun = 'mds', with = 'umatrix') %>%
    plot

  components(test_combi, kdim = 2, red_fun = 'mds', with = 'distance') %>%
    plot

  var(test_combi)

  cv(test_combi, type = 'som')
  cv(test_combi, type = 'propagation')

# SOM combi clustering plots -------

  plot(test_combi, type = 'diagnostic')

  plot(test_combi, type = 'training')

  plot(test_combi,
       type = 'components',
       red_fun = 'umap',
       with = 'data',
       kdim = 2)

  plot(test_combi,
       type = 'data')

# SOM combi semi-supervised clustering ------

  test_new_combi <- predict(object = test_combi,
                            newdata = new_data,
                            type = 'propagation',
                            kNN = 10,
                            resolve_ties = TRUE)

  test_new_neuro_combi <- predict(object = test_combi,
                                  newdata = new_data,
                                  type = 'som')

  plot(test_new_combi,
       type = 'components',
       red_fun = 'umap',
       with = 'data',
       k = 2)

  plot(test_new_neuro_combi,
       type = 'components',
       red_fun = 'umap',
       with = 'data',
       k = 2)

  plot(test_combi,
       type = 'components',
       red_fun = 'umap',
       with = 'data',
       k = 2)

# Feature importance ------

  test_impact <-
    impact(test_pam,
           n_iter = 100,
           .parallel = TRUE)

  test_som_impact <-
    impact(test_combi,
           n_iter = 100,
           .parallel = TRUE)

  plot(test_som_impact)
  summary(test_som_impact)

# Feature plotting ------

  ft_clust <- hcluster(data = t(test_data), distance_method = 'euclidean', k = 2)

  plot(ft_clust)

  plot_clust_hm(x_object = test_combi,
                y_object = ft_clust)

# Cross distance -------

  cross_distance(test_data, new_data)

  cross_distance(new_data, method = 'cosine')

  cross_distance(test_pam)

  test_cross <- cross_distance(test_combi, test_new_combi)

  summary(cross_distance(test_combi, test_new_combi))
  summary(cross_distance(test_combi, test_new_neuro_combi))

  summary(cross_distance(test_pam))

  plot(test_cross)
  plot(test_cross, type = 'mean')
  plot(cross_distance(test_combi, test_new_neuro_combi), type = 'mean')

  plot(cross_distance(test_pam))
  plot(cross_distance(test_pam), type = 'mean')

  plot(test_cross, type = 'hist')
  plot(cross_distance(test_hcl), type = 'hist', color = 'black')

# Silhouettes ------

  test_hcl <-
    rename(test_hcl,
           nm = c('2' = 'cluster_2',
                  '1' = 'cluster_1',
                  '3' = 'cluster_3'))

  test_combi <-
    rename(test_combi,
           nm = c('3' = 'cluster_C',
                  '2' = 'cluster_B',
                  '1' = 'cluster_A'))



  silhouette(test_hcl)

  summary(silhouette(test_combi))

  plot(silhouette(test_combi), fill_by = 'sign')
  plot(silhouette(test_hcl), fill_by = 'neighbor')

# Batch


# Multi-layer SOM -------

  biopsy_data <- MASS::biopsy %>%
    filter(complete.cases(.))

  rownames(biopsy_data) <- paste0('sample_', 1:nrow(biopsy_data))

  biopsy_data <- sample(1:nrow(biopsy_data), size = 383, replace = FALSE) %>%
    list(train = biopsy_data[., ],
         test = biopsy_data[-., ])

  biopsy_data <- biopsy_data[c('train', 'test')]

  biopsy_lst <- biopsy_data %>%
    map(function(dat) list(cell_features = paste0('V', 1:5),
                           nucleus_features = paste0('V', 6:9)) %>%
          map(~dat[, .x]) %>%
          map(as.matrix))

  train_batch <- som_cluster(data = biopsy_lst$train,
                             distance_method = c('manhattan', 'euclidean'),
                             xdim = 8,
                             ydim = 8,
                             topo = 'hexagonal',
                             neighbourhood.fct = 'gaussian',
                             toroidal = FALSE,
                             rlen = 1000,
                             mode = 'batch')

  nobs(train_batch)
  ngroups(train_batch)
  var(train_batch)
  silhouette(train_batch) %>%
    summary

  model.frame(train_batch)

  extract(train_batch, type = 'distance')
  extract(train_batch, type = 'assignment')
  extract(train_batch, type = 'data')
  extract(train_batch, type = 'object')

  dist(train_batch, type = 'umatrix')
  dist(train_batch)

  components(train_batch, kdim = 2, red_fun = 'mds', with = 'umatrix')
  components(train_batch, kdim = 2, red_fun = 'pca', with = 'distance')
  components(train_batch, kdim = 2, red_fun = 'umap', with = 'umatrix')
  components(train_batch, kdim = 2, red_fun = 'umap', with = 'data')

  plot(train_batch, type = 'diagnostic')
  plot(train_batch, type = 'components', with = 'data', kdim = 2, red_fun = 'mds')
  plot(train_batch, type = 'components', with = 'umatrix', kdim = 2, red_fun = 'umap')

  plot(train_batch, 'training')

  train_batch_cv <- cv(train_batch,
                       nfolds = 10,
                       type = 'som',
                       .parallel = TRUE)

  train_batch_importance <- impact(train_batch, n_iter = 10, .parallel = TRUE)

  plot(train_batch_importance)

  test_batch <- predict(train_batch,
                        newdata = biopsy_lst$test,
                        type = 'som')

  nobs(test_batch, 'heat_map')
  ngroups(test_batch)
  var(test_batch)
  silhouette(test_batch) %>%
    summary

  model.frame(train_batch)

  dist(test_batch)

  components(test_batch, kdim = 3, red_fun = 'mds', with = 'distance')
  components(test_batch, kdim = 3, red_fun = 'umap', with = 'distance')
  components(test_batch, kdim = 2, red_fun = 'pca', with = 'distance')

  components(test_batch, kdim = 2, red_fun = 'umap', with = 'data')

  plot(test_batch, type = 'components', with = 'distance', kdim = 2, red_fun = 'umap')

  test_homo_cross <- cross_distance(train_batch,
                                    .parallel = TRUE)

  test_hetero_cross <- cross_distance(train_batch,
                                      test_batch,
                                      .parallel = TRUE)

  plot(test_homo_cross, type = 'mean', show_txt = FALSE, line_color = NA)

  plot(test_hetero_cross, type = 'mean', show_txt = FALSE)

  plot_clust_hm(train_batch)
  plot_clust_hm(test_batch)

# Batch SOM with mtcars -------

  my_cars <- list(numeric = c('mpg', 'disp', 'hp', 'wt', 'qsec'),
                  ordinal = c('cyl', 'gear', 'carb'),
                  binary = c('vs', 'am')) %>%
    map(~mtcars[.x])

  my_cars[c('numeric', 'ordinal')] <-
    my_cars[c('numeric', 'ordinal')] %>%
    map(center_data, 'median')

  my_cars <- my_cars %>%
    map(as.matrix)

  car_batch <-
    som_cluster(data = my_cars,
                distance_method = c('euclidean', 'manhattan', 'tanimoto'),
                xdim = 4,
                ydim = 3,
                topo = 'hexagonal',
                neighbourhood.fct = 'gaussian',
                toroidal = FALSE,
                mode = 'online',
                rlen = 1000)

  nobs(car_batch)
  ngroups(car_batch)
  var(car_batch)
  silhouette(car_batch) %>%
    summary

  silhouette(car_batch) %>%
    plot

  plot(car_batch, 'heat_map')

  car_batch %>%
    predict(newdata = my_cars,
            type = 'som') %>%
    plot('heat_map')

  plot(car_batch,
       type = 'components',
       with = 'umatrix',
       kdim = 2,
       red_fun = 'umap',
       n_neighbors = 9)

  plot(car_batch, 'training') +
    facet_grid(`SOM layer` ~ .,
               scales = 'free')

  batch_cross <- cross_distance(car_batch)

  summary(batch_cross)
  plot(batch_cross)

  cv(car_batch, nfolds = 5, type = 'som')

  perm_car_data <- impact(car_batch, n_iter = 50, .parallel = TRUE)

  plot(perm_car_data)

  plot_clust_hm(car_batch)

# Minimal analysis subclass -------

  ## working with a U matrix of the biopsy clustering object

  biopsy_umatrix <- dist(train_batch, 'umatrix')

  biopsy_umatrix_hcl <- hcluster(biopsy_umatrix, k = 2)
  biopsy_umatrix_kmeans <- kcluster(biopsy_umatrix, k = 2, clust_fun = 'kmeans')
  biopsy_umatrix_pam <- kcluster(biopsy_umatrix, k = 2, clust_fun = 'pam')
  biopsy_umatrix_dbscan <- dbscan_cluster(biopsy_umatrix, eps = 0.82)

  biopsy_umatrix_hcl <-
    rename(biopsy_umatrix_hcl,
           nm = c('2' = 'cluster A',
                  '1' = 'cluster B'))

  model.frame(biopsy_umatrix_hcl)
  model.frame(biopsy_umatrix_kmeans)

  extract(biopsy_umatrix_pam, 'assignment')

  dist(biopsy_umatrix_hcl)
  dist(biopsy_umatrix_kmeans)

  ngroups(biopsy_umatrix_hcl)
  ngroups(biopsy_umatrix_kmeans)
  ngroups(biopsy_umatrix_dbscan)

  nobs(biopsy_umatrix_hcl)

  var(biopsy_umatrix_hcl)
  var(biopsy_umatrix_kmeans)
  var(biopsy_umatrix_pam)

  silhouette(biopsy_umatrix_hcl) %>% summary
  silhouette(biopsy_umatrix_kmeans) %>% plot
  silhouette(biopsy_umatrix_pam) %>% plot

  plot(biopsy_umatrix_hcl)
  plot(biopsy_umatrix_hcl, 'heat_map')

  plot(biopsy_umatrix_kmeans)
  plot(biopsy_umatrix_kmeans, 'heat_map')

  plot(biopsy_umatrix_pam)
  plot(biopsy_umatrix_pam, 'heat_map')

  components(biopsy_umatrix_hcl, kdim = 2, red_fun = 'umap') %>% plot
  components(biopsy_umatrix_kmeans, kdim = 2, red_fun = 'mds') %>% plot
  components(biopsy_umatrix_pam, kdim = 2, red_fun = 'pca') %>% plot

  cross_distance(biopsy_umatrix_hcl) %>% plot('mean')
  cross_distance(biopsy_umatrix_kmeans) %>% plot('mean')
  cross_distance(biopsy_umatrix_pam) %>% plot('histogram')

  impact(biopsy_umatrix_hcl)
  cv(biopsy_umatrix_kmeans)
  plot_clust_hm(biopsy_umatrix_pam)

  plot(biopsy_umatrix_hcl,
       type = 'components',
       kdim = 2,
       red_fun = 'umap')

# Extra SOM options for combi clustering -------

  train_combi_biopsy_extra <-
    combi_cluster(data = reduce(biopsy_lst$train, cbind),
                  distance_som = 'euclidean',
                  xdim = 8,
                  ydim = 8,
                  topo = 'hexagonal',
                  neighbourhood.fct = 'gaussian',
                  toroidal = FALSE,
                  rlen = 1000,
                  som_args = list(mode = 'online',
                                  alpha = c(0.01, 0.1)),
                  node_clust_fun = kcluster,
                  clust_fun = 'kmeans',
                  distance_nodes = 'euclidean',
                  k = 2,
                  seed = 12345)

  plot(train_combi_biopsy_extra, 'training')

  plot(train_combi_biopsy_extra)
  plot(train_combi_biopsy_extra, 'components', kdim = 2, with = 'umatrix', red_fun = 'umap')

  var(train_combi_biopsy_extra)
  silhouette(train_combi_biopsy_extra) %>%
    summary

  train_combi_biopsy_extra_cv <-
    cv(train_combi_biopsy_extra,
       nfolds = 10,
       type = 'som',
       kNN = 15,
       simple_vote = FALSE,
       .parallel = TRUE)

  train_combi_biopsy_varimp <-
    impact(train_combi_biopsy_extra,
           n_iter = 10,
           .parallel = TRUE)

  train_combi_biopsy_varimp %>% plot

  test_combi_biopsy_extra <-
    predict(train_combi_biopsy_extra,
            newdata = reduce(biopsy_lst$test, cbind),
            type = 'som')

  var(test_combi_biopsy_extra)

  plot(test_combi_biopsy_extra, 'components', with = 'data', kdim = 2, red_fun = 'umap')

  silhouette(test_combi_biopsy_extra) %>% summary
  silhouette(train_combi_biopsy_extra) %>% summary

  cross_distance(train_combi_biopsy_extra, .parallel = TRUE) %>% plot('mean')
  cross_distance(test_combi_biopsy_extra, .parallel = TRUE) %>% plot('mean')

  cross_distance(train_combi_biopsy_extra,
                 test_combi_biopsy_extra,
                 .parallel = TRUE) %>%
    plot('mean')

  plot_clust_hm(train_combi_biopsy_extra)
  plot_clust_hm(test_combi_biopsy_extra)

# Combined multi-layer SOM - unsupervised clustering -------

  train_batch_combi <-
    multi_cluster(data = biopsy_lst$train,
                  distance_method = c('manhattan', 'euclidean'),
                  xdim = 8,
                  ydim = 8,
                  topo = 'hexagonal',
                  neighbourhood.fct = 'gaussian',
                  toroidal = FALSE,
                  rlen = 1000,
                  som_args = list(mode = 'batch',
                                  alpha = c(0.01, 0.1)),
                  node_clust_fun = kcluster,
                  clust_fun = 'kmeans',
                  k = 2)

  test_batch_combi <- predict(train_batch_combi,
                              newdata = biopsy_lst$test)

  dist(train_batch_combi, 'umatrix')
  dist(test_batch_combi)

  nobs(train_batch_combi)
  nobs(test_batch_combi)

  ngroups(train_batch_combi)
  ngroups(test_batch_combi)

  var(train_batch_combi)
  var(test_batch_combi)

  silhouette(train_batch_combi) %>% plot
  silhouette(test_batch_combi) %>% plot

  silhouette(train_batch_combi) %>% summary
  silhouette(test_batch_combi) %>% summary

  components(train_batch_combi,
             kdim = 2,
             with = 'data',
             red_fun = 'umap')

  components(test_batch_combi,
             kdim = 2,
             with = 'data',
             red_fun = 'umap')

  components(train_batch_combi, kdim = 2, red_fun = 'mds', with = 'distance')
  components(test_batch_combi, kdim = 2, red_fun = 'mds', with = 'distance')

  plot(train_batch_combi)
  plot(train_batch_combi, 'training')

  plot(train_batch_combi,
       type = 'components',
       with = 'distance',
       kdim = 2,
       jitter_width = 0.1,
       jitter_height = 0.1)

  plot(test_batch_combi,
       type = 'components',
       with = 'distance',
       kdim = 2,
       jitter_width = 0.1,
       jitter_height = 0.1)

  plot(train_batch_combi,
       type = 'components',
       with = 'data',
       kdim = 2,
       jitter_width = 0.1,
       jitter_height = 0.1)

  plot(test_batch_combi,
       type = 'components',
       with = 'data',
       kdim = 2,
       jitter_width = 0.1,
       jitter_height = 0.1)

  plot_clust_hm(train_batch_combi)
  plot_clust_hm(test_batch_combi)

  impact(train_batch_combi, n_iter = 5, .parallel = TRUE) %>%
    plot

  cv(train_batch_combi, nfolds = 10, type = 'som', .parallel = TRUE)

  ## cross-distances

  cross_distance(test_batch_combi) %>% plot('mean')
  cross_distance(train_batch_combi) %>% plot('mean')

  cross_distance(train_batch_combi, test_batch_combi) %>% plot('mean')
  cross_distance(test_batch_combi, train_batch_combi) %>% plot('mean')


# END -----