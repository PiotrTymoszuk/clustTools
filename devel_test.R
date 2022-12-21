# Tests the functionality during development

  library(clustTools)
  library(somKernels)

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

  plot_dendro(test_hcl,
              k = 3,
              labels = FALSE,
              cluster_colors = c('firebrick',
                                 'steelblue',
                                 'darkolivegreen',
                                 'black'))

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

  test_pred <- predict(object = test_som,
                       newdata = new_data,
                       type = 'propagation',
                       simple_vote = FALSE)

  plot(test_pred, type = 'components', with = 'data', red_fun = 'pca')

  plot(test_som, type = 'components', with = 'data', red_fun = 'pca')

  test_pred_red <- propagate(object = test_pam,
                             newdata = reduce_data(data = new_data,
                                                   distance_method = 'cosine',
                                                   kdim = 3,
                                                   red_fun = 'umap'))

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

  test_cv <- cv_cluster(data = test_red$umap,
                        nfolds = 10,
                        kNN = 5,
                        clustering_fun = kcluster,
                        seed = 1234,
                        distance_method = 'cosine',
                        kernel_fun = identity,
                        clust_fun = 'pam',
                        k = 2,
                        .parallel = FALSE)

  cv(test_hcl,
     nfolds = 10,
     kNN = 5)

  cv(test_pam,
     nfolds = 10,
     kNN = 5)

  cv(test_som,
     nfolds = 10,
     kNN = 5)

  cv(test_dbscan,
     nfolds = 10,
     kNN = 5)

# Variable importance and its OOP -----

  test_impact <- importance_cluster(data = test_red$umap,
                                    clustering_fun = dbscan_cluster,
                                    seed = 1234,
                                    distance_method = 'euclidean',
                                    eps = 1,
                                    .parallel = FALSE)

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
                              k = 2,
                              hc_method = 'complete')

  extract(test_combi, 'clust_object')

  is_combi_analysis(test_combi)

  dist(test_combi)

  components(test_combi, kdim = 3, red_fun = 'umap', with = 'data')

  var(test_combi)

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

  plot(test_new_combi, type = 'components', red_fun = 'pca', with = 'data', k = 2)

  plot(test_combi, type = 'components', red_fun = 'pca', with = 'data', k = 2)

# Cross-validation -------

  cv_cluster(data = test_data,
             nfolds = 10,
             kNN = 5,
             clustering_fun = combi_cluster,
             distance_som = 'manhattan',
             xdim = 5,
             ydim = 4,
             topo = 'hexagonal',
             neighbourhood.fct = 'gaussian',
             toroidal = FALSE,
             rlen = 1000,
             node_clust_fun = hcluster,
             kernel_fun = identity,
             k = 2,
             hc_method = 'complete')

  cv(test_combi,
     nfolds = 10,
     kNN = 5,
     .parallel = FALSE)

# Feature importance ------

  importance_cluster(data = test_data,
                     clustering_fun = combi_cluster,
                     distance_som = 'manhattan',
                     xdim = 5,
                     ydim = 4,
                     topo = 'hexagonal',
                     neighbourhood.fct = 'gaussian',
                     toroidal = FALSE,
                     rlen = 1000,
                     node_clust_fun = hcluster,
                     k = 2,
                     hc_method = 'complete')

  plot(impact(test_combi))

# Feature plotting ------

  ft_clust <- hcluster(data = t(test_data), distance_method = 'euclidean', k = 2)

  plot(ft_clust)

  plot_clust_hm(x_object = test_combi,
                y_object = ft_clust)














