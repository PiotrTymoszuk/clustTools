# Fitting a clustering structure with the hard-threshold KMEANS algorithm

# tools --------

  library(kohonen)

  library(clustTools)
  library(tidyverse)
  library(clusterHD)
  library(caret)
  library(somKernels)

  extract <- clustTools::extract
  map <- purrr::map

  library(patchwork)

# analysis data set -------

  ## the wines dataset pre-processing: elimination of invariant features
  ## (low variance to mean ratio) and mean centering

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
    center_data('mean')

  ## appending the parental data frame with wine classification
  ## it will be used for the final validation of the results

  my_wines <- my_wines %>%
    mutate(vintage = vintages)

  ## training: 100 randomly chosen observations
  ## the rest used for validation
  ## created with caret's createDataPartition() to keep
  ## the vintage distribution

  set.seed(12345)

  train_ids <- createDataPartition(my_wines$vintage, p = 100/177)[[1]]

  wine_lst <-
    list(train = my_wines[train_ids, ],
         test = my_wines[-train_ids, ]) %>%
    map(select, -vintage)

# variable distribution ------

  ## variable violin plots

  clust_variables <- names(wine_lst$train)

  var_distr_plot <- wine_lst$train %>%
    pivot_longer(cols = all_of(clust_variables),
                 names_to = 'variable',
                 values_to = 'z_score') %>%
    ggplot(aes(x = variable,
               y = z_score)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),
                scale = 'width',
                fill = 'cornsilk') +
    geom_point(shape = 16,
               size = 1,
               alpha = 0.5,
               position = position_jitter(width = 0.1))

# constructing clustering structures in the training data set -------

  train_clusters <- list()

  ## canonical KMEANS

  train_clusters$kmeans <- kcluster(data = wine_lst$train,
                                    distance_method = 'squared_euclidean',
                                    clust_fun = 'kmeans',
                                    k = 3)

  ## regularized KMEANS, finding lambda by CV-tuning

  train_tune <- tune_htk(data = wine_lst$train,
                         k = 3,
                         lambdas = seq(0, 1, by = 0.025),
                         select_stat = 'silhouette',
                         type = 'cv',
                         nfolds = 10,
                         kNN = 11,
                         .parallel = TRUE)

  summary(train_tune)

  plot(train_tune)

  train_clusters$htk <- extract(train_tune, 'analysis')

  ## renaming the clusters after their predominant vintages (checked elsewhere)

  train_clusters$kmeans <- train_clusters$kmeans %>%
    rename(c('2' = 'Barbera',
             '1' = 'Barolo',
             '3' = 'Grignolino'))

  train_clusters$htk <- train_clusters$htk %>%
    rename(c('3' = 'Barbera',
             '1' = 'Barolo',
             '2' = 'Grignolino'))

# performance in the training data and cross-validation -------

  train_clusters %>%
    map(summary)

  train_clusters %>%
    map(cv,
        type = 'propagation',
        nfolds = 10,
        kNN = 11,
        active_variables = TRUE) %>%
    map(summary) %>%
    map(select, ends_with('mean'))

# performance in the test data -------

  test_clusters <- train_clusters %>%
    map(predict.clust_analysis,
        newdata = wine_lst$test,
        type = 'propagation',
        kNN = 11,
        active_variables = TRUE)

  test_clusters %>%
    map(summary)

  ## cross distances

  map2(train_clusters,
       test_clusters,
       cross_distance) %>%
    map(plot, 'mean') %>%
    map(~.x +
          scale_fill_gradient2(low = 'firebrick',
                               mid = 'white',
                               high = 'steelblue',
                               midpoint = 20,
                               limits = c(7, 32)))

# summary of the performance stats -----

  clust_performance <-
    c(train_clusters, test_clusters) %>%
    map_dfr(summary) %>%
    mutate(dataset = c(rep('training', 2),
                       rep('test', 2)),
           algorithm = rep(c('KMEANS', 'HTKmeans'), 2))

  clust_performance %>%
    pivot_longer(cols = c(sil_width,
                          frac_misclassified,
                          frac_var,
                          frac_np),
                 names_to = 'statistic',
                 values_to = 'value') %>%
    ggplot(aes(x = value,
               y = algorithm,
               fill = dataset)) +
    geom_bar(stat = 'identity',
             color = 'black',
             position = position_dodge(0.9)) +
    facet_wrap(facets = vars(statistic),
               scales = 'free',
               labeller = as_labeller(c(sil_width = 'silhouette\nwidth',
                                        frac_misclassified = 'fraction\nmisclassified',
                                        frac_var = 'explained\nvariance',
                                        frac_np = 'fraction\npreserved neighbors'))) +
    labs(title = 'Clustering of the wines dataset')


# UMAP plots ------

  ## common UMAP layout

  umap_train <- train_clusters[[1]] %>%
    components(with = 'data',
               distance_method = 'cosine',
               kdim = 2,
               red_fun = 'umap',
               random_state = 12345)

  ## plotting the cluster assignment on the common UMAP layout

  kmeans_plots <- c(train_clusters,
                    test_clusters) %>%
    set_names(c('train_kmeans', 'train_htk',
                'test_kmeans', 'test_htk')) %>%
    map(plot,
        type = 'components',
        with = 'data',
        red_fun = 'umap',
        train_object = umap_train)

  ## and adjustments of plot elements using the ggplot interface

  kmeans_plots <-
    map2(kmeans_plots,
         c('Wines, training, KMEANS',
           'Wines, training, HTKmeans',
           'Wines, test, KMEANS',
           'Wines, test, HTKmeans'),
         ~.x +
           labs(title = .y) +
           theme(plot.tag = element_blank(),
                 legend.position = 'none'))

  kmeans_plots$train_kmeans +
    kmeans_plots$test_kmeans +
    kmeans_plots$train_htk +
    kmeans_plots$test_htk

# Heat maps of clustering factors -----

  plot_clust_hm(train_clusters$kmeans)
  plot_clust_hm(train_clusters$htk)

# variable importance -------

  var_imp <- train_clusters %>%
    map(impact,
        n_iter = 50,
        .parallel = TRUE)

  ## plotting and customizing the title

  var_imp_plots <- var_imp %>%
    map(plot) %>%
    map2(.,
         paste('Variable importance, wines',
               c('KMEANS', 'HTKmeans')),
         ~.x +
           geom_vline(xintercept = 0,
                      linetype = 'dashed') +
           labs(title = .y))

  var_imp_plots$kmeans +
    var_imp_plots$htk

# Accuracy of vintage prediction ------

  vintage_assignment <- my_wines %>%
    rownames_to_column('observation') %>%
    select(observation, vintage)

  ## assignment of samples to the clusters and vintages

  kmeans_assignment <- list(train = train_clusters$kmeans,
                            test = test_clusters$kmeans) %>%
    map(extract, 'assignment') %>%
    map(left_join, vintage_assignment, by = 'observation')

  htk_assignment <- list(train = train_clusters$htk,
                         test = test_clusters$htk) %>%
    map(extract, 'assignment') %>%
    map(left_join, vintage_assignment, by = 'observation')

  ## frequencies of vintages in the clusters

  kmeans_wine_counts <- kmeans_assignment %>%
    map(count, clust_id, vintage)

  htk_wine_counts <- htk_assignment %>%
    map(count, clust_id, vintage)

  ## receiver-operating characteristic

  kmeans_roc <- kmeans_assignment %>%
    map(transmute,
        obs = vintage,
        pred = clust_id) %>%
    map(as.data.frame) %>%
    map(multiClassSummary,
        lev = c('Barbera', 'Barolo', 'Grignolino'))

  htk_roc <- htk_assignment %>%
    map(transmute,
        obs = vintage,
        pred = clust_id) %>%
    map(as.data.frame) %>%
    map(multiClassSummary,
        lev = c('Barbera', 'Barolo', 'Grignolino'))

# summary of roc stats -------

  roc_summary <- c(kmeans_roc, htk_roc) %>%
    reduce(rbind) %>%
    as_tibble %>%
    mutate(dataset = rep(c('training', 'test'), 2),
           algorithm = c(rep('KMEANS', 2),
                         rep('HTKmeans', 2)))

  roc_summary %>%
    pivot_longer(cols = c(Accuracy,
                          Kappa,
                          Mean_Sensitivity,
                          Mean_Specificity),
                 names_to = 'statistic',
                 values_to = 'value') %>%
    ggplot(aes(x = value,
               y = algorithm,
               fill = dataset)) +
    geom_bar(stat = 'identity',
             color = 'black',
             position = position_dodge(0.9)) +
    facet_wrap(facets = vars(statistic),
               scales = 'free') +
    labs(title = 'Vintage prediction')


# combi clustering --------

  ## training, comparing kmeans with regularized kmeans

  train_combi <- list()

  train_combi$kmeans <- combi_cluster(data = wine_lst$train,
                                      distance_som = 'cosine',
                                      xdim = 5,
                                      ydim = 4,
                                      topo = 'hexagonal',
                                      neighbourhood.fct = 'gaussian',
                                      toroidal = FALSE,
                                      rlen = 3000,
                                      node_clust_fun = kcluster,
                                      distance_nodes = 'squared_euclidean',
                                      k = 3,
                                      clust_fun = 'kmeans')

  train_combi$htk <- combi_cluster(data = wine_lst$train,
                                   distance_som = 'cosine',
                                   xdim = 5,
                                   ydim = 4,
                                   topo = 'hexagonal',
                                   neighbourhood.fct = 'gaussian',
                                   toroidal = FALSE,
                                   rlen = 3000,
                                   node_clust_fun = htk_cluster,
                                   lambdas = seq(0, 1, by = 0.025),
                                   standardize = FALSE,
                                   k = 3)

  ## comparing performance in the training data set and cross-validation

  train_combi %>%
    map(summary)

  train_combi %>%
    map(cv,
        nfolds = 10,
        type = 'propagation',
        active_variables = FALSE,
        .parallel = TRUE) %>%
    map(summary) %>%
    map(select, ends_with('mean'))

  ## comparing performance in the test data set

  train_combi %>%
    map(predict,
        newdata = wine_lst$test,
        type = 'propagation',
        active_variables = FALSE) %>%
    map(summary)

  ## variable importance

  train_combi %>%
    map(impact,
        n_iter = 10,
        .parallel = TRUE)


# END ------
