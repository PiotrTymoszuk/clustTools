# Working with the MASS's Cars93 data set

# packages and data -------

  ## Cars93 retrieved from MASS
  ## the package order matters!

  library(MASS)
  library(kohonen)

  library(tidyverse)
  library(clustTools)
  library(somKernels)

  ## patchwork for plot panels

  library(patchwork)

  ## the Cars93 data set. No hold-out, cross-validation user instead,
  ## since it is a pure explanatory analysis
  ##
  ## we'll check if multi-layer SOM clustering can recapitulate
  ## car classification provided by the `Type` variable

  my_cars <- MASS::Cars93 %>%
    select(Make, Type, Price,
           MPG.city, MPG.highway,
           AirBags, DriveTrain,
           Cylinders, EngineSize, Horsepower,
           Man.trans.avail, Man.trans.avail,
           Fuel.tank.capacity, Passengers,
           Length, Weight, Origin) %>%
    as_tibble

  ## formatting factor variables

  my_cars <- my_cars %>%
    mutate(AirBags = car::recode(AirBags,
                                 "'None' = '0';
                                 'Driver only' = '1';
                                 'Driver & Passenger' = '2'",
                                 as.factor = FALSE),
           AirBags = as.numeric(AirBags),
           Man.trans.avail = car::recode(Man.trans.avail,
                                         "'No' = '0';
                                         'Yes' = '1'",
                                         as.factor = FALSE),
           Man.trans.avail = as.numeric(Man.trans.avail),
           Origin = car::recode(Origin,
                                "'USA' = '1';
                                'non-USA' = '0'",
                                as.factor = FALSE),
           Origin = as.numeric(Origin),
           Cylinders = as.numeric(Cylinders))

  my_cars <- my_cars %>%
    filter(complete.cases(.))

  ## data layers: numeric, ordinal and binary
  ## based on the variable class. Numeric and ordinal variables
  ## are normalized and median-centered

  car_lst <-
    list(numeric = c('Price', 'MPG.city', 'MPG.highway', 'EngineSize',
                     'Horsepower', 'Fuel.tank.capacity', 'Length',
                     'Weight'),
         ordinal = c('AirBags', 'Cylinders', 'Passengers'),
         binary = c('Man.trans.avail', 'Origin')) %>%
    map(~my_cars[c(.x, 'Make')]) %>%
    map(column_to_rownames, 'Make')

  car_lst[c("numeric", "ordinal")] <-
    car_lst[c("numeric", "ordinal")] %>%
    map(center_data, type = 'median')

# fitting the SOM -------

  ## the layers will be handled by three various distances:
  ## Euclidean for numeric variables, Manhattan for ordinal variables
  ## and Tanimoto for binary features

  car_som <-
    multi_cluster(data = car_lst,
                  distance_method = c('cosine', 'manhattan', 'tanimoto'),
                  xdim = 5,
                  ydim = 6,
                  topo = 'hexagonal',
                  neighbourhood.fct = 'gaussian',
                  toroidal = FALSE,
                  rlen = 3000,
                  som_args = list(mode = 'online',
                                  alpha = c(0.01, 0.2),
                                  user.weights = c(1, 1, 0.5),
                                  normalizeDataLayers = FALSE),
                  node_clust_fun = kcluster,
                  k = 3,
                  clust_fun = 'pam')

# convergence and cluster number --------

  ## SOM convergence

  car_convergence <- plot(car_som, type = 'training')

  car_convergence <- car_convergence$observation +
    facet_grid(`SOM layer` ~ .,
               scales = 'free') +
    theme(plot.tag = element_blank(),
          legend.position = 'none')

  ## WSS and silhouette plots

  car_clust_number <- plot(car_som, type = 'diagnostic')$node %>%
    map(~.x + theme(plot.tag = element_blank()))

  car_clust_number$wss +
    car_clust_number$silhouette

# Cross-validated numeric stats ------

  car_cross_validation <-
    cv(car_som,
       nfolds = 10,
       type = 'som') %>%
    summary

  car_cross_validation %>%
    select(ends_with('mean'))

# Silhouette and heat map plots -------

  ## silhouette values for particular clusters

  car_sil_plot <- car_som %>%
    silhouette %>%
    plot

  ## heat map of pairwise distances

  car_heat_map <-
    plot(car_som, type = 'heat_map')$node +
    labs(title = 'Distance between SOM nodes') +
    theme(plot.tag = element_blank())

# Heat map visualization of the clustering features ------

  ## plot_clust_hm() returns a list
  ## of heat maps for multi-layer data
  ##
  ## you may adjust them separately

  car_feature_hm <- car_som %>%
    plot_clust_hm

  car_feature_hm <- car_feature_hm %>%
    map(~.x +
          theme(legend.position = 'bottom',
                axis.text.x = element_blank(),
                plot.tag = element_blank()))

  car_feature_hm$numeric <-
    car_feature_hm$numeric +
    scale_fill_gradient2(low = 'steelblue',
                         mid = 'black',
                         high = 'firebrick',
                         midpoint = 0,
                         limits = c(-3, 3),
                         oob = scales::squish) +
    labs(title = 'Cars93: numeric features')

  car_feature_hm$ordinal <-
    car_feature_hm$ordinal +
    scale_fill_gradient2(low = 'steelblue',
                         mid = 'black',
                         high = 'firebrick',
                         midpoint = 0,
                         limits = c(-2.5, 2.5),
                         oob = scales::squish) +
    guides(fill = 'legend') +
    labs(title = 'Cars93: ordinal features')

  car_feature_hm$binary <-
    car_feature_hm$binary +
    scale_fill_gradient(low = 'steelblue',
                        high = 'firebrick') +
    guides(fill = 'legend') +
    labs(title = 'Cars93: binary features')

  car_feature_hm$numeric +
    car_feature_hm$ordinal +
    car_feature_hm$binary +
    plot_layout(ncol = 2)

# Automobile category in the clusters ------

  car_assignment <- car_som %>%
    extract('assignment') %>%
    mutate(Make = observation) %>%
    left_join(my_cars[c('Make', 'Type')],
              by = 'Make')

  car_counts <-
    count(car_assignment, clust_id, Type) %>%
    group_by(clust_id) %>%
    mutate(percent = n/sum(n) * 100) %>%
    ungroup

  car_count_plot <- car_counts %>%
    ggplot(aes(x = clust_id,
               y = percent,
               fill = Type)) +
    geom_bar(position = 'stack',
             stat = 'identity',
             color = 'black') +
    scale_fill_viridis_d() +
    labs(title = 'Car type distribution in the clusters',
         x = 'Cluster of Cars93')

# Variable importance -------

  car_importance <- impact(car_som,
                           n_iter = 50,
                           .parallel = TRUE)

  plot(car_importance) +
    geom_vline(xintercept = 0,
               linetype = 'dashed')

# END -----
