# GitHub examples for dimensionality reduction

# tools -------

  library(tidyverse)
  library(clustTools)

  library(patchwork)

# Example data: the biopsy data set --------

  ## the MASS biopsy data set.
  ## Including only complete cases in the analysis
  ## The ID is not unique, since multiple samples
  ## are available for the same patient, so I'm redefining it
  ## no scaling of the clustering features V1 - V9, they are
  ## on the same scale

  my_biopsy <- MASS::biopsy %>%
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

# Devel -------

  my_biopsy %>%
    select(starts_with('V')) %>%
    reduce_data(kdim = 3, red_fun = 'pca') %>%
    np %>%
    summary

  biopsy_lst$train %>%
    select(starts_with('V')) %>%
    reduce_data(kdim = 3, red_fun = 'umap') %>%
    np %>%
    summary
