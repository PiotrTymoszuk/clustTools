[![R](https://github.com/PiotrTymoszuk/clustTools/actions/workflows/r.yml/badge.svg)](https://github.com/PiotrTymoszuk/clustTools/actions/workflows/r.yml)

<img src="https://github.com/PiotrTymoszuk/clustTools/assets/80723424/fcf39384-39f7-43f8-a396-b546a3218398" width="20%" height="20%" align = "right">

# clustTools
Comprehensive dimensionality reduction and cluster analysis toolset

## Description

The `clustTools` package provides a medley of functions used for seemless integration of various dimensionality reduction methods (multi-dimensional scaling/MDS, principal component analysis/PCA, uniform manifold approximation and projection/UMAP or factor analysis/FA), clustering (hierarchical clustering, K-means, PAM and density DBSCAN clustering) and self-orgnizing map (SOM) analyses. In addition, a set of functions is provided for visualization, quality control and cross-validation of the clustering results and semi-supervised clusetring.

## Installation

You may easily fetch the package and its dependency `somKernels` with `devtools`: 

```r

devtools::install_github('PiotrTymoszuk/clustTools')
devtools::install_github('PiotrTymoszuk/somKernels')

```
## Terms of use

The package is available under a [GPL-3 license](https://github.com/PiotrTymoszuk/clustTools/blob/main/LICENSE).

## Contact

The package maintainer is [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com).

## Acknowledgements

`clustTools` uses tools provided by the [rlang](https://rlang.r-lib.org/), [tidyverse](https://www.tidyverse.org/), [stringi](https://stringi.gagolewski.com/), [caret](https://topepo.github.io/caret/), [coxed](https://cran.r-project.org/web/packages/coxed/index.html), [dbscan](https://cran.r-project.org/web/packages/dbscan/index.html), [dendextend](https://github.com/talgalili/dendextend), [factoextra](https://cran.r-project.org/web/packages/factoextra/index.html), [furrr](https://furrr.futureverse.org/), [future](https://future.futureverse.org/), [ggrepel](https://ggrepel.slowkow.com/), [kohonen](https://cran.r-project.org/web/packages/kohonen/index.html), [nomclust](https://cran.r-project.org/web/packages/nomclust/index.html), [pcaPP](https://github.com/cran/pcaPP), [philentropy](https://github.com/drostlab/philentropy), [umap](https://github.com/tkonopka/umap), [utils](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/00Index.html), [cluster](https://cran.r-project.org/web/packages/cluster/index.html), and [generics](https://github.com/r-lib/generics). Many thanks to their developers, maintainers and contributors.


## Usage

### Hierarchical and mean/medoid clustering

<details>

`clustTools` represents a one-stop shop for construction, diagnostic and validation of clustering solutions, which are scattered between many excellent R packages. Let's take a look at its basic functionalities by semi-supervised clustering of the `biopsy` data set insluded in the R's `MASS` package. This portion of data stores results of a clinical study on breast lesion biopsies, which are classified as benign or malignant based on 9 morphological and cytological variables assessed by a pathologist on a 1 - 10 scale each. The `biopsy` data set will be split into a training and test portion. I'll use the training portion for choice of the best performing clustering algorithm and cluster number. The test portion will be used for the final validation of the clustering results. The algorithms of interest are Ward's hierarchical clustering (internally implemented by `stats`), kmeans (`stats` package) and PAM (partitioning around medoids, package `cluster`). 

```r

  ## data source, the loading order matters
  ## load MASS as the first package

  library(MASS)

  library(tidyverse)
  library(clustTools)
  library(somKernels)

  ## patchwork package for plot panels

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

```

```r

> my_biopsy
# A tibble: 683 × 12
   ID                V1    V2    V3    V4    V5    V6    V7    V8    V9 class     observation
   <chr>          <int> <int> <int> <int> <int> <int> <int> <int> <int> <fct>     <chr>      
 1 obs_1_1000025      5     1     1     1     2     1     3     1     1 benign    obs_1      
 2 obs_2_1002945      5     4     4     5     7    10     3     2     1 benign    obs_2      
 3 obs_3_1015425      3     1     1     1     2     2     3     1     1 benign    obs_3      
 4 obs_4_1016277      6     8     8     1     3     4     3     7     1 benign    obs_4      
 5 obs_5_1017023      4     1     1     3     2     1     3     1     1 benign    obs_5      
 6 obs_6_1017122      8    10    10     8     7    10     9     7     1 malignant obs_6      
 7 obs_7_1018099      1     1     1     1     2    10     3     1     1 benign    obs_7      
 8 obs_8_1018561      2     1     2     1     2     1     3     1     1 benign    obs_8      
 9 obs_9_1033078      2     1     1     1     2     1     1     1     5 benign    obs_9      
10 obs_10_1033078     4     2     1     1     2     1     2     1     1 benign    obs_10     
# … with 673 more rows
# ℹ Use `print(n = ...)` to see more rows

```

```r

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

```
I'm checking first the clustering tendency of the training and test data subsets with the [Hopkins statistic](https://en.wikipedia.org/wiki/Hopkins_statistic) with the `get_clust_tendency()` functions implemented in the great package `factoextra`. Of note, Hopkins statistic ranges from 0 (completely homogenous distribution) to 1 (highly clustered distribution). Values around 0.5 indicate a random normal-like ditribution. With the Hopkins statistic values of approximately 0.7, the `biopsy` data set demonstrates a weak-to-moderate clustering tendency:

```r

biopsy_clust_tendency <- biopsy_lst %>% 
    map(get_clust_tendency, 
        n = 200)
```
```r
> biopsy_clust_tendency$train$hopkins_stat
[1] 0.7148951
> biopsy_clust_tendency$test$hopkins_stat
[1] 0.7307125
```

I'm using the `clustTools` functions `hcluster()` and `kcluster()` to cluster the training portion of the data set with the hierarchical clusetring, kmeans and PAM. I'm also trying two distance measures: Euclidean and cosine. The initially guessed cluster number is set to `k = 2` as expected for the data set containing two types of samples. Importantly, those clustering functions take numeric data frames or matrices as the first argument and return objects of `clust_analysis` class, which can be, irrespectively of the clustering algorithm, subjected to a consistent quality control and validation.

```r

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

```

```r

> clust_objects$pam_cosine
Clustering analysis with PAM and cosine distance method.
Cluster assignment:
# A tibble: 383 × 2
   observation     clust_id
   <chr>           <fct>   
 1 obs_284_601265  1       
 2 obs_101_1167471 1       
 3 obs_623_1277792 1       
 4 obs_645_1339781 2       
 5 obs_400_1238948 2       
 6 obs_98_1166630  2       
 7 obs_103_1168736 2       
 8 obs_602_1042252 1       
 9 obs_326_806423  2       
10 obs_79_1137156  1       
# … with 373 more rows
# ℹ Use `print(n = ...)` to see more rows

```
The `clust_analysis` object bundles the clustering data (called by `model.frame()`), distance matrix (`dist()`), a data frame with the cluster assignment (`extract(type = 'assignment')`) and other internal elements allowing for re-fitting of the clustering structure:

```r

> head(model.frame(clust_objects$pam_cosine))
                V1 V2 V3 V4 V5 V6 V7 V8 V9
obs_284_601265  10  4  4  6  2 10  2  3  1
obs_101_1167471  4  1  2  1  2  1  3  1  1
obs_623_1277792  4  1  1  1  2  1  1  1  1
obs_645_1339781  1  1  1  1  2  1  2  1  1
obs_400_1238948  8  5  6  2  3 10  6  6  1
obs_98_1166630   7  5  6 10  5 10  7  9  4

> as.matrix(dist(clust_objects$pam_cosine))[1:5, 1:5]
                obs_284_601265 obs_101_1167471 obs_623_1277792 obs_645_1339781 obs_400_1238948
obs_284_601265      0.00000000      0.21342699      0.15789388       0.2976895      0.08462418
obs_101_1167471     0.21342699      0.00000000      0.06341419       0.1204068      0.15371658
obs_623_1277792     0.15789388      0.06341419      0.00000000       0.2049536      0.19244959
obs_645_1339781     0.29768946      0.12040679      0.20495361       0.0000000      0.18009751
obs_400_1238948     0.08462418      0.15371658      0.19244959       0.1800975      0.00000000

> extract(clust_objects$pam_cosine, type = 'assignment')
# A tibble: 383 × 2
   observation     clust_id
   <chr>           <fct>   
 1 obs_284_601265  1       
 2 obs_101_1167471 1       
 3 obs_623_1277792 1       
 4 obs_645_1339781 2       
 5 obs_400_1238948 2       
 6 obs_98_1166630  2       
 7 obs_103_1168736 2       
 8 obs_602_1042252 1       
 9 obs_326_806423  2       
10 obs_79_1137156  1       
# … with 373 more rows
# ℹ Use `print(n = ...)` to see more rows

```
</details>

### Cluster number

<details>

The package offers a visual help for verifying your cluster number guess. By calling `plot(type = 'diagnostic')` for a `clust_analysis` object, plots of dendrograms (hierarchical clustering only), within-cluster sum of squares and mean silhouette statistic values for clustering structures with varyiing cluster numbers can retrieved. All plots generated with the package's tools are `ggplot` graphical objects and can be easily customized by the user:

```r

  ## plot(type = 'diagnostic') returns
  ## a list of diagnostic plots for determination
  ## of the optimal cluster number
  
  cl_number_plots <- clust_objects %>% 
    map(plot, 
        type = 'diagnostic') %>% 
    map(map, 
        ~.x + 
          theme(plot.tag.position = 'bottom'))

```

```r
## for hierarchical clustering dendrograms are available

> cl_number_plots$hcl_euclidean$dendrogram

```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/b5086ac2-06fd-4256-948a-4dcba8367cc7)

```r
## within-cluster sum of squares and mean silhouette plots
## are generated for most algorithms

> cl_number_plots$pam_euclidean$wss + cl_number_plots$pam_euclidean$silhouette
```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/c4fee14a-00d9-4062-b4ba-31bb9c662feb)

As indicated by the two main branches of the dendrogram, the bend of the within-cluster sum of square curve, and the peak of the mean [silhouette statistic](https://en.wikipedia.org/wiki/Silhouette_(clustering)), the two-cluster solution seem to be a resonable guess.

</details>

## Quality control of clustering solutions

<details>

The package offers basically two numeric measures of cluster quality:

* [silhouette statistic](https://en.wikipedia.org/wiki/Silhouette_(clustering)), which gauges the quality of discrimination between the clusters

* _explained clustering variance_ defined as a ratio of the total between-cluster sum of squares to the total sum of squares. As such, explained clustering variance works in a similar way to R^2 or ANOVA

Function `silhouette()` computes silhouette values for each observation, mean values for the whole object and particular clusters can be retrieved with the `summary` method. Genarally, silhouettes range from -1 to 1, with high values indicative of good separation of the given observation or cluster from other clusers:

```r

> clust_objects$kmeans_cosine %>% 
+ silhouette %>% 
+ summary
# A tibble: 3 × 13
  clust_id     n n_negative perc_negative  mean    sd median   q025    q25   q75  q975    min   max
  <fct>    <int>      <int>         <dbl> <dbl> <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl>  <dbl> <dbl>
1 global     383         37          9.66 0.357 0.251  0.313 -0.115 0.213  0.574 0.736 -0.202 0.764
2 1          165          0          0    0.567 0.172  0.638  0.219 0.433  0.697 0.756  0.156 0.764
3 2          218         37         17.0  0.198 0.171  0.219 -0.135 0.0717 0.302 0.482 -0.202 0.522

```
For kmeans clustering of the `biopsy` training subset with cosine distance, low mean silhouette for the second large cluster indicate a possibly poor distinctness. This may be investigated in a plot in more details:

```r
> clust_objects$kmeans_cosine %>% 
+ silhouette %>% 
+ plot
```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/22525c18-e5ad-4a05-b449-b62031868e2e)

In particular, 18% of the cluster 2 observations had negative statistic values, which indicate that they are most likely misclassified. 

Let's compare mean silhouettes and clusetring variances for different clustering solutions. The later can be computed for `clust_analysis` objects by calling `var()` with the `frac_var` element of the result storing the fraction of explained clustering variance which concerns us most.

```r
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

```
In the `biopsy` data set, Euclidean distance yields a generally better separation between the clusters as the cosine metics as unequivocally demonstrated by mean silhouettes. The PAM/Euclidean solution has the lowest count of potentially misclassified observations with negative silhouette widths.

```r
> cl_silhouettes
# A tibble: 6 × 14
  clustering_algorithm clust_id     n n_negative perc_negative  mean    sd median    q025   q25   q75  q975     min   max
  <chr>                <fct>    <int>      <int>         <dbl> <dbl> <dbl>  <dbl>   <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>
1 hcl_euclidean        global     383         25          6.53 0.567 0.317  0.747 -0.253  0.339 0.816 0.849 -0.602  0.849
2 hcl_cosine           global     383         65         17.0  0.318 0.323  0.298 -0.356  0.136 0.592 0.784 -0.521  0.805
3 kmeans_euclidean     global     383         17          4.44 0.583 0.282  0.748 -0.106  0.365 0.810 0.843 -0.291  0.843
4 kmeans_cosine        global     383         37          9.66 0.357 0.251  0.313 -0.115  0.213 0.574 0.736 -0.202  0.764
5 pam_euclidean        global     383          6          1.57 0.588 0.252  0.738  0.0642 0.372 0.795 0.825 -0.120  0.825
6 pam_cosine           global     383         11          2.87 0.359 0.202  0.332 -0.0143 0.248 0.524 0.676 -0.0963 0.695

```

In turn, kmeans with Euclidean distance provides the most informative clustering procedure as demonstrated by the highest explained clustering variance:

```r
>   cl_variances
   hcl_euclidean       hcl_cosine kmeans_euclidean    kmeans_cosine    pam_euclidean       pam_cosine 
       0.8082272        0.2898587        0.8303221        0.3172726        0.8173912        0.2914502 
```

The `clustTools` package offers also a rich set of visual diagnostic tools. Heat maps generated with `plot(type = 'heat_map')` and `cross_distance() %>% plot(type = 'mean')` are useful at assessing distinctness of the clusters. Let's check them out for the PAM algorithm - eye-balling of the plots confirms the poor separation between the clusters by PAM/cosine:

```r
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
```

```r
> pam_heat_maps$pam_euclidean + pam_heat_maps$pam_cosine

```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/ea7c2f76-8f8a-4970-9803-2d89d826a4d4)

```r

 pam_mean_heat_maps <- clust_objects[c("pam_euclidean", "pam_cosine")] %>% 
    map(cross_distance) %>% 
    map(plot, type = 'mean') %>% 
    map(~.x + theme(legend.position = 'bottom'))
    
```

![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/9e55beb1-49f0-4eba-993e-84657e210bfb)
  
</details>

## Cross-validation

<details>

Unsupervised clustering objects resemble multi-parameter classification machine learning models in multiple aspects. By principle, they can be cross-validated just like any other machine learning models. In brief, in each cross-validation fold, the training portion is used for fitting of the clustering structure with the same algorithm as the parent clustering solution. The test portion of the fold is used for calculation of numeric stats and assessment of cluster assignment accuracy as compared with the parental clustering solution. Such approach has few advantages as compared with computation of numeric statistics only for the training data subset like a more robust handling of atypical observations and assessment of over-fitting. 
The default method of assignment of test portion observations to the clusters defined in the trainig portion of the cross-validation fold is so-called k-nearest neighbor label propagation. This prediction algorithm is resonably fast and agnostic to the clustering procedure. 

Here, we'll cross-validate the `biopsy` training subset clustering objects generated with the best performing Euclidean distance. This is accomplished with the `cv()` function, which returns out-of-fold predictions, numeric statistics for the folds along with means and confidence intervals of cluster assignment accuracy, classification error, explained clustering variance and mean silhouette. The later statistics can be retrieved from the `cv()` output by calling `summary()`:

```r

  euclidean_cv_stats <- 
    clust_objects[c("hcl_euclidean", "kmeans_euclidean", "pam_euclidean")] %>% 
    map(cv) %>% 
    map(summary)

```

```r

> euclidean_cv_stats %>% 
+ map(select, ends_with('mean'))
$hcl_euclidean
# A tibble: 1 × 4
  accuracy_mean error_mean frac_var_mean sil_width_mean
          <dbl>      <dbl>         <dbl>          <dbl>
1         0.763      0.237         0.805          0.586

$kmeans_euclidean
# A tibble: 1 × 4
  accuracy_mean error_mean frac_var_mean sil_width_mean
          <dbl>      <dbl>         <dbl>          <dbl>
1         0.395      0.605         0.813          0.588

$pam_euclidean
# A tibble: 1 × 4
  accuracy_mean error_mean frac_var_mean sil_width_mean
          <dbl>      <dbl>         <dbl>          <dbl>
1         0.790      0.210         0.801          0.586

```

In case of kmeans clustering and of other algorithms with stochastic determination of the initial cluster centers, the accuracy and classification error are expected to be poor. Still, we can use fractions of explained variations (`frac_var`) and mean silhouette widths (`sil_width`), to compare the clustering procedures. While all of them have similar discriminatory power, the kmeans/Euclidean colution displays the largest fraction of explained variance.
  
</details>

## Semi-supervised clustering

<details>

As for cross-validation, there's a possibility to apply the k-nearest neighbor (kNN) label propagation algorithm to predict the cluster assignment in a hold-out subset or external validation data set. Here we'll predict the cluster assignment defined by the kmeans/Euclidean clustering in the training portion for the test subset of the `biopsy` data set and check quality of the cluster assignment with explained clustering variance and silhouette width. 
The prediction is done with the `predict()` function, which takes a clustering analysis object as the first argument and a numeric data frame or matrix as `newdata`; the `type` argument is set to 'propagation' for the kNN procedure. The output is an instance of `clust_analysis` class, which allows for direct comparison of the training and test data clustering.

```r

  ## predictions for KMEANS/Euclidean
  
  kmeans_clusters <- list()
  
  kmeans_clusters$training <- clust_objects$kmeans_euclidean
  
  kmeans_clusters$test <- 
    predict(clust_objects$kmeans_euclidean, 
            newdata = biopsy_lst$test, 
            type = 'propagation')

```

Let's compare explained variances and mean silhouette witdths of the training clustering and predictions: 

```r

  ## explained variance and silhouette width
  
  kmeans_variance <- kmeans_clusters %>% 
    map(var) %>% 
    map_dbl(~.x$frac_var)
  
  kmeans_silhouette <- kmeans_clusters %>% 
    map(silhouette) %>% 
    map(summary)

```
```r
>   kmeans_variance
 training      test 
0.8303221 0.8385571

>   kmeans_silhouette
$training
# A tibble: 3 × 13
  clust_id     n n_negative perc_negative  mean    sd median   q025   q25   q75  q975    min   max
  <fct>    <int>      <int>         <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl>
1 global     383         17          4.44 0.583 0.282  0.748 -0.106 0.365 0.810 0.843 -0.291 0.843
2 1          249          0          0    0.764 0.102  0.805  0.443 0.750 0.819 0.843  0.282 0.843
3 2          134         17         12.7  0.248 0.184  0.297 -0.182 0.166 0.388 0.465 -0.291 0.473

$test
# A tibble: 3 × 13
  clust_id     n n_negative perc_negative  mean    sd median    q025   q25   q75  q975    min   max
  <fct>    <int>      <int>         <dbl> <dbl> <dbl>  <dbl>   <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl>
1 global     300          5          1.67 0.608 0.262  0.771  0.0694 0.391 0.815 0.849 -0.230 0.849
2 1          198          0          0    0.774 0.104  0.806  0.415  0.775 0.822 0.849  0.195 0.849
3 2          102          5          4.90 0.287 0.156  0.311 -0.141  0.192 0.402 0.480 -0.230 0.500

```

The cluster prediction in the test portion of the `biopsy` data set seems to be equally informative as in the parental clustering solution. Moreover, the discriminatory performance measured by silhouettes is even marginally better for the test data than for the training subset. The distinctess of the clusters can be assessed by heat maps of pairwise distances between observations:

```r

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

```

```r

> kmeans_heat_maps$training + kmeans_heat_maps$test

```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/e9725bda-1939-4a0a-85a9-7168b56e3200)

By calling `cross_distance()` for the training and test data clustering structures, we can also assess how the training and test clusters relate to each other:

```r
kmeans_cross_dists <- 
    cross_distance(kmeans_clusters$training, 
                   kmeans_clusters$test)

```

```r
> kmeans_cross_dists %>% 
+ plot(type ='mean')
```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/4d0a1091-bb61-4e0e-bd22-d7b28e80dfa6)

Finally, let's compare the distribution of malignant and benign samples between the clusters. This analysis reveals, that nearly all benign samples were assigned to the cluster 1 and malingnant samples were located in the cluster 2:

```r

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

```
```r
> kmeans_clust_counts
$training
# A tibble: 4 × 3
  clust_id class         n
  <fct>    <fct>     <int>
1 1        benign      240
2 1        malignant     9
3 2        benign        7
4 2        malignant   127

$test
# A tibble: 4 × 3
  clust_id class         n
  <fct>    <fct>     <int>
1 1        benign      193
2 1        malignant     5
3 2        benign        4
4 2        malignant    98

```
We can corroborate it further by computing accuracy, Cohen's kappa and receiver-operating characteristic with `caret`:

```r
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
```
The kmeans/Euclidean clustering of the `biopsy` data set allows an excellent discrimination of the malignant and benign breast biopsy samples:

```r

>  kmeans_roc_stats
$training
         Accuracy             Kappa                F1       Sensitivity       Specificity    Pos_Pred_Value    Neg_Pred_Value 
        0.9582245         0.9084854         0.9677419         0.9716599         0.9338235         0.9638554         0.9477612 
        Precision            Recall    Detection_Rate Balanced_Accuracy 
        0.9638554         0.9716599         0.6266319         0.9527417 

$test
         Accuracy             Kappa                F1       Sensitivity       Specificity    Pos_Pred_Value    Neg_Pred_Value 
        0.9700000         0.9333136         0.9772152         0.9796954         0.9514563         0.9747475         0.9607843 
        Precision            Recall    Detection_Rate Balanced_Accuracy 
        0.9747475         0.9796954         0.6433333         0.9655759 

```
  
</details>

## Advanced visualization options

<details>

Visualization of clustering structure with help of dimensionality reduction methods has a long tradition. By calling `plot(type = 'components')` for a clustering analysis objects, the `clustTools` package generates a scatter plot of observations's scores associated with the first two dimensions of prlincipal component analysis (PCA), multi-dimensional scaling (MDS) or Uniform Manifold Approximation and Projection (UMAP). By default, the dimensionality reduction is done for the distance matrix, by specifying `with = 'data'`, the user can request it for the genuine data set.

```r

## plots of MDS of the distance matrix, as an alternative 
  ## of distance heat map
  
  kmeans_mds_dist <- kmeans_clusters %>% 
    map(plot,
        type = 'components', 
        red_fun = 'mds') %>% 
    map(~.x + 
          theme(plot.tag = element_blank(), 
                legend.position = 'bottom'))

```
```r
> kmeans_mds_dist$training + kmeans_mds_dist$test
```

![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/501958ae-017c-41b4-a1c1-fcd1d51535b0)

```r

  ## UMAP of the data set
  
  kmeans_umap_data <- kmeans_clusters %>% 
    map(plot, 
        type = 'components', 
        red_fun = 'umap',
        with = 'data') %>% 
    map(~.x + 
          theme(plot.tag = element_blank(), 
                legend.position = 'bottom'))

```
```r
> kmeans_umap_data$training + kmeans_umap_data$test
```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/d105f4b6-8f66-4421-bbc4-8705ef513d2f)

By calling `plot_clust_hm()` for a clusterin analysis object, a heat map representation of clustering features is plotted. The color scale may be easily customized e.g. with ggplot's `scale_fill_gradient2()`:

```r
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
```
```r
> kmeans_hm_variables$training + kmeans_hm_variables$test
```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/12709b76-99eb-4106-8afb-0ed8caad70c2)
  
</details>

## Clustering variable importance

<details>
Permutation variable importance as proposed by Breiman for machine learning algorithms can be applied directly to clusterin analyses. The principle is quite simple: we're comparing quality fo clustering of the genuine clustering structure with a clustering objects fitted to the data set with one of the variables reshuffled by random. This procedure is repeated for all clusetring features and may run in several iterations exclude random unexpecte results. In `clustTools`, permutation variable importance is computed with the `impact()` function for explained clustering variance as loss function, i.e. the metric used to compare the input and reshuffled clustering structure. Of importance, this procedue can be done only for clustering objects fitted to the training data set and not for predictions. The computation for multiple iterations can be accelerated by launching he function in parallel. The results can be visualized by calling `plot()` and their wrpa-up retirieved by `summary()`:

```r

kmeans_variable_importance <-
    impact(kmeans_clusters$training,
           n_iter = 50,
           .parallel = TRUE)

```
```r
> summary(kmeans_variable_importance)
# A tibble: 9 × 9
  variable   mean      sd median    q25    q75     min    max n_iter
  <chr>     <dbl>   <dbl>  <dbl>  <dbl>  <dbl>   <dbl>  <dbl>  <int>
1 V1       0.0507 0.00502 0.0501 0.0481 0.0539 0.0386  0.0600     50
2 V2       0.0749 0.00897 0.0735 0.0700 0.0809 0.0468  0.0964     50
3 V3       0.0596 0.00790 0.0589 0.0538 0.0644 0.0461  0.0811     50
4 V4       0.0555 0.00867 0.0546 0.0490 0.0616 0.0333  0.0736     50
5 V5       0.0290 0.00567 0.0283 0.0257 0.0327 0.0156  0.0440     50
6 V6       0.157  0.00880 0.157  0.151  0.163  0.134   0.174      50
7 V7       0.0404 0.00516 0.0402 0.0365 0.0441 0.0307  0.0514     50
8 V8       0.0857 0.0107  0.0842 0.0785 0.0920 0.0621  0.109      50
9 V9       0.0166 0.00723 0.0154 0.0124 0.0219 0.00289 0.0339     50
```

```r
> plot(kmeans_variable_importance)
```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/2225235a-26df-4c8b-a4bf-3f6d4feddb25)

As inferred from the summary table and the plot, the variable `V6` is by far the most influential clustering variable.

</details>

## Self-organizing maps



## References

1. Murtagh F, Contreras P. Algorithms for hierarchical clustering: An overview. Wiley Interdiscip Rev Data Min Knowl Discov (2012) 2:86–97. doi:10.1002/widm.53
2. Hartigan JA, Wong MA. Algorithm AS 136: A K-Means Clustering Algorithm. Appl Stat (1979) 28:100. doi:10.2307/2346830
3. Schubert E, Rousseeuw PJ. Faster k-Medoids Clustering: Improving the PAM, CLARA, and CLARANS Algorithms. in Lecture Notes in Computer Science (including subseries Lecture Notes in Artificial Intelligence and Lecture Notes in Bioinformatics) (Springer), 171–187. doi:10.1007/978-3-030-32047-8_16
4. Kassambara A, Mundt F. factoextra: Extract and Visualize the Results of Multivariate Data Analyses. (2020) Available at: https://cran.r-project.org/web/packages/factoextra/index.html [Accessed May 14, 2022]
5. Rousseeuw PJ. Silhouettes: A graphical aid to the interpretation and validation of cluster analysis. J Comput Appl Math (1987) 20:53–65. doi:10.1016/0377-0427(87)90125-7
6. Lange T, Roth V, Braun ML, Buhmann JM. Stability-based validation of clustering solutions. Neural Comput (2004) 16:1299–1323. doi:10.1162/089976604773717621
7. Leng M, Wang J, Cheng J, Zhou H, Chen X. Adaptive semi-supervised clustering algorithm with label propagation. J Softw Eng (2014) 8:14–22. doi:10.3923/jse.2014.14.22
8. Kuhn M. Building predictive models in R using the caret package. J Stat Softw (2008) 28:1–26. doi:10.18637/jss.v028.i05
9. Croux C, Filzmoser P, Oliveira MR. Algorithms for Projection-Pursuit robust principal component analysis. Chemom Intell Lab Syst (2007) 87:218–225. doi:10.1016/j.chemolab.2007.01.004
10. McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. (2018) Available at: https://arxiv.org/abs/1802.03426v3 [Accessed February 21, 2022]
11. Konopka T. umap: Uniform Manifold Approximation and Projection. (2022) Available at: https://cran.r-project.org/web/packages/umap/index.html [Accessed June 1, 2022]
12. Breiman L. Random forests. Mach Learn (2001) 45:5–32. doi:10.1023/A:1010933404324
