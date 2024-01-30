[![R](https://github.com/PiotrTymoszuk/clustTools/actions/workflows/r.yml/badge.svg)](https://github.com/PiotrTymoszuk/clustTools/actions/workflows/r.yml)

<img src="https://github.com/PiotrTymoszuk/clustTools/assets/80723424/fcf39384-39f7-43f8-a396-b546a3218398" width="20%" height="20%" align = "right">

# clustTools
Comprehensive dimensionality reduction and cluster analysis toolset

## Description

The `clustTools` package provides a medley of functions used for seemless integration of various dimensionality reduction methods (multi-dimensional scaling/MDS, principal component analysis/PCA, uniform manifold approximation and projection/UMAP or factor analysis/FA), clustering (hierarchical clustering, K-means, PAM and density DBSCAN clustering) and self-orgnizing map (SOM) analyses. In addition, a set of functions is provided for visualization, quality control and cross-validation of the clustering results and semi-supervised clustering.

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

`clustTools` represents a one-stop shop for construction, diagnostic and validation of clustering solutions, which are scattered between many excellent R packages. Let's take a look at its basic functionalities by semi-supervised clustering of the `biopsy` data set included in the R's `MASS` package. This portion of data stores results of a clinical study on breast lesion biopsies, which are classified as benign or malignant based on 9 morphological and cytological variables assessed by a pathologist on a 1 - 10 scale each. The `biopsy` data set will be split into a training and test portion. I'll use the training portion for choice of the best performing clustering algorithm and cluster number. The test portion will be used for the final validation of the clustering results. The algorithms of interest are Ward's hierarchical clustering (internally implemented by `stats`), kmeans (`stats` package) and PAM (partitioning around medoids, package `cluster`). 

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
I'm checking first the clustering tendency of the training and test data subsets with the [Hopkins statistic](https://en.wikipedia.org/wiki/Hopkins_statistic) computed with the `get_clust_tendency()` function implemented in the package `factoextra`. Of note, Hopkins statistic ranges from 0 (completely homogenous distribution) to 1 (highly clustered distribution). Values around 0.5 indicate a random normal-like ditribution. With the Hopkins statistic values of approximately 0.7, the `biopsy` data set demonstrates a weak-to-moderate clustering tendency:

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

I'm using the `clustTools` functions `hcluster()` and `kcluster()` to cluster the training portion of the data set with hierarchical clusetring, kmeans and PAM. I'm also testing two distance measures: Euclidean and cosine. The initially guessed cluster number is set to `k = 2` as expected for the data set containing two types of samples. Importantly, those clustering functions take numeric data frames or matrices as the first argument and return objects of `clust_analysis` class, which can be, irrespective of the clustering algorithm, subjected to consistent quality control and validation.

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

The package offers a visual helper for verifying the initial cluster number guess. By calling `plot(type = 'diagnostic')` for a `clust_analysis` object, plots of dendrograms (hierarchical clustering only), within-cluster sum of squares and mean silhouette statistic values for clustering structures with varying cluster numbers can retrieved. All plots generated with the package's tools are `ggplot` graphical objects and can be easily customized by the user:

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

As indicated by the two main branches of the dendrogram, the bend of the within-cluster sum of square curve, and the peak of the mean [silhouette statistic](https://en.wikipedia.org/wiki/Silhouette_(clustering)), the two-cluster solution seem resonable.

</details>

### Quality control of clustering solutions

<details>

The package offers basically two numeric measures of cluster quality:

* [silhouette statistic](https://en.wikipedia.org/wiki/Silhouette_(clustering)), which gauges discrimination between the clusters

* _explained clustering variance_ defined as a ratio of the total between-cluster sum of squares to the total sum of squares. As such, explained clustering variance works in a similar way to R^2 or ANOVA's $\eta^2$ as a metric of explanatory performance

The `silhouette()` function computes silhouette values for each observation. Mean silhouette values for the whole object and particular clusters can be retrieved with the `summary()` method. Generally, silhouettes range from -1 to 1, with high values indicative of good separation of the given observation or cluster from other clusers. By contrast, negative silhouette widths suggest missclassified observations:

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

In particular, 18% of the cluster 2 observations had negative silhouette values, which indicates that they are most likely placed in a wrong cluster. 

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
In the `biopsy` data set, Euclidean distance yields a generally better separation between the clusters than the cosine metics as unequivocally demonstrated by mean silhouette value. The PAM/Euclidean solution has the lowest count of potentially misclassified observations with negative silhouette widths.

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

### Cross-validation

<details>

Unsupervised clustering objects resemble multi-parameter classification machine learning models in multiple aspects. By principle, they can be cross-validated just like any other machine learning model. In brief, in each cross-validation fold, the training portion is used for fitting of the clustering structure with the same algorithm as the parent clustering solution. The test portion of the fold is used for calculation of numeric stats and assessment of cluster assignment accuracy as compared with the parental clustering solution. Such approach has few advantages as compared with computation of numeric statistics only for the training data subset. for instance it is more robust at handling atypical observations and assessment of over-fitting. 
The default method of assignment of the test portion observations to the clusters is the so-called k-nearest neighbor (kNN) label propagation. This prediction algorithm is resonably fast and agnostic to the clustering procedure. 

Here, I'll cross-validate the `biopsy` clustering objects generated with the best performing Euclidean distance. This is accomplished with the `cv()` function, which returns out-of-fold predictions, numeric statistics for the folds along with means and confidence intervals of cluster assignment accuracy, classification error, explained clustering variance and mean silhouette. The means of statistics across all folds can be retrieved from the `cv()` output by calling `summary()`:

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

In case of kmeans clustering and other algorithms with stochastic determination of the initial cluster centers, the accuracy and classification error are expected to be poor. Still, we can use fractions of explained variations (`frac_var`) and mean silhouette widths (`sil_width`), to compare the clustering procedures. While all of them have similar discriminatory power, the kmeans/Euclidean colution displays the largest fraction of explained variance.
  
</details>

### Semi-supervised clustering

<details>

As for cross-validation, there's a possibility to apply the kNN label propagation algorithm for prediction of the cluster assignment in a hold-out subset or external validation data set. Here we'll predict the cluster assignment defined by the kmeans/Euclidean clustering in the training portion for the test subset of the `biopsy` data set and check quality of the cluster assignment with explained clustering variance and silhouette width. 
The prediction is done with the `predict()` function, which takes a clustering analysis object as the first argument and a numeric data frame or matrix as `newdata`; the `type` argument is set to 'propagation' for the kNN procedure. The output is an instance of `clust_analysis` class, which allows for a direct comparison of the training and test data clustering.

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

Finally, let's compare the distribution of malignant and benign samples between the clusters. This analysis reveals, that nearly all benign samples were assigned to the cluster 1 and malignant samples were located in the cluster 2:

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
We can corroborate it further by computing accuracy, Cohen's kappa and receiver-operating characteristic with `multiClassSummary()` provided by the `caret` package:

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

### Advanced visualization options

<details>

Visualization of a clustering structure with help of dimensionality reduction methods has a long tradition. By calling `plot(type = 'components')` for a clustering analysis objects, the `clustTools` package generates a scatter plot of observations's scores associated with the first two dimensions of prlincipal component analysis (PCA), multi-dimensional scaling (MDS) or Uniform Manifold Approximation and Projection (UMAP). By default, the dimensionality reduction is done for the distance matrix. By specifying `with = 'data'`, the user can request it for the genuine data set.

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

By calling `plot_clust_hm()` for a clustering analysis object, a heat map representation of the clustering features is plotted. The color scale may be easily customized e.g. with ggplot's `scale_fill_gradient2()`:

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

### Clustering variable importance

<details>
Permutation variable importance as proposed by Breiman for machine learning algorithms can be applied directly to clustering analyses. The principle is quite simple: we're comparing quality of the genuine clustering structure with a clustering object fitted to the data set with one of the variables reshuffled by random. This procedure is repeated for all clusetring features and may run in several iterations to exclude random unexpected results. In `clustTools`, permutation variable importance is computed with the `impact()` function for explained clustering variance as the loss function, i.e. the metric used to compare the input and reshuffled clustering structure. Of importance, this procedue can be done only for clustering objects fitted to the training data set and not for predictions. The computation in multiple iterations can be accelerated by launching he function in parallel. The results can be visualized by calling `plot()` and their wrap-up retrieved by `summary()`:

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

As inferred from the summary table and the plot, `V6` is by far the most influential clustering variable.

</details>

### Density clustering

<details>

Density clustering is implemented in the `clustTools` package by the DBSCAN algorithm. While its applicability to high dimension data is limited, it can robustly cope with data subjected previously to dimensionality reduction e.g. with PCA or UMAP. To streamline the reduction - clustering workflow, nearly all clustering functions of the package can be coupled with the `reduce_data()` tool. In the following example, we'll subject the `biopsy` data set to UMAP and subsequently cluster the UMAP score layout with DBSCAN. Unlike for hierarchical clustering or mean/medoid clustering, the user does not have to specify the cluster number. The clustering behavior is controlled by the `eps` parameter defining the distance threshold for noise observations and `minPts` specifying the size of the nearest neighborhood:

```r

  ## generating UMAP layouts for the training and test
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

```
In a similar manner to hierarchical or kmeans clustering, the function `plot(type = 'diagnostic')` may be employed for setting the optimal `eps` value. Following the logics of Belyadi and colleagues, noisy observations tend towards rapidly expanding k-nearest neighbor distance. Hence the optimal `eps` threshold should be placed just below the steep increase of the nearest neighbor distance plot. With `eps = 1`, we set the threshold obviously too high:

```r
> plot(biopsy_dbscan$training)
$knn_dist
```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/0f1ecd4a-1364-4da5-b71c-03c07ae36ab4)


`eps` value of 0.7 seems to do the job

```r
  biopsy_dbscan$training <-
    dbscan_cluster(data = biopsy_density$train,
                   distance_method = 'manhattan', 
                   eps = 0.7,
                   minPts = 7)

```

```r
> plot(biopsy_dbscan$training)
$knn_dist
```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/175191fb-039e-4599-8299-23f36f1c2b27)

Numeric statistics, cross-validation and semi-supervised clustering are accomplished in the usual way. Let's check the cluster predictions for the test portion of the `biopsy` data set:

```r
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

```

```r
>   biopsy_dbscan_variance
 training      test 
0.9922007 0.2608777

>   biopsy_dbscan_silhouette 
$training
# A tibble: 5 × 13
  clust_id     n n_negative perc_negative   mean       sd median   q025    q25   q75  q975    min   max
  <fct>    <int>      <int>         <dbl>  <dbl>    <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl>  <dbl> <dbl>
1 global     383         62          16.2 0.481   0.482    0.660 -0.507  0.236 0.938 0.959 -0.911 0.960
2 0            1          0           0   0      NA        0      0      0     0     0      0     0    
3 1          122          0           0   0.947   0.00912  0.947  0.928  0.939 0.955 0.960  0.926 0.960
4 2          177         62          35.0 0.0448  0.360    0.148 -0.583 -0.391 0.364 0.419 -0.911 0.420
5 3           83          0           0   0.733   0.0662   0.763  0.594  0.675 0.783 0.806  0.564 0.807

$test
# A tibble: 5 × 13
  clust_id     n n_negative perc_negative     mean      sd median   q025    q25    q75   q975     min      max
  <fct>    <int>      <int>         <dbl>    <dbl>   <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl>    <dbl>
1 global     300        146          48.7   0.0143  0.479   0.262 -0.730 -0.427  0.439  0.855  -0.785    0.885
2 0            0          0         NaN   NaN      NA      NA     NA     NA     NA     NA     Inf     -Inf    
3 1           13          0           0     0.843   0.0643  0.882  0.705  0.806  0.883  0.885   0.666    0.885
4 2          146        146         100    -0.453   0.159  -0.442 -0.759 -0.578 -0.325 -0.223  -0.785   -0.216
5 3          141          0           0     0.422   0.0751  0.430  0.262  0.380  0.461  0.529   0.228    0.529
```

Interestingly, the UMAP - DBSCAN procedure overfits massively as demonstrated by the huge differences in clustering variances and mean silhouettes between the training and test data set. The failure of the clustering solution to predict the test subset assignment is also obvius, when the UMAP layout is visualized in a scatter plot. To this end we'll call `plot(type = 'data')`, which plots the first two variables of the clustering data set - in this case the first two UMAP dimensions.

```r

  biopsy_dbscan_umap_layout <- biopsy_dbscan %>%
    map(plot, 
        type = 'data') %>% 
    map(~.x + 
          theme(plot.tag = element_blank(), 
                legend.position = 'bottom'))

```
```r
> biopsy_dbscan_umap_layout$training + biopsy_dbscan_umap_layout$test
```

![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/55ee42a3-d0d6-41f1-b173-07afb3fb1c64)

</details>

### Self-organizing maps

<details>
Self-organizing maps or SOM represent a Swiss army knife for dimensionality reduction and clustering. The algorithm based on the neural network principle was proposed by Teuvo Kohonen and implemented in R by the excellent package `kohonen` by Wehrens and colleagues. The output of the SOM algorithm is a set of vectors, so called 'codebook vectors' storing positions of SOM nodes (called also 'winning units') and, hence, a simplified, 'reduced' form of the input data set. As proposed by Juha Vesanto some 20 years ago, the matrix of distances between the codebook vectors, so called 'U matrix' may be tackled by traditional unsupervised clustering algorithms such as hierarchical or kmeans clustering. From my own experience, such combined SOM - unsupervised clustering procedure is especially useful at handling high dimension data such as gene expression matrices of flow cytometry measurmenents. 

The `clustTool` package integrates the SOM concept, SOM tools provided by the `kohonen` package and unsupervised clustering functions in a comprehensive workflow. In the current example, I'll use the combined SOM - PAM clustering to classify Italian wines listed in the UCI's `wines` data set based on their physical and chemical properties. The preprocessing will include elimination of non-variant features and normalization with median centering of the clustering variables (function `center_data()` from the `clustTools` package). The training subset will consist of 100 randomly selected observations, and will be used for tuning of the clustering algorithm. The test portion will be used solely for validation of the clustering structure.

```r
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

```

As before, we're checking the clustering tendency with Hopkins statistic using `get_clust_tendency()`. With the statistic values in the 0.65 - 0.7 range, only weak clustering tendency of both training and test data subsets can be inferred. 

```r
clust_tendency <- wine_lst %>%
    map(get_clust_tendency,
        n = 60)
```
```r
> clust_tendency$train$hopkins_stat
[1] 0.6983034
> clust_tendency$test$hopkins_stat
[1] 0.6566695
```
We are constructing three SOM-clusetring objects employing $5 \times 4$ hexagonal SOM grids to generate U matrices with Euclidean, Manhattan and cosine distances. The U-matrices will be subsequently clustered with the PAM procedure with 3 clusters as an intial guess. The construction step is accomplished by the `combi_cluster()` function returning a `combi_analysis` object. As we'll see in a moment, `combi_analysis` objects share the same quality control, validation, prediction, visualization and variable importance framework with other clustering analyses done with the `clustTools` package.

```r

algos <- list()

  algos$som_pam_euclidean <-
    combi_cluster(data = wine_lst$train,
                  distance_som = 'euclidean',
                  xdim = 5,
                  ydim = 4,
                  topo = 'hexagonal',
                  neighbourhood.fct = 'gaussian',
                  toroidal = FALSE,
                  rlen = 2000,
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
                  rlen = 2000,
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
                  rlen = 2000,
                  node_clust_fun = kcluster,
                  k = 3,
                  clust_fun = 'pam')

```
The first step of quality control is to check if the SOM converged. This can be easily done by plotting mean distances of the observations to the SOM node as a function of algorithm iterations. A substantial reduction of the distance followed by a plateau suggests convergence of the algorithm. Such distance plots are generated by `plot(type = 'training')`. As presented in the graphic below, SOM with all of Euclidean, Manhattan and cosine distances converged within 2000 iterations of the algorithm. If this is not the case, you may consider increasing the `rlen` parameter values in the `combi_cluster()` function.

```r

  ## the training plots are stored by the 'observation' element
  ## there are warnings generated that the training plots can not be generated
  ## for SOM nodes, but you can ignore them
  ##
  ## The plot() function returns a list of ggplot pbjects
  ## which are easy to customize by the user

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

```
```r
  som_convergence_plots$som_pam_euclidean +
    som_convergence_plots$som_pam_manhattan +
    som_convergence_plots$som_pam_cosine +
    plot_layout(ncol = 2)
```

![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/5a86e00c-af71-42a1-8456-b973e79ba50d)

In the next step, we're verifying, if our cluster number guess holds by comparing mean slhouette widths for clustering of the SOM nodes with varying cluster numbers. This can be done by calling `plot(type = 'diagnostic')` for a `combi_analysis` object. As for the remaining clustering analysis types, this function returns a list of diagnostic plots for the observation clustering by SOM (element 'observation') and for the SOM node clustering by PAM (element 'node' of the list). For deetrmination of the cluster number for node clustering, the node diagnostic plots are needed.

```r
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
```

```r
node_cluster_number_plots$som_pam_euclidean +
    node_cluster_number_plots$som_pam_manhattan +
    node_cluster_number_plots$som_pam_cosine +
    plot_layout(ncol = 2)
```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/0bf99fe0-d3a7-4292-9151-3d3290918020)

In the graphic panel above, the red dashed lines indicate the guessed number of clusters and the blue ones represent the optimal cluster number based on the peak of the silhouette statistic. With the cosine distance, our cluster number guess seems to be resonable. In turn, k = 2 clusters will likely work better for the clustering solutions with Euclidean and Manhattan distances. 

In order to select the algorithm with the largest explanatory value and the best separation between the clusters, we'll compare the explained clustering variance and mean silhouette statistic in 5-fold cross-validation. Generally, cross-validation for `combi_analysis` objects is performed in a similar way as for 'simple' clustering analyses, by calling the `cv()` function. Yet, because SOM is principally a neuronal network, we can use this trained structure to predict the node assignment for new data. Hence, we'll set the `type` argument to `'som'`:

```r
  ## SOM prediction method is recommended

  algos_cv <- algos %>%
    map(cv, type = 'som') %>%
    map(summary) %>%
    map(select, ends_with('mean'))
```
```r
> algos_cv
$som_pam_euclidean
# A tibble: 1 × 4
  accuracy_mean error_mean frac_var_mean sil_width_mean
          <dbl>      <dbl>         <dbl>          <dbl>
1         0.402      0.598         0.561          0.247

$som_pam_manhattan
# A tibble: 1 × 4
  accuracy_mean error_mean frac_var_mean sil_width_mean
          <dbl>      <dbl>         <dbl>          <dbl>
1         0.116      0.884         0.630          0.289

$som_pam_cosine
# A tibble: 1 × 4
  accuracy_mean error_mean frac_var_mean sil_width_mean
          <dbl>      <dbl>         <dbl>          <dbl>
1        0.0591      0.941         0.708          0.440

```
Similar to other clustering algorithms using an intial random placement of cluster centers or nodes, there's little concordance of the cluster assignment between the genuine clustering structure and the cluster assignment in the cross-validation folds. This results in a poor accuracy and large classification error. Still, the explained clustering variance and mean silhouette width may be used to compare the algorithms. In such comparison, the SOM/PAM/cosine procedure is characterized by the largest fraction of explained variance and the best cluster separation measured by silhouette width. 

The SOM/PAM/cosine clustering solution will be subsequently used to for semi-supervised clustering of the test portion of the `wines` data set. As above, the prediction will be done with the trained SOM network with `predict(type = 'som')`:

```r

  ## working with the best performing SOM/PAM/cosine algorithm

  cosine_clusters <- list()

  cosine_clusters$train <- algos$som_pam_cosine

  ## prediction of the cluster assignment for the test
  ## subset of the wines data set. Using the recommended
  ## SOM prediction method

  cosine_clusters$test <- predict(cosine_clusters$train,
                                  newdata = wine_lst$test,
                                  type = 'som')

```
The `predict()` function applied to `combi_analysis` objects returns `clust_analysis` objects, for which clustering variance, silhouette and plots can be retrieved as presented before for 'simple' hierarchical, kmeans and PAM clustering. 
For instance, the cluster structures for the training and test subsets display quite comparable clustering variances and mean silhouette values, as well as the numbers of potentially misclassified observations with negative silhouette values:

```r
  ## comparison of variances and silhouette widths
  ## in the training and test data portions

  cosine_variance <- cosine_clusters %>%
    map(var) %>%
    map_dbl(~.x$frac_var)

  cosine_silhouettes <- cosine_clusters %>%
    map(silhouette) %>%
    map(summary)
```

```r
> cosine_variance
    train      test 
0.6981659 0.7477879

> cosine_silhouettes
$train
# A tibble: 4 × 13
  clust_id     n n_negative perc_negative  mean    sd median    q025   q25   q75  q975    min   max
  <fct>    <int>      <int>         <dbl> <dbl> <dbl>  <dbl>   <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl>
1 global     102          5          4.90 0.476 0.237  0.560 -0.148  0.355 0.655 0.737 -0.263 0.747
2 1           42          2          4.76 0.517 0.231  0.612 -0.0248 0.388 0.685 0.745 -0.263 0.747
3 2           23          1          4.35 0.335 0.207  0.384 -0.0614 0.176 0.505 0.622 -0.137 0.647
4 3           37          2          5.41 0.516 0.232  0.602 -0.167  0.492 0.647 0.723 -0.229 0.734

$test
# A tibble: 4 × 13
  clust_id     n n_negative perc_negative  mean    sd median    q025   q25   q75  q975    min   max
  <fct>    <int>      <int>         <dbl> <dbl> <dbl>  <dbl>   <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl>
1 global      74          2          2.70 0.522 0.233  0.580 -0.0171 0.376 0.703 0.759 -0.284 0.766
2 1           27          0          0    0.610 0.155  0.676  0.273  0.545 0.722 0.764  0.216 0.766
3 2           21          2          9.52 0.296 0.241  0.372 -0.284  0.181 0.433 0.559 -0.284 0.562
4 3           26          0          0    0.613 0.164  0.695  0.270  0.549 0.731 0.743  0.209 0.751
```
We can easily visualize pairwise distances between the observations and SOM nodes (i.e. the U matrix) in form of heat maps by calling `plot(type = 'heat_map')`:

```r
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
```

```r
  cosine_umatrix_hm +
    cosine_test_hm
```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/86006d9b-8036-4531-b8be-ae03ab72a66d)

`plot()` can be also employed to obtain plots of UMAP layouts of the training and test data subsets:

```r
  cosine_umap_plots <- cosine_clusters %>%
    map(plot,
        type = 'components',
        kdim = 2,
        with = 'data',
        red_fun = 'umap')

  ## for the training data subset:
  ## a list of UMAP layout plots is returned, we need only
  ## the UMAP layout for the observations with color coding
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
```
```r
  cosine_umap_plots$train +
    cosine_umap_plots$test
```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/e2e3e989-d897-4e71-9264-e82e04395b53)

Heat map representations of the clustering variable levels can be obtained with `plot_clust_hm()`:

```r

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

```

```r
  cosine_feature_hm$train + 
    cosine_feature_hm$test
```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/02e02181-51d5-43bf-842c-ef35cdb22b71)

By calling `impact()`, the user can compute permutation variable importance the same way as for simpler clustering analyses - but only for the training clustering structure. In such analysis, content of flavonoids, proanthocyanins and proline, total phenol level and optical density ratio are indetified as the most important clustering variables for explanatory performance of the algorithm:

```r
  cosine_importance <- impact.combi_analysis(cosine_clusters$train,
                                             n_iter = 50,
                                             .parallel = TRUE)

  plot(cosine_importance)

```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/7ee0a38c-18c1-4ccc-9462-43246d3f7d9b)

In the final analysis step of the `wine` data, we'll check for vintage classes in the clusters. As presented below, the cluster 1 consist primarily of Barolo wines, the cluster 2 includes almost exclusively Grignolino. Barbera and, to a lesser extent, Grignolino wines populate the cluster 3. 

```r
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
```
```r
> cosine_counts
$train
# A tibble: 5 × 3
  clust_id vintage        n
  <fct>    <fct>      <int>
1 1        Barolo        33
2 1        Grignolino     9
3 2        Grignolino    23
4 3        Barbera       28
5 3        Grignolino     9

$test
# A tibble: 7 × 3
  clust_id vintage        n
  <fct>    <fct>      <int>
1 1        Barolo        23
2 1        Grignolino     4
3 2        Barolo         2
4 2        Grignolino    19
5 3        Barbera       20
6 3        Grignolino     6
7 NA       Grignolino     1
```
Interestingly, one Grignolino sample in the training subset cuold not be assigned to any cluster defined in the training portion of the data.

Let's compare the observed vintage assignment with the predominant cluster vintages in a more formal inter-rater and reveiver operating characteristic analysis done with caret's `multiClassSummary()`:

```r
  ## kappa and ROC analysis
  ## renaming of the clusters after
  ## the predominant wine type

  cosine_assignment <- cosine_assignment %>%
    map(mutate,
        obs = vintage,
        pred = car::recode(clust_id,
                           "'1' = 'Barolo';
                           '2' = 'Grignolino';
                           '3' = 'Barbera'"),
        pred = factor(pred,
                      levels = c('Barbera', 'Barolo', 'Grignolino')))

  cosine_roc <- cosine_assignment %>%
    map(as.data.frame) %>%
    map(multiClassSummary,
        lev = c('Barbera', 'Barolo', 'Grignolino'))
```
```r
> cosine_roc
$train
              Accuracy                  Kappa                Mean_F1       Mean_Sensitivity       Mean_Specificity 
             0.8235294              0.7391675              0.8200962              0.8536585              0.9159812 
   Mean_Pos_Pred_Value    Mean_Neg_Pred_Value         Mean_Precision            Mean_Recall    Mean_Detection_Rate 
             0.8474903              0.9240506              0.8474903              0.8536585              0.2745098 
Mean_Balanced_Accuracy 
             0.8848199 

$test
              Accuracy                  Kappa                Mean_F1       Mean_Sensitivity       Mean_Specificity 
             0.8378378              0.7581699              0.8380602              0.8583908              0.9209373 
   Mean_Pos_Pred_Value    Mean_Neg_Pred_Value         Mean_Precision            Mean_Recall    Mean_Detection_Rate 
             0.8419482              0.9229225              0.8419482              0.8583908              0.2792793 
Mean_Balanced_Accuracy 
             0.8896640 
```
  
</details>

### Handling multi-layer data with self-organizing maps

<details>
As discussed in mode details by the authors of the `kohonen` package (see Wehrens et al. in the reference list), SOM are well suited to handle multi-layer data sets. Such layers may consist of variables in various formats (e.g. binary and ordinal layer) or features concerning diverse properties of the investigated objects (e.g. position and RGB color of pixels in a picture). In the following example we will apply this SOM-based multi-layer approach to cluster the `Cars93` data set provided by the R package `MASS`. In brief, the data layers will consist of (1) numeric features such as horsepower or weight, (2) ordinal/count variables like number of passengers or airbags, and (3) binary such as manual transmission or US origin. Since the data set is pretty small, we won't validate the results in a 'clean' hold-out setting and stay with cross-validation.

```r

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

```
The data set has a moderate clustering tendency as investigated by `get_clust_tendency()`. 

Fitting of SOM with subsequent clustering of the U matrix, i.e. matrix of distances between the SOM nodes is done with a dedicated functon of the `clustTools` package: `multi_cluster()`. Its syntax is quite similar to `combi_cluster()`. Note, that the function allows for definition of separate distance measures for each data layer via the `distance_method` argument. In our case, the numeric layer will be clustered with cosine, the ordinal layer with Manhattan and the binary layer with Tanimoto distance. Additional arguments that control the SOM behavior are passed by a named list via the `som_args` argument: here, we are defining a custom learning rate `alpha`, weights for the layers via `user.weights` and telling the SOM algorithm to apply the user's layer weights directly to the data layers by specifying `normalizeDataLayers = FALSE`. Clustering of the U matrix is done with PAM implemented by the `kcluster()` function.

```r

  ## the layers will be handled by three various distances:
  ## cosine for numeric variables, Manhattan for ordinal variables
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

```
As demonstrated before, we're checking convergence of the SOM algorithm by calling `plot(type = 'training')`. We're customizing the training plot to present the training curve for each layer in a separate facet:

```r
  car_convergence <- plot(car_som, type = 'training')

  car_convergence <- car_convergence$observation +
    facet_grid(`SOM layer` ~ .,
               scales = 'free') +
    theme(plot.tag = element_blank(),
          legend.position = 'none')

  car_convergence
```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/0ac5ef8c-2ddb-4fe0-b623-e27e251e4c06)


Based on plots of within-cluster sum of squares and silhouette widths for varying number of clusters of SOM nodes (generated by `plot(type = 'diagnostic')`), 2 - 3 clusters are proposed for the optimal solution. Hence our initial guess (k = 3) turned to be quite correct. 

```r

 ## WSS and silhouette plots

  car_clust_number <- plot(car_som, type = 'diagnostic')$node %>%
    map(~.x + theme(plot.tag = element_blank()))

  car_clust_number$wss +
    car_clust_number$silhouette

```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/862639b5-9e7d-44c3-90c6-2e24bbb87b18)

Cross-validated explained clustering variance and silhouette width will be obtained by calling `cv(type = 'som')` for the multi-layer cluster analysis object:

```r
  car_cross_validation <-
    cv(car_som,
       nfolds = 10,
       type = 'som') %>%
    summary

  car_cross_validation %>%
    select(ends_with('mean'))
```
```r
# A tibble: 1 × 4
  accuracy_mean error_mean frac_var_mean sil_width_mean
          <dbl>      <dbl>         <dbl>          <dbl>
1         0.123      0.877         0.701          0.398
```
Again, although the overal accuracy is quite poor because of partly stochastic process of SOM fitting, the clustering structure could explain more than 70% of the `Cars93` data set variance. The cross-validated mean silhouette width suggests a moderate degree of cluster separation. Let's take a more detailed look at separation of particular clusters by plotting silhouette widths for observations and a heat map of the U matrix. Those plots indicate possibly problematic separation of the clusters 1 and 3:

```r

  ## silhouette values for particular clusters

  car_sil_plot <- car_som %>%
    silhouette %>%
    plot

  ## heat map of pairwise distances

  car_heat_map <- 
    plot(car_som, type = 'heat_map')$node + 
    labs(title = 'Distance between SOM nodes') + 
    theme(plot.tag = element_blank())

```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/382091b5-9f67-44e5-8fee-816f5cc4918a)

![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/1d59c272-7c27-48da-bd8e-7a48f981b5cb)

To get some insight into levels of the cluster-defining variables, `plot_clust_hm()` is used. This function called for multi-layer analysis objects returns, unlike for single-layer analyses, a list of heat maps representing the data layers:

```r
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
```
```r
  car_feature_hm$numeric + 
    car_feature_hm$ordinal + 
    car_feature_hm$binary + 
    plot_layout(ncol = 2)
```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/1a776fa6-6238-4870-9494-cbac7fe739bd)

This graphic reveals that the largest cluster 2 groups low-motorized, small and economical cars with low airbag numbers, and designed for few passengers. By contrast, cluster 1 vehicles are expensive, mostly US origin, large and frequently equipped with automatic transmission and offer substantially more passenger seats. Cluster 3 cars are powerfully motorized, expensive, middle sized and designed for very few passengers. With those cluster characteristic in mind, we'll investigate the 'true' car type frequency within the clusters.

```r
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
```

```r
> car_counts
# A tibble: 10 × 4
   clust_id Type        n percent
   <fct>    <fct>   <int>   <dbl>
 1 1        Compact     2    6.06
 2 1        Large      11   33.3 
 3 1        Midsize    11   33.3 
 4 1        Van         9   27.3 
 5 2        Compact    14   28.6 
 6 2        Midsize     6   12.2 
 7 2        Small      21   42.9 
 8 2        Sporty      8   16.3 
 9 3        Midsize     5   45.5 
10 3        Sporty      6   54.5 

```
```r
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

  car_count_plot
```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/612639c7-fa5a-421d-822b-3fa15f6928e2)

In accordance with the clustering features, the cluster 2 consists predominantly of small and compact cars, while vans, midsize and large vehicles populate the cluster 2. The remaining cluster 3 consists of middle size and sporty models. 

Finally, let's investigate impact of particular clustering variables on the algorithm's explanatory performance by computing permutation variable importance with `impact()`:

```r
  car_importance <- impact(car_som,
                           n_iter = 50,
                           .parallel = TRUE)

  plot(car_importance) + 
    geom_vline(xintercept = 0, 
               linetype = 'dashed')
```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/55a7150f-7d4c-4b3d-ab74-954e14d86767)

Those results suggest cylinder number as by far the most influential clustering feature. Quite unexpectedly, the binary variables manual transmission and origin were identified as the least important variables. This can be however explained by setting the `user.weights` parameter during construction of the multi-layer analysis object. With this manual intervention, the binary data layer was weighted lower than numeric and ordinal features.
  
</details>

### Regularized clustering

<details>

At the moment, the `clustTools` package implements one algorithm of regularized unsupervised clustering, __the Hard-Threshold KMEANS or HTK clustering__, introduced recently by Raymaekers and Zamar. Here, I'll demonstrate its usage and potential advantages at dealing with sparse data with the `wines` data set provided in the `MASS` package.

```r

  ## required libraries and de-masking of some functions

  library(kohonen)

  library(clustTools)
  library(tidyverse)
  library(clusterHD)
  library(caret)
  library(somKernels)

  extract <- clustTools::extract
  map <- purrr::map
  library(patchwork)

```

As before, we'll start with selection of variables with a threshold of the mean-to-variance ratio statistic, normalize the clustering features and create the training and test subsets. For simplicity, I'm skippng an analysis of distribution of clustering variables.

```r
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
```

Next, we are going to train clustering structures with the canonical KMEANS and the regularized HTK algorithms. Note: for the KMEANS clustering structure, squared Euclidean distance is used, the same as implemented in the HTK algorithm by default. By this means, we can compare the two clustering procedures in a head-to-head manner. For the HTK clustering structure, choice of the penalty parameter `lambda` is of critical importance. Cross-validation tuning popularized by supervised machine learning is a handy solution to this problem. By calling `tune_htk()` with the argument specified below, we're looking for the optimal `lambda` for a k = 3 cluster solution with the maximum out-of-fold silhouette width in 10-fold cross-validation. The argument `kNN = 11` specifies the number of nearest neighbors for the nearest neighbor classifier used to assign CV fold observations to the clusters.

```r

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

```
By calling `summary()` for the tuning object, the full results of tuning can be accessed:

```r
>   summary(train_tune)
# A tibble: 41 × 7
   lambda object_id sil_width frac_misclassified frac_var frac_np n_active_vars
    <dbl> <chr>         <dbl>              <dbl>    <dbl>   <dbl>         <dbl>
 1  0     obj_1         0.382              0.151    0.667   0.453             9
 2  0.025 obj_2         0.382              0.151    0.667   0.453             9
 3  0.05  obj_3         0.382              0.151    0.667   0.453             9
 4  0.075 obj_4         0.382              0.151    0.667   0.453             9
 5  0.1   obj_5         0.382              0.151    0.667   0.453             9
 6  0.125 obj_6         0.382              0.151    0.667   0.453             9
 7  0.15  obj_7         0.382              0.151    0.667   0.453             9
 8  0.175 obj_8         0.382              0.151    0.667   0.453             8
 9  0.2   obj_9         0.382              0.151    0.667   0.453             8
10  0.225 obj_10        0.382              0.151    0.667   0.453             8
```
Quality control graphics such as clusetr performance statistics for subsequent lambda values and regularization paths can be retrieved by `plot()`. A short visual analysis of the output suggests that the HTK algorithm works best for `lambda = 0.18`, which results in elimination of one variable: the `ash alkalinity`:

```r
>   plot(train_tune)
$statistics

$regularization
```

![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/3c6c5fbd-119b-44b6-b89f-8d1087dadeaf)

![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/7ff210e5-51b0-491f-b284-aff40755e8be)

Next, we can get the best clustering structure wrapped in a `clust_analysis` object by calling `extract()` for the tuning result. I'm also renaming the clusters by their predominant vintage checked elswhere:

```r

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

```
As usual, performance statistics such as mean silhouette width, fraction of candidate misclassified observations (i.e. observations with negative silhouette width), fraction of explained clustering variance and nearest neighborhood preservation is returned by `summary()` for the KMEANS and the optimal HTK clustering solution in the training data portion. Of note, the HTK solution is characterized by better cluster separation and supreme neighborhood preservation. In turn, the explained clustering variance is worse, but we must keep in mind, that the HTK algorithm works with one variable less as compared with the canonical KMEANS:

```r
>   train_clusters %>%
+     map(summary)
$kmeans
# A tibble: 1 × 4
  sil_width frac_misclassified frac_var frac_np
      <dbl>              <dbl>    <dbl>   <dbl>
1     0.399             0.0686    0.626   0.904

$htk
# A tibble: 1 × 4
  sil_width frac_misclassified frac_var frac_np
      <dbl>              <dbl>    <dbl>   <dbl>
1     0.440             0.0196    0.572   0.922

```

We can verify it in the test portion of the wine data set. Prediction of the cluster assignment in the test subset is done by a nearest neighbor classifier which works in the same manner for the classical KMEANS and the HTK clustering:

```r
  test_clusters <- train_clusters %>%
    map(predict.clust_analysis,
        newdata = wine_lst$test,
        type = 'propagation',
        kNN = 11,
        active_variables = TRUE)
```
Also in the test subset, the HTK clustering solution beats the KMEANS algorithm in terms of cluster seperation and neighborhood preservation:

```r
>   test_clusters %>%
+     map(summary)
$kmeans
# A tibble: 1 × 4
  sil_width frac_misclassified frac_var frac_np
      <dbl>              <dbl>    <dbl>   <dbl>
1     0.391              0.107    0.604    0.84

$htk
# A tibble: 1 × 4
  sil_width frac_misclassified frac_var frac_np
      <dbl>              <dbl>    <dbl>   <dbl>
1     0.459             0.0667    0.623   0.877
```

The classification performance of the HTK algorithm can be visualized in an even more impressive manner in UMAP plots. For sake of clarity, we're going to plot the cluster assignment with the same UMAP layout, which is accomplished by providind a common 'training' UMAP structure:


```r

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

```
```r
>   kmeans_plots$train_kmeans +
+     kmeans_plots$test_kmeans +
+     kmeans_plots$train_htk +
+     kmeans_plots$test_htk
```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/8e3f2bfa-4932-442e-b4cb-60f41e45f5a2)

Finally, let's compare the wine data set clusters with the observed vintage with a formal ROC analysis. This analysis shows clearly, that the HTK algorithm predicts the vintage with grater sensitivity, specificity and concordance (Cohen's $\kappa$) that the KMEANS solution:

```r

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

## and a common ROC summary used later for plotting

  roc_summary <- c(kmeans_roc, htk_roc) %>%
    reduce(rbind) %>%
    as_tibble %>%
    mutate(dataset = rep(c('training', 'test'), 2),
           algorithm = c(rep('KMEANS', 2),
                         rep('HTKmeans', 2)))

```
```r

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


```
![image](https://github.com/PiotrTymoszuk/clustTools/assets/80723424/2524ec0c-cca9-4a87-8d32-f5d5057eb32f)

A quick look at the confusion matrix, which is left to an interested user, indicates that the better performance of the HTK algoritm manifests prmarily by almost perfect identification of the Grignolino wines, i.e. the vintage which is recognized by KMEANS with pretty low accuracy.
  
</details>

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
13. Hahsler M, Piekenbrock M, Doran D. Dbscan: Fast density-based clustering with R. J Stat Softw (2019) 91:1–30. doi:10.18637/jss.v091.i01
14. Belyadi H, Haghighat A, Nguyen H, Guerin A-J. IOP Conference Series: Earth and Environmental Science Determination of Optimal Epsilon (Eps) Value on DBSCAN Algorithm to Clustering Data on Peatland Hotspots in Sumatra Related content EPS conference comes to London-EPS rewards quasiparticle research-EP. IOP Conf Ser Earth Environ Sci (2016) 31: doi:10.1088/1755-1315/31/1/012012
15. Kohonen T. Self-Organizing Maps. Berlin, Heidelberg: Springer Berlin Heidelberg (1995). doi:10.1007/978-3-642-97610-0
16. Vesanto J, Alhoniemi E. Clustering of the self-organizing map. IEEE Trans Neural Networks (2000) 11:586–600. doi:10.1109/72.846731
17. Wehrens R, Kruisselbrink J. Flexible self-organizing maps in kohonen 3.0. J Stat Softw (2018) 87:1–18. doi:10.18637/jss.v087.i07
18. Raymaekers J, Zamar RH. Regularized K-means Through Hard-Thresholding. J Mach Learn Res (2022) 23:1–48. Available at: http://jmlr.org/papers/v23/21-0052.html [Accessed November 15, 2023]
