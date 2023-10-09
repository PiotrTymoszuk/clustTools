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

I'm using the `clustTools` functions `hcluster()` and `kcluster()` to cluster the training portion of the data set with the hierarchical clusetring, kmeans and PAM. I'm also trying two distance measures: Euclidean and cosine. The initially guessed cluster number is set to `k = 2` as expected for the data set containing two types of samples. Importantly, those clustering functions return objects of `clust_analysis` class, which can be, irrespectively of the clustering algorithm, subjected to a consistentquality control and validation.

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
In the `biopsy` data set, Euclidean distance yields a generally better separation between the clusters as the cosine metics as unequivocally demonstrated by mean silhouettes. The PAM/Euclidean solution has also the lowest count of potentially misclassified observations.

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

## References

1. Murtagh F, Contreras P. Algorithms for hierarchical clustering: An overview. Wiley Interdiscip Rev Data Min Knowl Discov (2012) 2:86–97. doi:10.1002/widm.53
2. Hartigan JA, Wong MA. Algorithm AS 136: A K-Means Clustering Algorithm. Appl Stat (1979) 28:100. doi:10.2307/2346830
3. Schubert E, Rousseeuw PJ. Faster k-Medoids Clustering: Improving the PAM, CLARA, and CLARANS Algorithms. in Lecture Notes in Computer Science (including subseries Lecture Notes in Artificial Intelligence and Lecture Notes in Bioinformatics) (Springer), 171–187. doi:10.1007/978-3-030-32047-8_16
4. Kassambara A, Mundt F. factoextra: Extract and Visualize the Results of Multivariate Data Analyses. (2020) Available at: https://cran.r-project.org/web/packages/factoextra/index.html [Accessed May 14, 2022]
5. Rousseeuw PJ. Silhouettes: A graphical aid to the interpretation and validation of cluster analysis. J Comput Appl Math (1987) 20:53–65. doi:10.1016/0377-0427(87)90125-7
