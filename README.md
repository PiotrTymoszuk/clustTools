[![R](https://github.com/PiotrTymoszuk/clustTools/actions/workflows/r.yml/badge.svg)](https://github.com/PiotrTymoszuk/clustTools/actions/workflows/r.yml)

# clustTools
Comprehensive dimensionality reduction and cluster analysis toolset

## Description

The `clustTools` package provides a medley of functions used for seemless integration of various dimensionality reduction methods (multi-dimensional scaling/MDS, principal component analysis/PCA, uniform manifold approximation and projection/UMAP or factor analysis/FA), clustering (hierarchical clustering, K-means, PAM and density DBSCAN clustering) and self-orgnizing map (SOM) analyses. In addition, a set of functions is provided for visualization, quality control and cross-validation of the clustering results and semi-superrvised clusetring.

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
