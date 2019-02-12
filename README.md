# TreeLS

High performance R functions for forest inventory based on **T**errestrial **L**aser **S**canning (but not only) point clouds.

## Description

This package is a refactor of the methods described in [this paper](https://www.researchgate.net/publication/321434623_Performance_of_stem_denoising_and_stem_modelling_algorithms_on_single_tree_point_clouds_from_terrestrial_laser_scanning).

The algorithms were rewritten in C++ and wrapped in R functions through `Rcpp`. The algorithms were reviewed and enhanced, new functionalities introduced and the rebuilt functions now work upon [`lidR`](https://github.com/Jean-Romain/lidR/)'s `LAS` objects infrastructure.

This is an ongoing project and new features will be introduced often. For any questions or comments please contact me through github. Suggestions, ideas and references of new algorithms are always welcome - as long as they fit into TreeLS' scope.

## Main functionalities
- Tree detection at plot level
- Stem points detection at single tree and plot levels

## Coming soon:
- RANSAC circle fitting
- RANSAC based individual tree modelling
- RANSAC plot-wise tree modelling
- `lidR` wrappers for writing TLS data with extra header fields
- Specific 3D visualization functions for TLS specific outputs
- Eigen decomposition feature detection for trees and stems
- Tree modelling based on robust cylinder fitting
- 3D interactive point cloud manipulation

## Installation

### Requirements
- `devtools`: run `install.packages('devtools', dependencies = TRUE)` from the R console
- Rcpp compiler:
    - on Windows: install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for your installed R version - make sure to add it to your system's *path*
    - on Mac: install Xcode
    - on Linux: be sure to have `r-base-dev` installed

### Install TreeLS development version

On the R console, run:
```
devtools::install_github('tiagodc/TreeLS', ref='devel')
```

#### Legacy code

For anyone still interested in the old implementations (fully developed in R code, slow but suitable for research) of this library, you can still use it. In order to do it, uninstall any recent instances of `TreeLS` and reinstall the legacy version:
```
devtools::install_github('tiagodc/TreeLS', ref='old')
```
Not all features from the old package were reimplemented in `Rcpp`, but I'm gonna get there eventually.

## Examples

...