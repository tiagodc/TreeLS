# TreeLS

R functions for processing individual tree TLS point clouds.

# 2019 UPDATE
I haven't touched this project in 2 years. I have reimplemented some of its functions (way) more efficiently in C++. I will refactor the whole package on the next weeks, incorporating the heavy functionality through `Rcpp`, as well as make the package compatible with the LAS infrastructure kindly provided by @Jean-Romain with the [`rlas`](https://github.com/Jean-Romain/rlas) and [`lidR`](https://github.com/Jean-Romain/lidR) packages.

The current package should work just fine for the main functionalities, although some dependendencies might be broken. I'll open a development branch to work on the *new*, improved version of `TreeLS`, releasing it as soon as I have a first stable version of the package.

## Installation
`devtools` is required*

#### On the R console, enter: `devtools::install_github('tiagodc/TreeLS')`

## Example

Load the package in the R environment.
```
require(TreeLS)
```

Two sample tree point clouds are provided: `spruce` and `pine`. To vizualize them:
```
rgl.points(pine, size=.5, col=cloud.col(pine, n=20))
```

Classification of trunk points can be done by 3 different methods: `pref_HT`, `pref_SD` and `pref_VN`. The fastest method is usually `pref_HT`.

```
trunk = pref_HT(pine)

#plot the trunk points
rgl.points(trunk, col='green')

#plot the rest of the tree
rgl.points(pine, size=.5)
```

Stem modelling is done over the previously extracted trunk points. Again, 3 methods can be used: `fit_RANSAC_circle`, `fit_RANSAC_cylinder` and `fit_IRTLS`. The fastest method is usually `fit_RANSAC_circle`.

```
stem = fit_RANSAC_circle(trunk)

#plot the stem model
mod = stem.model(stem, cyl.len=.3)

#plot the rest of the tree
rgl.points(pine, size=.5)
```

The output of the stem function is a list with two matrices: `$stem.points`, containing the stem point cloud, and `$circles`, containing estimations and statistics of modelled stem segments.

```
head(stem$stem.points)
print(stem$circles)
```
