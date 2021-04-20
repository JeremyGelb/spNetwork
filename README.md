
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spNetwork <img src='man/figures/spNetwork_logo.png' align="right" height="138.5" />

## A R package to perform spatial analysis on networks.

## Getting Started

A good start point for this package is the vignettes. They present the
main features of the package

### Installing

you can install this package with the following code in R. The packages
use mainly the following packages in its internal structure :

-   igraph
-   sp
-   rgeos
-   maptools
-   raster
-   future
-   future.apply
-   dplyr
-   tidyR
-   Rcpp
-   RcppArmadillo

``` r
devtools::install_github("JeremyGelb/spNetwork")
```

### Examples

We provide here some short examples of the main features

-   realizing a kernel network density estimate

``` r
library(spNetwork)
data(mtl_network)
data(bike_accidents)
lixels <- lixelize_lines(mtl_network,200,mindist = 50)
samples <- lines_center(lixels)
densities <- nkde(mtl_network,
                 events = bike_accidents,
                 w = rep(1,nrow(bike_accidents)),
                 samples = samples,
                 kernel_name = "quartic",
                 bw = 300, div= "bw",
                 method = "discontinuous", digits = 1, tol = 1,
                 grid_shape = c(1,1),
                 verbose=FALSE)
```

-   Building a spatial matrix based on network distance and use it to
    calculate the Moran I with spdep.

``` r
library(spNetwork)
library(spdep)
data(mtl_network)

conv_function <- function(x){
  if(x<=500){
    return(1/x**2)
  }else{
    return(0)
  }
}

listw <- line_ext_listw_gridded(mtl_network,maxdistance=500,
       dist_func = conv_function, matrice_type='W',
       grid_shape = c(5,5),
       mindist = 10)

moran.test(mtl_network$nbAccident,listw, zero.policy = T)
```

Note that you can use this in every spatial analysis you would like to
perform. With the converter function of spdep (like listw2mat), you can
convert the listw object into regular matrix if needed

### Work in progress

Currently, functions to calculate K-function and cross-K-function for
sets of points are available but still experimental.

We also added functions for bandwidth selection by likelihood
cross-validation

## Authors

-   **Jeremy Gelb** - *Creator and maintainer*

## Contribute

To contribute to `spNetwork`, please follow these
[guidelines](https://github.com/JeremyGelb/spNetwork/blob/master/CONTRIBUTING.md).

Please note that the `spNetwork` project is released with a [Contributor
Code of
Conduct](https://github.com/JeremyGelb/spNetwork/blob/master/CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

## License

`spNetwork` version 0.1.1 is licensed under [GPL2
License](https://github.com/JeremyGelb/spNetwork/blob/master/LICENSE.txt).

## Acknowledgments

-   Hat tip to Philippe Apparicio for its support during the development
-   Hat tip to Hadley Wickham and its helpful book *R packages*
