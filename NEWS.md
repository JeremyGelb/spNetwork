# spNetwork 0.4.1.9000

## New Features

* It is possible now to choose between two strategies to deal with densities at 0 when calculating LOO likelihood for several bandwidths. For more details, see the parameter *zero_strat* in the functions *bw_cvl_calc*, *bw_cvl_calc.mc*, *bw_cv_likelihood_calc*, *bw_cv_likelihood_calc.mc*, *bws_tnkde_cv_likelihood_calc* and *bws_tnkde_cv_likelihood_calc.mc*.

* It is possible now to set directly the size of local bandwidths for the functions *nkde*, *nkde.mc*, *tnkde* and *tnkde.mc*. This feature is currently experimental, please submit an issue if you encounter a bug.

# spNetwork 0.4.0.9000

## New Features

* the package has been reworked to use sf rather than rgeos, maptools and sp. All the test passed and the documentation has been modified.

# spNetwork 0.3.0.9000

## New Features

* Added a set of functions to perform Temporal Network Kernel Density Estimation (see vignette TNKDE)

## documentation

* A new vignette is presenting the Temporal Network Kernel Density Estimation
* Adding a CITATION file

## performance

* Performance was improved for data driven bandwidth selection for NKDE
* The function *lines_points_along* has been reworked with c++ and is now way faster
* The function *lines_center* has been reworked with c++ and is now way faster

## corrected bugs

* An error was raised when points were sharing the exact same locations in the function *network_knn*. It has been corrected an a warning was added to clarify the situation.

* A clearer error message is now raised when multiple geometries are passed to the  *nkde* and *nkde.mc* functions.

* When using the discontinuous NKDE, NA values were returned if an event was located exactly at the extremity of an edge with no connexions. This has been corrected.

# spNetwork 0.2.1

Minor release to fix a bug on SOLARIS for CRAN

# spNetwork 0.2.0

## New Features

* Added two function for data-driven bandwidth selection (see vignette NKDE)
* Added a function to simplify network before processing (*simplify_network*), reducing network complexity by removing useless edges or merge continuous edge can improve performance. This function is still experimental and must be used only for undirected networks.
* Added two functions to calculate k nearest neighbours: *network_knn* and *network_knn.mc*
* Added a function to calculate isochrones on a geographical network: *calc_isochrones* and a Vignette (*Calculating isochrones*) to give an example.
* Making some internal functions visible for the user to allow graph analysis on networks, see the new vignette *Building networks*.

## documentation

Several vignettes were added to present the new features. A pkgdown website was also built to give a more attractive way to browse the documentation.

We follow now the code coverage. We want to reach 70% coverage before the next CRAN release.

## corrected bugs

* issue 1 fixed by changing the function maptools::snapPointsToLines for a home made function. The package BH including the c++ boost library is used to realize fast spatial query.
* issue 2 fixed by handling the special case of lines of length 0.
* an error was raised when using *listw_network* with polygons and the method = "pointsalong". This was caused by a minor error in the code and works now as expected.
* an error was present in the Diggle correction factor, leading to overestimation of the correction. This error has been corrected.
* for the k functions, a minor error in the calculus was corrected, affecting only the scale of the result. These functions are now formally tested.
* In the nkde function, when using an integer matrix, an error was crashing the R session in some case. This has been corrected.

## performance

* Performance was improved for functions *listw_network*, *listw_network.mc*, *kfunctions*, *kfunctions.mc*, *cross_kfunctions* and *cross_kfunctions.mc*, by reducing the complexity of the network with no loss of accuracy.
* Several geometrical functions have been accelerated by c++ code using the boost library
* upgraded c++ compilation requirement to c++17
* the functions used for the nkde calculation were reworked. Sensible speed gain can be observed.

# spNetwork 0.1.1

* Accepted CRAN version (minor changes)

# spNetwork 0.1.0
  
* Initial release version
