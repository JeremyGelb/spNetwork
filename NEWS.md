# spNetwork 0.4.3.7

## corrected bugs

Corrected an error in the way the border correction was calculated for the TNKDE

# spNetwork 0.4.3.6.9000

## corrected bugs

Corrected a minor bug in the functions network_listw.mc and network_listw

# spNetwork 0.4.3.6

## corrected bugs

The new version of sf created an error because of a bad use of the argument by_element in the function st_distance. It is corrected and tested now. 

# spNetwork 0.4.3.5

## New Features

* It is now possible to create isochrones as donught instead of plain isochrones (no overlapping lines). This can be handy if the are used for mapping and to limit the size of the produced feature collections.

## corrected bugs

A bug occuring when splitting the data with a spatial grid has been corrected.

## Documentation

* There is a new vignette on the website about spatio-temporal dbscan with network distance.

# spNetwork 0.4.3.4

## corrected bugs

* corrected two minor bugs with the kfunctions, kfunctions.mc, cross_kfunctions and cross_kfunctions.mc ([issue 8](https://github.com/JeremyGelb/spNetwork/issues/8))

# spNetwork 0.4.3.3

## corrected bugs

* When using the nkde and nkde.mc functions, one could get an error if the samples or the events had a X or a Y column. This has been corrected.

# spNetwork 0.4.3.2

## other changes

* The package searchTrees is not on CRAN anymore. It has been replaced in spNetwork by home made c++ classes using the BH library.

# spNetwork 0.4.3.1

## other changes

* example datasets are now stored in .rda file and can be loaded with the function `data`

# spNetwork 0.4.1.9000

## New Features

* It is possible now to choose between two strategies to deal with density at 0 when calculating LOO likelihood for several bandwidths. For more details, see the parameter *zero_strat* in the functions *bw_cvl_calc*, *bw_cvl_calc.mc*, *bw_cv_likelihood_calc*, *bw_cv_likelihood_calc.mc*, *bws_tnkde_cv_likelihood_calc* and *bws_tnkde_cv_likelihood_calc.mc*.

* It is possible now to set directly the size of local bandwidths for the functions *nkde*, *nkde.mc*, *tnkde* and *tnkde.mc*. This feature is currently experimental, please submit an issue if you encounter a bug.

* One can now use the identity matrix type for neighbour weighting (*network_listw*). It is useful if one wants to apply different function to calculate the weights and to compare the results.

* Adaptive bandwidth for TNKDE can now be calculated separately or simultaneously. For more details, see the parameter *adaptive_separate*.

## documentation

* The vignette *Details about NKDE"* has been removed from the package because it could not be knitted on some CRAN machines. It is still accessible on the website.

* A new vignette *K-nearest neighbour adaptive bandwidth* has also be added to the website. It presents how spNetwork can use user defined local bandwidths.

## corrected bugs

* The new c++ function for adding vertex on lines was buggy and could break the geometry of the original lines. This has been fixed.
* On Linux with gcc-12, the package couldn't be compiled. This has been fixed.

# spNetwork 0.4.0.9000

## New Features

* the package has been reworked to use sf rather than rgeos, maptools and sp. All the tests passed and the documentation has been modified.

# spNetwork 0.3.0.9000

## New Features

* Added a set of functions to perform Temporal Network Kernel Density Estimation (see vignette TNKDE)

## documentation

* A new vignette is presenting the Temporal Network Kernel Density Estimation
* Adding a CITATION file

## performance

* Performance was improved for data-driven bandwidth selection for NKDE
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

* Added two functions for data-driven bandwidth selection (see vignette NKDE)
* Added a function to simplify network before processing (*simplify_network*), reducing network complexity by removing useless edges or merge continuous edges can improve performance. This function is still experimental and must be used only for undirected networks.
* Added two functions to calculate k nearest neighbours: *network_knn* and *network_knn.mc*
* Added a function to calculate isochrones on a geographical network: *calc_isochrones* and a Vignette (*Calculating isochrones*) to give an example.
* Making some internal functions visible for the user to allow graph analysis on networks, see the new vignette *Building networks*.

## documentation

Several vignettes were added to present the new features. A pkgdown website was also built to give a more attractive way to browse the documentation.

We follow now the code coverage. We want to reach 70% coverage before the next CRAN release.

## corrected bugs

* issue 1 fixed by changing the function maptools::snapPointsToLines for a home-made function. The package BH including the c++ boost library is used to realize fast spatial query.
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

