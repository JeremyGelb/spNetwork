# spNetwork 0.1.2.9000

## New Features

* Added two function for datadriven bandwidth selectionn (see vignette NKDE)
* Added a function to simplify network before processing (*simplify_network*), reducing network complexity by removing useless edges or merge continuous edge can improve performance.
* Added two functions to calculate k nearest neighbours: *network_knn* and *network_knn.mc*
* Added a function to calculate isochrones on a geographical network: *calc_isochrones* and a Vignette (*Calculating isochrones*) to give an example.
* Making some internal functions visible for the user to allow graph analysis on networks, see the new vignette *Building networks*.

## documentation

Several vignettes were added to present the new features. A pkgdown website was aslo built to give a more attractive way to browse the documentation.

## corrected bugs

* issue 1 fixed by changing the function maptools::snapPointsToLines for a home made function. The package BH including the c++ boost library is used to realize fast spatial query.
* issue 2 fixed by handling the special case of lines of length 0.
* an error was raised when using *listw_network* with polygons and the method = "pointsalong". This was caused by a minor error in the code and works now as expected.
* an error was present in the Diggle correction factor, leading to overestimation of the correction. This error has been corrected.

## performance

* Performance was improved for function *listw_network*, *listw_network.mc*, *kfunctions*, *kfunctions.mc*, *cross_kfunctions* and *cross_kfunctions.mc*, by reducing the complexity of the network built from the SpatialLinesDataFrame with no loss of accuracy.
* Several geometrical functions have been accelerated by c++ code using the boost library
* upgraded c++ compilation requirement to c++17

# spNetwork 0.1.1

* Accepted CRAN version (minor changes)

# spNetwork 0.1.0
  
* Initial release version
