# spNetwork 0.1.2.9000

## New Features

* Added two function for datadriven bandwidth selectionn (see vignette NKDE)
* Added a function to simplify network before processing (*simplify_network*), reducing network complexity by removing useless edges or merge continuous edge can improve performance.

## corrected bugs

* issue 1 fixed by changing the function maptools::snapPointsToLines for a home made function. There is room for improvement here because it requires a conversion to sf and then back to sp which can be unsafe.
* issue 2 fixed by handling the special case of lines of length 0.

# spNetwork 0.1.1

* Accepted CRAN version (minor changes)

# spNetwork 0.1.0
  
* Initial release version
