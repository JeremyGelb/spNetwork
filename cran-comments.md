# Version 0.1.0

### Test environments

* local R installation, R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel)

### R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Round 1 (after automatic checks)

**problem**:  
Possibly mis-spelled words in DESCRIPTION:  
  listw (10:23)  
  reticular (11:52)  
  spNetwork (8:18)  
  
**correction**:  
These words are not mis-spelled.

<br>
<br>
**problem**:  
Found the following (possibly) invalid URLs:  
URL: http://github.com/JeremyGelb/spNetwork (moved to https://github.com/JeremyGelb/spNetwork)  
URL: http://github.com/JeremyGelb/spNetwork/issues (moved to https://github.com/JeremyGelb/spNetwork/issues)  

**correction**:  
The *http* were replaced by *https*

<br>
<br>
**problem**:  
Package in Depends/Imports which should probably only be in LinkingTo: 'RcppArmadillo'

**correction**:  
RcppArmadillo is directly used in the package, this note is spurious

<br>
<br>
**problem**:  
checking whether package 'spNetwork' can be installed ... WARNING  
Found the following significant warnings:  
  code.cpp:601:18: warning: suggest parentheses around comparison in operand of '&' [-Wparentheses]  
  code.cpp:712:18: warning: suggest parentheses around comparison in operand of '&' [-Wparentheses]  
  code.cpp:974:20: warning: suggest parentheses around comparison in operand of '&' [-Wparentheses]  
  code.cpp:1083:20: warning: suggest parentheses around comparison in operand of '&' [-Wparentheses]  
  
**correction**:  
The parentheses were added as suggested

<br>
<br>
**problem**:  
checking data for non-ASCII characters ... WARNING  
AREA["Canada - Quebec and Ontario - 75<c2><b0>W to 72<c2><b0>W"]

**correction**: 
This warning is identified for each dataset provided with the package spNetwork. The error is caused by the way SpatialDataFrames from **sp** package store CRS information. I removed all the datasets from the data file and saved them in a classical GIS format (gpkg) in the folder inst/extdata. I edited all the examples, vignettes and tests to use load the data from this directory. However, this leads to a quite large directory (3.4 MB) and raise a new note.

<br>
<br>
**problem**:  
Examples with CPU (user + system) or elapsed time > 10s

**correction**: 
all the examples with long example time are not ran anymore during checks (\\dontrun{})

<br>
<br>
### R CMD check results

0 errors | 0 warnings | 3 note

checking CRAN incoming feasibility ... NOTE  
  Maintainer: 'Jeremy Gelb <jeremy.gelb@ucs.inrs.ca>'  
  
  New submission  
(this is a new submission, no problem here)

checking package dependencies ... NOTE  
  Package in Depends/Imports which should probably only be in LinkingTo: 'RcppArmadillo'  
  
(as stated above, if RcppArmadillo is only linked, errors occur)

checking installed package size ... NOTE  
    installed size is 11.1Mb  
    sub-directories of 1Mb or more:  
      extdata   3.6Mb  
      libs      6.4Mb  
      
(as stated above, libs directory is large because of the use of Rcpp. The directory extdata is large too because data used for examples raises warnings when included as .rda files. The datasets are SpatialDataFrames from sp package and they store their Coordinates Reference System with non-ASCII characters. They are thus stored in geopackage files. in extdata directory.)
