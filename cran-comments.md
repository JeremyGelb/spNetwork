# Version 0.1.0

### Test environments

* local R installation, R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel)

### R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Round 1 (after automatic checks) the 01/12/2020

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
      
(as stated above, libs directory is large because of the use of Rcpp. The directory extdata is large too because data used for examples raises warnings when included as .rda files. The datasets are SpatialDataFrames from sp package and they store their Coordinates Reference System with non-ASCII characters. They are thus stored in geopackage files, in extdata directory.)


## Round 2 (after first comments) the 07/12/2020

**Problem**:  
please omit the redudnant "The spNetwork package provides functions to perform ". Perhaps simply replace by "Performs"
  
**correction**: 
The DESCRIPTION has been modified as suggested:  
The spNetwork package performs spatial analysis on network.
  
**Problem**:  
Check: Overall checktime, Result: NOTE
   Overall checktime 28 min > 10 min

mainly from

* checking re-building of vignette outputs ... [21m] OK

Can this be reduced to few minuted to keep the overall check time < 10
min? For example, use few itarations, toy data and if not possible by
reducing the comlexity of the computations, please provide precomputed
results for the lengthy parts.  

**correction**:
The duration of vignettes building has been reduced by precomputing the results for
the two longest vignettes.  

The new durations are :  
KNetworkFunctions.Rmd : 18 seconds  
SpatialWeightMatrices.Rmd : 79 seconds  
NKDE.Rmd : 69 seconds  
NKDEdetailed.Rmd : 8 seconds

## Round 3 (after second comments) the 31/12/2020

" Size of tarball: 37777188 bytes

   Overall checktime 16 min > 10 min"

The size of the tarball has been reduced to 3.5 Mo (one data file was not compressed, my bad).
The Overall checktime is now only  4min on my laptop.

## Round 4 (after third comments, human validation) the 08/01/2021

**Problem**:  
Please omit the redundant "The spNetwork package" from your description.

**correction**:
The description has been modified as requested

**Problem**:
Please always write package names, software names and API (application
programming interface) names in single quotes in title and description.
e.g: --> 'spNetwork'

**correction**:
As requested, packages names and API are now written with single quotes in the description and the vignettes. The spNetwork name is still in bold in the vignettes.

**Problem**:
if there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <[https:...]https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

**correction**:
The main references were added in the description as suggested.

**Problem**:
Warning: Unexecutable code in man/cross_kfunctions.mc.Rd:
Warning: Unexecutable code in man/kfunctions.mc.Rd:
Warning: Unexecutable code in man/lixelize_lines.mc.Rd:
Warning: Unexecutable code in man/network_listw.mc.Rd:
Warning: Unexecutable code in man/nkde.mc.Rd:

**correction**:
The problem is caused by the use of dontshow for the multicore functions. I removed it. I think it is a good choice to show to the user how to come back to the sequential plan after multiprocessing.

**Problem**:
Please always make sure to reset to user's options(), working directory
or par() after you changed it in examples and vignettes and demos.
e.g.: inst/doc/NKDEdetailed.R, inst/doc/KNetworkFunctions.R
oldpar <- par(mfrow = c(1,2))
...
par(oldpar)

**correction**:
As suggested, the modifications applied to the user's option, par and working directory are now reset after vignettes. These parameters are never modified in examples.

# Version 0.1.1

## Round 1 (after automatic checks) the 18/01/2021

The package was accepted on CRAN today, but an error occurred on SOLARIS.
The error occurs when I set manually the CRS of sp objects in the test and in the vignette NKDEdetailed. I removed the lines specifying the projection (maybe an error with PROJ4 ?). Because the package has already been accepted on CRAN, I have to raise the version number.


# Version 0.2.0

submitted the 08/10/2021

This is the result of RCMD check before submission.

-- R CMD check results ------------------------------------------------------ spNetwork 0.2.0 ----
Duration: 8m 21.4s

> checking installed package size ... NOTE
    installed size is  9.5Mb
    sub-directories of 1Mb or more:
      doc       2.0Mb
      extdata   4.0Mb
      libs      2.8Mb

0 errors √ | 0 warnings √ | 1 note x
