// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include<sstream>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <queue>

#define _USE_MATH_DEFINES
#include <math.h>

using namespace Rcpp;
using namespace std;
using namespace arma;
// [[Rcpp::depends(RcppProgress)]]
//[[Rcpp::plugins(cpp11)]]
typedef arma::vec (*fptr)(arma::vec, double);
typedef double (*fptros)(double, double);
