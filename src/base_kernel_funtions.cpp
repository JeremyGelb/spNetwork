#include "spNetwork.h"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### base kernel functions ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title c++ quartic kernel
//' @name quartic_kernel_cpp
//' @param d a vector of distances for which the density must be calculated
//' @param bw a double representing the size of the kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
arma::vec quartic_kernel_cpp(arma::vec d, double bw){
  arma::vec k = ((15.0/16.0)*arma::pow((1-arma::pow((d/bw),2)),2)) / bw;
  k.elem((arma::find(d>=bw))).fill(0.0);
  return k;
}

//' @title c++ quartic kernel for one distance
//' @name quartic_kernelos
//' @param d a double, the distances for which the density must be calculated
//' @param bw a double representing the size of the kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
double quartic_kernelos(double d, double bw){
  double u = d/bw ;
  double k;
  if(d > bw){
    k = 0;
  }else{
    double p2 = 1 - (u*u);
    k = ((15.0/16.0)*(p2*p2)) / bw;
  };
  return k;
}


//' @title c++ triangle kernel
//' @name triangle_kernel_cpp
//' @param d a vector of distances for which the density must be calculated
//' @param bw a double representing the size of the kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
arma::vec triangle_kernel_cpp(arma::vec d, double bw){
  arma::vec k = (1.0 - arma::abs((d/bw))) / bw;
  k.elem((arma::find(d>=bw))).fill(0.0);
  return k;
}

//' @title c++ triangle kernel for one distance
//' @name triangle_kernelos
//' @param d a double, the distances for which the density must be calculated
//' @param bw a double representing the size of the kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
double triangle_kernelos(double d, double bw){
  double u = d/bw ;
  double k;
  if(d>bw){
    k = 0;
  }else{
    k = (1.0 - u) / bw;
  };
  return k;
}

//' @title c++ uniform kernel
//' @name uniform_kernel_cpp
//' @param d a vector of distances for which the density must be calculated
//' @param bw a double representing the size of the kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
arma::vec uniform_kernel_cpp(arma::vec d, double bw){
  arma::vec k = d;
  k.fill(1.0/(bw*2.0));
  k.elem((arma::find(d>=bw))).fill(0.0);
  return k;
}

//' @title c++ uniform kernel for one distance
//' @name uniform_kernelos
//' @param d a double, the distances for which the density must be calculated
//' @param bw a double representing the size of the kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
double uniform_kernelos(double d, double bw){
  double k;
  if (d>bw){
    k=0;
  }else{
    k = 1.0/(bw*2.0);
  };
  return k;
}

//' @title c++ epanechnikov kernel
//' @name epanechnikov_kernel_cpp
//' @param d a vector of distances for which the density must be calculated
//' @param bw a double representing the size of the kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
arma::vec epanechnikov_kernel_cpp(arma::vec d, double bw){
  arma::vec k = ((3.0/4.0) * (1.0-arma::pow((d/bw),2))) / bw;
  k.elem((arma::find(d>=bw))).fill(0.0);
  return k;
}

//' @title c++ epanechnikov kernel for one distance
//' @name epanechnikov_kernelos
//' @param d a double, the distances for which the density must be calculated
//' @param bw a double representing the size of the kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
double epanechnikov_kernelos(double d, double bw){
  double u = d/bw ;
  double k;
  if (d > bw){
    k=0;
  }else{
    k = ((3.0/4.0) * (1.0-(u*u))) / bw;
  };
  return k;
}

//' @title c++ triweight kernel
//' @name triweight_kernel_cpp
//' @param d a vector of distances for which the density must be calculated
//' @param bw a double representing the size of the kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
arma::vec triweight_kernel_cpp(arma::vec d, double bw){
  arma::vec k = ((35.0/32.0) * arma::pow((1.0-arma::pow((d/bw),2)),3)) / bw ;
  k.elem((arma::find(d>=bw))).fill(0.0);
  return k;
}

//' @title c++ triweight kernel for one distance
//' @name triweight_kernelos
//' @param d a double, the distances for which the density must be calculated
//' @param bw a double representing the size of the kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
double triweight_kernelos(double d, double bw){
  double u = d/bw ;
  double k;
  if (d > bw){
    k=0;
  }else{
    double u2 = 1.0-(u*u);
    k = ((35.0/32.0) * (u2 * u2 * u2)) / bw;
  };
  return k;
}


//' @title c++ tricube kernel
//' @name tricube_kernel_cpp
//' @param d a vector of distances for which the density must be calculated
//' @param bw a double representing the size of the kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
arma::vec tricube_kernel_cpp(arma::vec d, double bw){
  arma::vec k = ((70.0/81.0) * arma::pow((1.0-pow(arma::abs((d/bw)),3)),3)) / bw;
  k.elem((arma::find(d>=bw))).fill(0.0);
  return k;
}

//' @title c++ tricube kernel for one distance
//' @name tricube_kernelos
//' @param d a double, the distances for which the density must be calculated
//' @param bw a double representing the size of the kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
double tricube_kernelos(double d, double bw){
  double u = d/bw ;
  double k;
  if(d > bw){
    k = 0;
  }else{
    double u3 = 1-u*u*u;
    k = ((70.0/81.0) * (u3*u3*u3)) / bw;
  };
  return k;
}

//' @title c++ cosine kernel
//' @name cosine_kernel_cpp
//' @param d a vector of distances for which the density must be calculated
//' @param bw a double representing the size of the kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
arma::vec cosine_kernel_cpp(arma::vec d, double bw){
  arma::vec k = ((M_PI/4.0) * arma::cos((M_PI/2.0)*(d/bw))) / bw;
  k.elem((arma::find(d>=bw))).fill(0.0);
  return k;
}

//' @title c++ cosine kernel for one distance
//' @name cosine_kernelos
//' @param d a double, the distances for which the density must be calculated
//' @param bw a double representing the size of the kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
double cosine_kernelos(double d, double bw){
  double u = d/bw ;
  double k;
  if (d > bw){
    k = 0;
  }else{
    double u2 = cos((M_PI/2.0)*u);
    k = ((M_PI/4.0) * u2) / bw;
  };
  return k;
}

//' @title c++ gaussian kernel
//' @name gaussian_kernel_cpp
//' @param d a vector of distances for which the density must be calculated
//' @param bw a double representing the size of the kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
arma::vec gaussian_kernel_cpp(arma::vec d, double bw){
  arma::vec u = d/bw ;
  double t1 = 1.0/(sqrt(2.0*M_PI));
  arma::vec k = (t1 * (arma::exp(-1.0 * (1.0/2.0) * arma::pow(u,2)))) / bw;
  arma::uvec test = arma::find(d>=bw);
  k.elem(test).fill(0.0);
  return k;
}

//' @title c++ gaussian kernel for one distance
//' @name gaussian_kernelos
//' @param d a double, the distances for which the density must be calculated
//' @param bw a double representing the size of the kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
double gaussian_kernelos(double d, double bw){
  double u2 = (d/bw) * (d/bw);
  double t1 = 1.0/(sqrt(2.0*M_PI));
  double k;
  if(d > bw){
    k=0;
  }else{
    k = (t1 * exp(-1.0 * 0.5 * u2)) / bw;
  };
  return k;
}

//' @title c++ scale gaussian kernel
//' @name gaussian_kernel_scaled_cpp
//' @param d a vector of distances for which the density must be calculated
//' @param bw a double representing the size of the kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
arma::vec gaussian_kernel_scaled_cpp(arma::vec d, double bw){
  double bw2 = bw/3.0;
  arma::vec u = d/bw2 ;
  double t1 = 1.0/(sqrt(2.0*M_PI));
  arma::vec k = (t1 * (arma::exp(-1.0 * (1.0/2.0) * arma::pow(u,2)))) / bw2;
  arma::uvec test = arma::find(d>=bw);
  k.elem(test).fill(0.0);
  return k;
}

//' @title c++ scaled gaussian kernel for one distance
//' @name gaussian_kernel_scaledos
//' @param d a double, the distances for which the density must be calculated
//' @param bw a double representing the size of the kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
double gaussian_kernel_scaledos(double d, double bw){
  double bw2 = bw/3.0;
  double u2 = (d/bw2) * (d/bw2) ;
  double t1 = 1.0/(sqrt(2.0*M_PI));
  double k;
  if(d > bw){
    k=0;
  }else{
    k = (t1 * exp(-1.0 * 0.5 * u2)) / bw2;
  };
  return k;
}


// a function to select the right kernel according to its name (vectorized version)
fptr select_kernel(std::string c) { //# nocov start
  if (c.compare("gaussian")==0 ) {
    return gaussian_kernel_cpp;
  }
  if (c.compare("scaled gaussian")==0 ) {
    return gaussian_kernel_scaled_cpp;
  }
  if (c.compare("cosine")==0 ) {
    return cosine_kernel_cpp;
  }
  if (c.compare("tricube")==0 ) {
    return tricube_kernel_cpp;
  }
  if (c.compare("triweight")==0 ) {
    return triweight_kernel_cpp;
  }
  if (c.compare("epanechnikov")==0 ) {
    return epanechnikov_kernel_cpp;
  }
  if (c.compare("triangle")==0 ) {
    return triangle_kernel_cpp;
  }
  if (c.compare("uniform")==0 ) {
    return uniform_kernel_cpp;
  }
  if (c.compare("quartic")==0 ) {
    return quartic_kernel_cpp;
  }
  return quartic_kernel_cpp;
}

// a function to select the right kernel according to its name (one shot version)
fptros select_kernelos(std::string c) {
  if (c.compare("gaussian")==0 ) {
    return gaussian_kernelos;
  }
  if (c.compare("scaled gaussian")==0 ) {
    return gaussian_kernel_scaledos;
  }
  if (c.compare("cosine")==0 ) {
    return cosine_kernelos;
  }
  if (c.compare("tricube")==0 ) {
    return tricube_kernelos;
  }
  if (c.compare("triweight")==0 ) {
    return triweight_kernelos;
  }
  if (c.compare("epanechnikov")==0 ) {
    return epanechnikov_kernelos;
  }
  if (c.compare("triangle")==0 ) {
    return triangle_kernelos;
  }
  if (c.compare("uniform")==0 ) {
    return uniform_kernelos;
  }
  if (c.compare("quartic")==0 ) {
    return quartic_kernelos;
  }
  return quartic_kernelos;
} //# nocov end


//#####################################################################################
// #####################  some common functions  ####################################
//#####################################################################################


// A simple function to find duplicates in an integer vector
std::map<int, int> count_values_intvec(IntegerVector x){
  // Creating a map
  std::map<int, int> counter;
  // Filling the map
  for (int i = 0; i != x.length(); ++i) {
    counter[ x[i] ] += 1;
  }
  return counter;
}

// A simple function to check if a map contains a value
bool map_contains_key(std::map<int, arma::vec>& map, int key){
  if (map.find(key) == map.end()) return false;
  return true;
}

