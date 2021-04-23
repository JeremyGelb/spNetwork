#include "spNetwork.h"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### base kernel functions ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arma::vec quartic_kernel(arma::vec d, double bw){
  arma::vec u = d/bw ;
  arma::vec k = ((15.0/16.0)*arma::pow((1-arma::pow(u,2)),2)) / bw;
  arma::uvec test = arma::find(d>=bw);
  k.elem(test).fill(0.0);
  return k;
}

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



arma::vec triangle_kernel(arma::vec d, double bw){
  arma::vec u = d/bw ;
  arma::vec k = (1.0 - arma::abs(u)) / bw;
  arma::uvec test = arma::find(d>=bw);
  k.elem(test).fill(0.0);
  return k;
}

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


arma::vec uniform_kernel(arma::vec d, double bw){
  arma::vec k = d;
  k.fill(1.0/(bw*2.0));
  arma::uvec test = arma::find(d>=bw);
  k.elem(test).fill(0.0);
  return k;
}

double uniform_kernelos(double d, double bw){
  double k;
  if (d>bw){
    k=0;
  }else{
    k = 1.0/(bw*2.0);
  };
  return k;
}

arma::vec epanechnikov_kernel(arma::vec d, double bw){
  arma::vec u = d/bw ;
  arma::vec k = ((3.0/4.0) * (1.0-arma::pow(u,2))) / bw;
  arma::uvec test = arma::find(d>=bw);
  k.elem(test).fill(0.0);
  return k;
}

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


arma::vec triweight_kernel(arma::vec d, double bw){
  arma::vec u = d/bw ;
  arma::vec k = ((35.0/32.0) * arma::pow((1.0-arma::pow(u,2)),3)) / bw ;
  arma::uvec test = arma::find(d>=bw);
  k.elem(test).fill(0.0);
  return k;
}

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


arma::vec tricube_kernel(arma::vec d, double bw){
  arma::vec u = d/bw ;
  arma::vec k = ((70.0/81.0) * arma::pow((1.0-pow(arma::abs(u),3)),3)) / bw;
  arma::uvec test = arma::find(d>=bw);
  k.elem(test).fill(0.0);
  return k;
}

double tricube_kernelos(double d, double bw){
  double u = d/bw ;
  double k;
  if(d > bw){
    k = 0;
  }else{
    double u2 = (1.0-u) * (1.0-u) * (1.0-u);
    k = ((70.0/81.0) * (u2*u2*u2)) / bw;
  };
  return k;
}


arma::vec cosine_kernel(arma::vec d, double bw){
  arma::vec u = d/bw ;
  arma::vec k = ((M_PI/4.0) * arma::cos((M_PI/2.0)*u)) / bw;
  arma::uvec test = arma::find(d>=bw);
  k.elem(test).fill(0.0);
  return k;
}

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

arma::vec gaussian_kernel(arma::vec d, double bw){
  arma::vec u = d/bw ;
  double t1 = 1.0/(sqrt(2.0*M_PI));
  arma::vec k = (t1 * (arma::exp(-1.0 * (1.0/2.0) * arma::pow(u,2)))) / bw;
  arma::uvec test = arma::find(d>=bw);
  k.elem(test).fill(0.0);
  return k;
}

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

arma::vec gaussian_kernel_scaled(arma::vec d, double bw){
  double bw2 = bw/3.0;
  arma::vec u = d/bw2 ;
  double t1 = 1.0/(sqrt(2.0*M_PI));
  arma::vec k = (t1 * (arma::exp(-1.0 * (1.0/2.0) * arma::pow(u,2)))) / bw2;
  arma::uvec test = arma::find(d>=bw);
  k.elem(test).fill(0.0);
  return k;
}

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
fptr select_kernel(std::string c) {
  if (c.compare("gaussian")==0 ) {
    return gaussian_kernel;
  }
  if (c.compare("scaled gaussian")==0 ) {
    return gaussian_kernel_scaled;
  }
  if (c.compare("cosine")==0 ) {
    return cosine_kernel;
  }
  if (c.compare("tricube")==0 ) {
    return tricube_kernel;
  }
  if (c.compare("triweight")==0 ) {
    return triweight_kernel;
  }
  if (c.compare("epanechnikov")==0 ) {
    return epanechnikov_kernel;
  }
  if (c.compare("triangle")==0 ) {
    return triangle_kernel;
  }
  if (c.compare("uniform")==0 ) {
    return uniform_kernel;
  }
  if (c.compare("quartic")==0 ) {
    return quartic_kernel;
  }
  return quartic_kernel;
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
}
