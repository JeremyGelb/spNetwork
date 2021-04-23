#ifndef base_kernel_functions
#define base_kernel_functions

arma::vec quartic_kernel(arma::vec d, double bw);
double quartic_kernelos(double d, double bw);


arma::vec triangle_kernel(arma::vec d, double bw);
double triangle_kernelos(double d, double bw);

arma::vec uniform_kernel(arma::vec d, double bw);
double uniform_kernelos(double d, double bw);

arma::vec epanechnikov_kernel(arma::vec d, double bw);
double epanechnikov_kernelos(double d, double bw);


arma::vec triweight_kernel(arma::vec d, double bw);
double triweight_kernelos(double d, double bw);

arma::vec tricube_kernel(arma::vec d, double bw);
double tricube_kernelos(double d, double bw);

arma::vec cosine_kernel(arma::vec d, double bw);
double cosine_kernelos(double d, double bw);

arma::vec gaussian_kernel(arma::vec d, double bw);
double gaussian_kernelos(double d, double bw);

arma::vec gaussian_kernel_scaled(arma::vec d, double bw);
double gaussian_kernel_scaledos(double d, double bw);


// a function to select the right kernel according to its name (vectorized version)
fptr select_kernel(std::string c);

// a function to select the right kernel according to its name (one shot version)
fptros select_kernelos(std::string c);

#endif
