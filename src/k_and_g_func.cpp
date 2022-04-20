#include "spNetwork.h"
#include "matrices_functions.h"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### base k function ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title c++ k function
//' @name kfunc_cpp
//' @description c++ k function (INTERNAL)
//' @param dist_mat A square matrix with the distances between points
//' @param start A float, the start value for evaluating the k-function
//' @param end A float, the last value for evaluating the k-function
//' @param step A float, the jump between two evaluations of the k-function
//' @param Lt The total length of the network
//' @param n The number of points
//' @param w The weight of the points (coincident points)
//' @return A numeric vector with the values of the k function evaluated at the required distances
//' @export
// [[Rcpp::export]]
NumericVector kfunc_cpp(arma::mat dist_mat,float start,float end, float step, float Lt, int n, arma::colvec w){

  std::vector<double> breaks = seq_num2(start,end,step);
  float t1 = (n-1)/Lt;

  NumericVector k_values(breaks.size());

  for(int i = 0; i < breaks.size(); ++i) {
    float dist = breaks[i];
    arma::mat int_mat = arma::conv_to<arma::mat>::from((dist_mat <= dist));
    int_mat.each_col() %= w;
    int_mat.diag().zeros();
    k_values(i) = arma::accu(int_mat) * t1;
  }
  return k_values;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### base g function ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title c++ g function
//' @name gfunc_cpp
//' @description c++ g function (INTERNAL)
//' @param dist_mat A square matrix with the distances between points
//' @param start A float, the start value for evaluating the g-function
//' @param end A float, the last value for evaluating the g-function
//' @param step A float, the jump between two evaluations of the k-function
//' @param width The width of each donut
//' @param Lt The total length of the network
//' @param n The number of points
//' @param w The weight of the points (coincident points)
//' @return A numeric vector with the values of the g function evaluated at the required distances
//' @export
// [[Rcpp::export]]
NumericVector gfunc_cpp(arma::mat dist_mat,float start,float end, float step, float width, float Lt, int n, arma::colvec w){

  std::vector<double> breaks = seq_num2(start,end,step);
  float t1 = (n-1)/Lt;
  width = width/2.0;
  NumericVector k_values(breaks.size());

  for(int i = 0; i < breaks.size(); ++i) {
    float dist = breaks[i];
    arma::mat int_mat = arma::conv_to<arma::mat>::from(( (dist_mat<= (dist+width)) && (dist_mat>=(dist-width)) ));
    int_mat.each_col() %= w;
    int_mat.diag().zeros();
    k_values(i) = arma::accu(int_mat) * t1;
  }

  return k_values;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### base cross k function ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title c++ cross k function
//' @name cross_kfunc_cpp
//' @param dist_mat A square matrix with the distances between points
//' @param start A float, the start value for evaluating the k-function
//' @param end A float, the last value for evaluating the k-function
//' @param step A float, the jump between two evaluations of the k-function
//' @param Lt The total length of the network
//' @param na The number of points in set A
//' @param nb The number of points in set B
//' @param wa The weight of the points in set A (coincident points)
//' @param wb The weight of the points in set B (coincident points)
// [[Rcpp::export]]
NumericVector cross_kfunc_cpp(arma::mat dist_mat,float start,float end, float step, float Lt, int na, int nb, arma::rowvec wa, arma::colvec wb){

  // note : in the matrix, the rows are the b points
  // and the columns are the a points
  std::vector<double> breaks = seq_num2(start,end,step);
  float t1 = (na)/Lt;

  NumericVector k_values(breaks.size());

  for(int i = 0; i < breaks.size(); ++i) {
    float dist = breaks[i];
    // applying the a weight (row wise)
    arma::mat int_mat = arma::conv_to<arma::mat>::from((dist_mat <= dist));
    int_mat.each_row() %= wa;
    int_mat.each_col() %= wb;
    k_values(i) = arma::accu(int_mat) * t1;
  }
  return k_values;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### base cross g function ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title c++ cross g function
//' @name cross_gfunc_cpp
//' @description c++ cross g function (INTERNAL)
//' @param dist_mat A matrix with the distances between points
//' @param start A float, the start value for evaluating the g-function
//' @param end A float, the last value for evaluating the g-function
//' @param step A float, the jump between two evaluations of the k-function
//' @param width The width of each donut
//' @param Lt The total length of the network
//' @param na The number of points in set A
//' @param nb The number of points in set B
//' @param wa The weight of the points in set A (coincident points)
//' @param wb The weight of the points in set B (coincident points)
// [[Rcpp::export]]
NumericVector cross_gfunc_cpp(arma::mat dist_mat,float start,float end, float step, float width, float Lt, int na, int nb, arma::rowvec wa, arma::colvec wb){

  // note : in the matrix, the rows are the b points
  // and the columns are the a points
  std::vector<double> breaks = seq_num2(start,end,step);
  float t1 = (na)/Lt;
  width = width/2.0;

  NumericVector k_values(breaks.size());

  for(int i = 0; i < breaks.size(); ++i) {
    float dist = breaks[i];
    // applying the a weight (row wise)
    arma::mat int_mat = arma::conv_to<arma::mat>::from( (dist_mat <= (dist+ width)) && (dist_mat >= (dist - width)) );
    int_mat.each_row() %= wa;
    int_mat.each_col() %= wb;
    k_values(i) = arma::accu(int_mat) * t1;
  }
  return k_values;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### base k space-time function ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title c++ k space-time function
//' @name k_nt_func_cpp
//' @param dist_mat_net A square matrix with the distances between points (network)
//' @param dist_mat_time A square matrix with the distances between points (time)
//' @param start_net A float, the start value for evaluating the k-function (network)
//' @param end_net A float, the last value for evaluating the k-function (network)
//' @param step_net A float, the jump between two evaluations of the k-function (network)
//' @param start_time A float, the start value for evaluating the k-function (time)
//' @param end_time A float, the last value for evaluating the k-function (time)
//' @param step_time A float, the jump between two evaluations of the k-function (time)
//' @param Lt The total length of the network
//' @param Tt The total duration of study area
//' @param n The number of points
//' @param w The weight of the points (coincident points)
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix k_nt_func_cpp(arma::mat dist_mat_net, arma::mat dist_mat_time,
                             float start_net,float end_net, float step_net,
                             float start_time,float end_time, float step_time,
                             float Lt, float Tt, int n, arma::colvec w){

  std::vector<double> breaks_net = seq_num2(start_net,end_net,step_net);
  std::vector<double> breaks_time = seq_num2(start_time,end_time,step_time);
  float t1 = (n-1)/(Lt * Tt);

  NumericMatrix k_values(breaks_net.size(),breaks_time.size());
  arma::mat int_mat(dist_mat_net.n_rows, dist_mat_net.n_cols);
  int_mat.zeros();
  arma::umat mat1(dist_mat_net.n_rows, dist_mat_net.n_cols);
  mat1.zeros();

  // pre-calculating the conditions in time
  std::vector<arma::umat> umat_times;
  for(int i = 0; i < breaks_time.size(); ++i){
    umat_times.push_back(dist_mat_time <= breaks_time[i]);
  }

  for(int i = 0; i < breaks_net.size(); ++i) {
    //mat1 = dist_mat_net <= breaks_net[i];
    mat1.elem(arma::find(dist_mat_net <= breaks_net[i])).ones();
    for(int j = 0; j < breaks_time.size(); ++j){
      int_mat.elem( arma::find((mat1) && (umat_times[j]))).ones();
      int_mat.each_col() %= w;
      int_mat.diag().zeros();
      k_values(i,j) = arma::accu(int_mat) * t1;
      int_mat.zeros();
    }
    mat1.zeros();
  }
  return k_values;
}


// [[Rcpp::export]]
IntegerMatrix k_nt_func_cpp2(arma::imat dist_mat_net, arma::imat dist_mat_time,
                            int start_net, int end_net, int step_net,
                            int start_time,int end_time, int step_time,
                            int Lt, int Tt, int n, arma::icolvec w){

  std::vector<int> breaks_net = seq_num2f(start_net,end_net,step_net);
  std::vector<int> breaks_time = seq_num2f(start_time,end_time,step_time);
  int t1 = (n-1)/(Lt * Tt);

  IntegerMatrix k_values(breaks_net.size(),breaks_time.size());
  arma::imat int_mat(dist_mat_net.n_rows, dist_mat_net.n_cols);
  int_mat.zeros();
  arma::umat mat1(dist_mat_net.n_rows, dist_mat_net.n_cols);
  mat1.zeros();

  // pre-calculating the conditions in time
  std::vector<arma::umat> umat_times;
  for(int i = 0; i < breaks_time.size(); ++i){
    umat_times.push_back(dist_mat_time <= breaks_time[i]);
  }

  for(int i = 0; i < breaks_net.size(); ++i) {
    //mat1 = dist_mat_net <= breaks_net[i];
    mat1.elem(arma::find(dist_mat_net <= breaks_net[i])).ones();
    for(int j = 0; j < breaks_time.size(); ++j){
      int_mat.elem( arma::find((mat1) && (umat_times[j]))).ones();
      int_mat.each_col() %= w;
      int_mat.diag().zeros();
      k_values(i,j) = arma::accu(int_mat) * t1;
      int_mat.zeros();
    }
    mat1.zeros();
  }
  return k_values;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### base g space-time function ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title c++ g space-time function
//' @name g_nt_func_cpp
//' @param dist_mat A square matrix with the distances between points on the network
//' @param dist_mat_time A square matrix with the distances between points in time
//' @param start_net A float, the start value for evaluating the g-function on the network
//' @param end_net A float, the last value for evaluating the g-function on the network
//' @param step_net A float, the jump between two evaluations of the g-function on the network
//' @param width_net The width of each donut on the network
//' @param start_time A float, the start value for evaluating the g-function in time
//' @param end_time A float, the last value for evaluating the g-function in time
//' @param step_time A float, the jump between two evaluations of the g-function in time
//' @param width_time The width of each donut in time
//' @param Lt The total length of the network
//' @param n The number of points
//' @param w The weight of the points (coincident points)
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix g_nt_func_cpp(arma::mat dist_mat_net, arma::mat dist_mat_time,
                            float start_net,float end_net, float step_net, float width_net,
                            float start_time,float end_time, float step_time, float width_time,
                            float Lt, float Tt, float n, arma::colvec w){

  std::vector<double> breaks_net = seq_num2(start_net,end_net,step_net);
  std::vector<double> breaks_time = seq_num2(start_time,end_time,step_time);
  float t1 = (n-1.0)/(Lt * Tt);

  NumericMatrix k_values(breaks_net.size(),breaks_time.size());
  arma::mat int_mat(dist_mat_net.n_rows, dist_mat_net.n_cols);
  int_mat.zeros();
  arma::umat mat1(dist_mat_net.n_rows, dist_mat_net.n_cols);
  mat1.zeros();

  // pre-calculating the conditions in time
  std::vector<arma::umat> umat_times;
  for(int i = 0; i < breaks_time.size(); ++i){
    umat_times.push_back(
      ( (dist_mat_time<= (breaks_time[i]+width_time)) && (dist_mat_time>=(breaks_time[i]-width_time)) )
    );
  }

  for(int i = 0; i < breaks_net.size(); ++i) {
    mat1.elem(arma::find(
        ( (dist_mat_net<= (breaks_net[i]+width_net)) && (dist_mat_net>=(breaks_net[i]-width_net)) )
    )).ones();
    for(int j = 0; j < breaks_time.size(); ++j){
      float dtime = breaks_time[j];
      int_mat.elem(arma::find(((mat1) && (umat_times[j])))).ones();
      int_mat.each_col() %= w;
      int_mat.diag().zeros();
      k_values(i,j) = arma::accu(int_mat) * t1;
      int_mat.zeros();
    }
    mat1.zeros();
  }
  return k_values;
}
