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
  // Iterating on breaks
  for(int i = 0; i < breaks.size(); ++i) {
    float dist = breaks[i];
    arma::mat int_mat = arma::conv_to<arma::mat>::from((dist_mat <= dist));
    int_mat.each_col() %= w;
    int_mat.diag().zeros();
    k_values(i) = arma::accu(int_mat) * t1;
  }
  return k_values;
}


//' @title c++ k function counting worker
//' @name kfunc_counting
//' @description c++ k function counting (INTERNAL)
//' @param dist_mat A matrix with the distances between points
//' @param wc The weight of the points represented by the columns (destinations)
//' @param wr The weight of the points represented by the rows (origins)
//' @param breaks A numeric vector with the distance to consider
//' @param cross A boolean indicating if we are calculating a cross k function or not (default is FALSE)
//' @return A numeric matrix with the countings of the k function evaluated at the required distances
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericMatrix kfunc_counting(arma::mat dist_mat, arma::rowvec wc, NumericVector wr, NumericVector breaks, bool cross = false){


  float end = max(breaks);

  NumericMatrix counting_k(dist_mat.n_rows, breaks.size());

  int malus = 1;
  if(cross){
    malus = 0;
  }

  // we iterate over each row of dist_mat
  for(int i = 0; i < dist_mat.n_rows; ++i) {

    float wi = wr(i);
    arma::rowvec row = dist_mat.row(i);

    // we do a first test to remove all the distances that are way to big
    arma::uvec test = arma::find(row <= end);

    arma::vec ok_d = row.elem(test) ;
    arma::vec ok_w = wc.elem(test);


    for(int z = 0; z < breaks.size(); ++z) {

      if(ok_d.n_elem == 0){
        break;
      }

      float dist = breaks[z];

      // for each iteration we do a fist cleaning test
      arma::uvec test3 = ok_d <= (dist);

      ok_d = ok_d.elem(arma::find(test3)) ;
      ok_w = ok_w.elem(arma::find(test3)) ;

      // minus one here is important to remove self weight
      counting_k(i,z) = (arma::sum(ok_w)-malus) * wi;

    }

  }

  return counting_k;
}


//' @title c++ k function 2
//' @name kfunc_cpp2
//' @description c++ k function (INTERNAL)
//' @param dist_mat A square matrix with the distances between points
//' @param start A float, the start value for evaluating the k-function
//' @param end A float, the last value for evaluating the k-function
//' @param step A float, the jump between two evaluations of the k-function
//' @param Lt The total length of the network
//' @param n The number of points
//' @param wc The weight of the points represented by the columns (destinations)
//' @param wr The weight of the points represented by the rows (origins)
//' @param cross A boolean indicating if we are calculating a cross k function or not (default is FALSE)
//' @return A numeric vector with the values of the k function evaluated at the required distances
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericVector kfunc_cpp2(arma::mat dist_mat,float start,float end, float step, float Lt, int n, arma::rowvec wc, NumericVector wr, bool cross = false){

  std::vector<float> breaks0 = seq_num3(start,end,step);
  std::reverse(breaks0.begin(), breaks0.end());
  NumericVector breaks = wrap(breaks0);
  float t1;
  if(cross){
    t1 = 1.0/((n)/Lt);
  }else{
    t1 = 1.0/((n-1.0)/Lt);
  }



  NumericMatrix counting = kfunc_counting(dist_mat, wc,wr, breaks, cross);

  float div = sum(wr);

  NumericVector k_values = rcppRev((colSums(counting) / div) * t1) ;

  return k_values;
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### base g function ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//' @title c++ g function counting worker
//' @name gfunc_counting
//' @description c++ k function counting (INTERNAL)
//' @param dist_mat A matrix with the distances between points
//' @param wc The weight of the points represented by the columns (destinations)
//' @param wr The weight of the points represented by the rows (origins)
//' @param breaks A numeric vector with the distance to consider
//' @param width The width of each donut
//' @return A numeric matrix with the countings of the g function evaluated at the required distances
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericMatrix gfunc_counting(arma::mat dist_mat, arma::colvec wc, NumericVector wr, NumericVector breaks, float width){


  NumericMatrix counting(dist_mat.n_rows, breaks.size());
  float end = max(breaks);

  // we iterate over each row of dist_mat
  for(int i = 0; i < dist_mat.n_rows; ++i) {

    arma::rowvec row = dist_mat.row(i);

    float wi = wr(i);

    // we do a first test to remove all the distances that are way to big
    arma::uvec test = arma::find(row <= end+width);

    arma::vec ok_d = row.elem(test) ;
    arma::vec ok_w = wc.elem(test);
    arma::vec ok_d2 = row.elem(test) ;
    arma::vec ok_w2 = wc.elem(test);

    for(int z = 0; z < breaks.size(); ++z) {

      float dist = breaks[z];

      // for each iterationm we do a the test in both direction separately
      arma::uvec test1 = ok_d <= (dist+width);
      arma::uvec test2 = ok_d >= (dist-width);

      // we can the use them to do the test conjointly
      arma::uvec trueTest = arma::find( (test1) && (test2) );


      // with the full test we can get the values
      ok_w2 = ok_w.elem(trueTest) ;
      // minus one here is important to remove self weight
      counting(i,z) = (arma::sum(ok_w2)-1) * wi;

      // with the first part of the test, we can reduce the next research
      ok_d = ok_d.elem(arma::find(test1)) ;
      ok_w = ok_w.elem(arma::find(test1)) ;
    }


  }

  return counting;
}


//' @title c++ g function
//' @name gfunc_cpp2
//' @description c++ g function (INTERNAL)
//' @param dist_mat A square matrix with the distances between points
//' @param start A float, the start value for evaluating the g-function
//' @param end A float, the last value for evaluating the g-function
//' @param step A float, the jump between two evaluations of the k-function
//' @param width The width of each donut
//' @param Lt The total length of the network
//' @param n The number of points
//' @param wc The weight of the points represented by the columns (destinations)
//' @param wr The weight of the points represented by the rows (origins)
//' @return A numeric vector with the values of the g function evaluated at the required distances
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericVector gfunc_cpp2(arma::mat dist_mat,float start,float end, float step, float width, float Lt, int n, arma::rowvec wc, NumericVector wr){


  std::vector<float> breaks0 = seq_num3(start,end,step);
  std::reverse(breaks0.begin(), breaks0.end());
  NumericVector breaks = wrap(breaks0);
  float t1 = 1.0/((n-1)/Lt);
  width = width/2.0;

  NumericMatrix counting = gfunc_counting(dist_mat, wc, wr, breaks, width);


  float div = sum(wr);

  NumericVector cppBreaks  = wrap(breaks) ;
  NumericVector k_values = rcppRev((colSums(counting) / div) * t1) ;

  return k_values;
}



//' @title c++ k and g function counting worker
//' @name kgfunc_counting
//' @description c++ k function counting (INTERNAL)
//' @param dist_mat A matrix with the distances between points
//' @param wc The weight of the points represented by the columns (destinations)
//' @param wr The weight of the points represented by the rows (origins)
//' @param breaks A numeric vector with the distance to consider
//' @param width The width of each donut
//' @param cross A boolean indicating if we are calculating a cross k function or not (default is FALSE)
//' @return A list  of two numeric matrices with the values of the k and g function evaluated at the required distances
//' @keywords internal
//' @export
// [[Rcpp::export]]
List kgfunc_counting(arma::mat dist_mat, arma::rowvec wc, NumericVector wr, NumericVector breaks, float width, double cross = false){

  // NumericMatrix counting(dist_mat.n_rows, breaks.size());
  float end = max(breaks);

  NumericMatrix counting_k(dist_mat.n_rows, breaks.size());
  NumericMatrix counting_g(dist_mat.n_rows, breaks.size());

  int cross_malus = 1;
  if(cross){
    cross_malus = 0;
  }

  // we iterate over each row of dist_mat
  for(int i = 0; i < dist_mat.n_rows; ++i) {

    float wi = wr(i);
    arma::rowvec row = dist_mat.row(i);

    // we do a first test to remove all the distances that are way to big
    arma::uvec test = arma::find(row <= end+width);

    arma::vec ok_d = row.elem(test) ;
    arma::vec ok_w = wc.elem(test);
    arma::vec ok_d2 = row.elem(test) ;
    arma::vec ok_w2 = wc.elem(test);
    arma::vec ok_w1 = wc.elem(test);


    for(int z = 0; z < breaks.size(); ++z) {

      if(ok_d.n_elem == 0){
        break;
      }

      float dist = breaks[z];
      // for each iteration we do a the test in both direction separately
      arma::uvec test1 = ok_d <= (dist+width);
      arma::uvec test2 = ok_d >= (dist-width);
      arma::uvec test3 = ok_d <= (dist);

      // for the g func, we must check that dist - width is > 0, otherwise we must apply
      // a malus to not do self counting
      int malus = 0;
      if(dist-width <= 0){
        malus = 1;
      }

      // we can the use them to do the test conjointly
      arma::uvec trueTest = arma::find((test1) && (test2));

      // with the part test I can get the value of the k function
      ok_w1 = ok_w.elem(arma::find(test3)) ;
      // minus one here is important to remove self weight
      counting_k(i,z) = (arma::sum(ok_w1)-cross_malus) * wi;


      // with the full test we can get the values for the g function
      ok_w2 = ok_w.elem(trueTest) ;
      // for g we do not include minus one because the donut form will prevent self count
      counting_g(i,z) = (arma::sum(ok_w2) - (malus * cross_malus)) * wi;

      // with the first part of the test, we can reduce the next research
      ok_d = ok_d.elem(arma::find(test1)) ;
      ok_w = ok_w.elem(arma::find(test1)) ;
    }

  }

  List results = List::create(counting_k, counting_g);

  return results;
}


//' @title c++ k and g function
//' @name kgfunc_cpp2
//' @description c++ g function (INTERNAL)
//' @param dist_mat A square matrix with the distances between points
//' @param start A float, the start value for evaluating the g-function
//' @param end A float, the last value for evaluating the g-function
//' @param step A float, the jump between two evaluations of the k-function
//' @param width The width of each donut
//' @param Lt The total length of the network
//' @param n The number of points
//' @param wc The weight of the points represented by the columns (destinations)
//' @param wr The weight of the points represented by the rows (origins)
//' @param cross A boolean indicating if we are calculating a cross k function or not (default is FALSE)
//' @return A numeric matrix with the values of the k (first col) and g (second col) function evaluated at the required distances
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericMatrix kgfunc_cpp2(arma::mat dist_mat,float start,float end, float step, float width, float Lt, int n, arma::rowvec wc, NumericVector wr, bool cross = false){


  std::vector<float> breaks0 = seq_num3(start,end,step);
  std::reverse(breaks0.begin(), breaks0.end());
  NumericVector breaks = wrap(breaks0);
  float t1;
  if(cross){
    t1 = 1.0/((n)/Lt);
  }else{
    t1 = 1.0/((n-1)/Lt);
  }
  width = width/2.0;

  // we start here by counting for each distance how many points are
  // reachable from each point
  // the produced matrice have a number of column equal to the number
  // of breaks, and a number of row equal to the number of point

  List elements = kgfunc_counting(dist_mat, wc, wr, breaks, width, cross);

  NumericMatrix counting_k = elements[0];
  NumericMatrix counting_g = elements[1];


  // we must apply here the weights to the rows when calculating the means
  float div = sum(wr);

  NumericVector k_values = rcppRev((colSums(counting_k) / div) * t1) ;
  NumericVector g_values = rcppRev((colSums(counting_g) / div) * t1) ;
  NumericMatrix final = cbind(k_values, g_values);

  return final;
}


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
  float t1 = 1.0/((na-1)/Lt);

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
  float t1 = 1.0/((na-1)/Lt);
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



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### base g space-time function ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title c++ g space-time function
//' @name g_nt_func_cpp
//' @param dist_mat_net A square matrix with the distances between points on the network
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
