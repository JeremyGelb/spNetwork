#include "spNetwork.h"

// a simple function to create a vector of values between a start and an end with defined step
// [[Rcpp::export]]
std::vector<double> seq_num2(double start, double end, double step){

  std::vector<double> values;
  double cumul = start - step;
  while(cumul+step <= end){
    cumul+=step;
    values.push_back(cumul);
  }
  return values;
}

// a simple function to create a vector of values between a start and an end with defined step
// [[Rcpp::export]]
std::vector<int> seq_num2f(int start, int end, int step){

  std::vector<int> values;
  int cumul = start - step;
  while(cumul+step <= end){
    cumul+=step;
    values.push_back(cumul);
  }
  return values;
}


// a simple function to find the index of the first occurence of value in a numeric vector
int get_first_index(NumericVector& v1, double x){
  int i;
  for( i = 0; i < v1.size(); ++i) {
    if(v1[i] == x){
      return i;
    }
  }
  return -1;
}

// a simple function to find the index of the first occurence of value in a numeric vector
int get_first_index_int(IntegerVector& v1, int x){
  int i;
  for( i = 0; i < v1.size(); ++i) {
    if(v1[i] == x){
      return i;
    }
  }
  return -1;
}

// a simple function to find the index of the first occurence of value in a numeric vector
std::vector<int> get_all_indeces(NumericVector& v1, double x){
  int i;
  std::vector<int> idxs;
  for( i = 0; i < v1.size(); ++i) {
    if(v1[i] == x){
      idxs.push_back(i);
    }
  }
  return idxs;
}

// a simple function to find the index of the first occurence of value in an integer vector
std::vector<int> get_all_indeces_int(IntegerVector& v1, int x){
  int i;
  std::vector<int> idxs;
  for( i = 0; i < v1.size(); ++i) {
    if(v1[i] == x){
      idxs.push_back(i);
    }
  }
  return idxs;
}



// short function to create a matrix from the neighbour_list
// note : might be better as a sparse matrix ?
IntegerMatrix make_matrix(DataFrame df, List neighbour_list){
  List edge_list;
  IntegerVector starts = df["start_oid"];
  IntegerVector ends = df["end_oid"];
  IntegerVector edge_id = df["graph_id"];

  //making all the elements empty
  int cnt1 = neighbour_list.length();
  IntegerMatrix edge_mat(cnt1+1,cnt1+1);

  //then filling it !
  int cnt = starts.length();
  for(int i=0; i < cnt; ++i){
    edge_mat(starts[i],ends[i]) = edge_id[i];
    edge_mat(ends[i],starts[i]) = edge_id[i];
  }
  return edge_mat;
}



arma::sp_mat make_matrix_sparse(DataFrame df, List neighbour_list){
  List edge_list;
  IntegerVector starts = df["start_oid"];
  IntegerVector ends = df["end_oid"];
  IntegerVector edge_id = df["graph_id"];

  //making all the elements empty
  int cnt1 = neighbour_list.length();
  arma::sp_mat edge_mat(cnt1+1,cnt1+1);

  //then filling it !
  int cnt = starts.length();
  for(int i=0; i < cnt; ++i){
    edge_mat(starts[i],ends[i]) = edge_id[i];
    edge_mat(ends[i],starts[i]) = edge_id[i];
  }
  return edge_mat;
}



arma::sp_imat make_imatrix_sparse(DataFrame df, List neighbour_list){
  List edge_list;
  IntegerVector starts = df["start_oid"];
  IntegerVector ends = df["end_oid"];
  IntegerVector edge_id = df["graph_id"];

  //making all the elements empty
  int cnt1 = neighbour_list.length();
  arma::sp_imat edge_mat(cnt1+1,cnt1+1);

  //then filling it !
  int cnt = starts.length();
  for(int i=0; i < cnt; ++i){
    edge_mat(starts[i],ends[i]) = edge_id[i];
    edge_mat(ends[i],starts[i]) = edge_id[i];
  }
  return edge_mat;
}


arma::sp_mat make_edge_weight_sparse(DataFrame& df, List& neighbour_list){
  List edge_list;
  IntegerVector starts = df["start_oid"];
  IntegerVector ends = df["end_oid"];
  NumericVector edge_weight = df["weight"];

  //making all the elements empty
  int cnt1 = neighbour_list.length();
  arma::sp_mat edge_mat(cnt1+1,cnt1+1);

  //then filling it !
  int cnt = starts.length();
  for(int i=0; i < cnt; ++i){
    edge_mat(starts[i],ends[i]) = edge_weight[i];
    edge_mat(ends[i],starts[i]) = edge_weight[i];
  }
  return edge_mat;
}


// a little function used in space-time kfunctions to extend a matrix by duplicating rows and cols
// for events that have been aggregated
// oids is a unique id for each original event
// locids is the location id of each original event (its location in the matrix)
// [[Rcpp::export]]
NumericMatrix extend_matrix_by_ids(NumericMatrix agg_mat, IntegerVector oids, IntegerVector locids){

  // creating the matrix that will receive the values
  NumericMatrix new_mat(oids.size(), oids.size());

  // filling it by row and col
  for(int i = 0; i < oids.size(); ++i){
    int locid_i = locids(i);
    for(int j = 0; j < oids.size(); ++j){
      int locid_j = locids(j);
      new_mat(i,j) = agg_mat(locid_i, locid_j);
    }
  }
  return new_mat;
}


// a little function used to reverse the order of the row in a matrix
// borrowed here : https://stackoverflow.com/questions/31946641/reverse-numeric-matrix-by-row
// TY Dirk !
// [[Rcpp::export]]
NumericMatrix reverseByRow(NumericMatrix inmat) {
  int r = inmat.nrow();
  NumericMatrix nw(r,inmat.ncol());
  for(int i = 0; i < r; i++){
    nw.row(i) = inmat.row(r-i-1);
  }
  return nw;
}


//' @title euclidean distance between rows of a matrix and a vector (arma mode)
//' @name calcEuclideanDistance3
//' @param y a matrix
//' @param x a vector (same length as ncol(matrix))
//' @return a vector (same length as nrow(matrix))
//' @export
//' @keywords internal
//'
// [[Rcpp::export]]
arma::colvec calcEuclideanDistance3(arma::mat y, arma::mat x){
   return arma::sqrt(arma::sum(arma::pow(y.each_row() - x,2),1));
 }
