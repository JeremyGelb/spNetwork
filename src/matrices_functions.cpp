#include "spNetwork.h"

// a simple function to create a vector of values between a start and an end with defined step
std::vector<double> seq_num(double start, double end, double step){

  std::vector<double> values;
  double cumul = 0;
  while(cumul+step < end){
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
