#include "spNetwork.h"

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
