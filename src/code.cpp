// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include<sstream>
#include <progress.hpp>
#include <progress_bar.hpp>
//#include <unordered_map>
#include <queue>

#define _USE_MATH_DEFINES
#include <math.h>

using namespace Rcpp;
using namespace std;
using namespace arma;
// [[Rcpp::depends(RcppProgress)]]
//[[Rcpp::plugins(cpp11)]]

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### General helper functions ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

typedef arma::vec (*fptr)(arma::vec, double);
typedef double (*fptros)(double, double);

void rcpp_rcout_num(NumericVector v){
  // printing value of vector
  Rcout << "" << v << "\n";
}

void rcpp_rcout_int(IntegerVector v){
  // printing value of vector
  Rcout << "The value of v : " << v << "\n";
}


// //this dict work as expected !
// unordered_map<int,unordered_map<int,int>> make_map(DataFrame df){
//   unordered_map<int,unordered_map<int,int>> mymap;
//   IntegerVector starts = df["start_oid"];
//   IntegerVector ends = df["end_oid"];
//   IntegerVector edge_id = df["graph_id"];
//   int cnt = starts.length();
//   for(int i=0; i < cnt; ++i){
//     mymap[starts[i]][ends[i]]=edge_id[i];
//     mymap[ends[i]][starts[i]]=edge_id[i];
//   }
//   return mymap;
// }

// short function to create a matrix from the neighbour_list
// note : might be better as a sparse matrix ?
IntegerMatrix make_matrix(DataFrame df, List neighbour_list){
  List edge_list;
  IntegerVector starts = df["start_oid"];
  IntegerVector ends = df["end_oid"];
  IntegerVector edge_id = df["graph_id"];

  //making all the elements empty
  int cnt1 = neighbour_list.length();
  IntegerMatrix edge_mat(cnt1,cnt1);

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


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### NKDE continuous functions ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//' @title The worker function to calculate continuous NKDE (with ARMADILLO and sparse matrix)
//' @name continuousWorker_sparse
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param samples_k a numeric vector of the actual kernel values, updates at
//' each recursion
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param v the actual node to consider for the recursion (int)
//' @param v1 the connected node to consider for the recursion (int)
//' @param l1 the edge connecting v and v1 (int)
//' @param d the actual distance traveled before the recursion
//' @param alpha the actual alpha value before the recursion
//' @param bw the kernel bandwidth
//' @param line_weights a vector with the length of the edges
//' @param samples_edgeid a vector associating each sample to an edge
//' @param samples_x a vector with x coordinates of each sample
//' @param samples_y a vector with y coordinates of each sample
//' @param nodes_x a vector with x coordinates of each node
//' @param nodes_y a vector with y coordinates of each node
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a vector with the kernel values calculated for each samples from
//' the first node given
arma::vec esc_kernel_rcpp_arma_sparse(fptr kernel_func, arma::vec samples_k, List neighbour_list, arma::sp_mat edge_mat, int v, int v1, int l1, double d,double alpha, double bw, NumericVector line_weights, arma::vec samples_edgeid, arma::vec samples_x, arma::vec samples_y, arma::vec nodes_x, arma::vec nodes_y , int depth, int max_depth){

  //step1 find the index of the right samples
  arma::uvec test = arma::find(samples_edgeid==l1);
  arma::vec sampling_x = samples_x.elem(test);
  arma::vec sampling_y = samples_y.elem(test);
  //extracting the X and Y coordinates of the starting node
  int v_m1 = v-1;
  double node_x = nodes_x[v_m1];
  double node_y = nodes_y[v_m1];
  //step2 calculating the distances for each sample
  arma::vec x_dists = arma::sqrt(arma::pow((sampling_x - node_x),2) + arma::pow((sampling_y - node_y),2)) + d;
  //step3 calculating the values of the new kernel
  arma::vec new_k = kernel_func(x_dists,bw)*alpha;
  //note:here we seems to always have 0 values
  samples_k.elem(test) += new_k;
  //mettre a jour d
  double d2 = line_weights[l1-1] + d;
  if((bw>d2) && (depth < max_depth)){
    //on veut trouver toutes les lignes emannant de v (Lv)
    IntegerVector v_neighbours = neighbour_list[v1-1];
    int n = v_neighbours.length();
    int new_depth;
    //updating depth
    if(n>2){
      new_depth = depth+1;
    }else{
      new_depth = depth+0;
    }
    if(n>1){
      for(int j=0; j < n; ++j){
        //int li = Lv[j];
        int vi = v_neighbours[j];
        //int li = edge_dict[v1][vi];
        int li = edge_mat(v1,vi);
        if(li==l1){
          //we must backfire only if we have a true intersection
          if(n>2){
            double p2 = (2.0-n)/n;
            double n_alpha = alpha * p2;
            samples_k = esc_kernel_rcpp_arma_sparse(kernel_func, samples_k, neighbour_list, edge_mat, v1,vi, li, d2, n_alpha, bw, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, new_depth, max_depth);
          }
        }else{
          double n_alpha = alpha * (2.0/n);
          samples_k = esc_kernel_rcpp_arma_sparse(kernel_func, samples_k, neighbour_list, edge_mat, v1,vi, li, d2, n_alpha, bw, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, new_depth, max_depth);
        };
      };
    };
  };
  return samples_k;
}


//' @title The worker function to calculate continuous NKDE (with ARMADILLO and Integer matrix)
//' @name continuousWorker_int
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param samples_k a numeric vector of the actual kernel values, updates at
//' each recursion
//' @param neighbour_list a List, providing an IntegerVector for each node with
//' its neighbours
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param v the actual node to consider for the recursion (int)
//' @param v1 the connected node to consider for the recursion (int)
//' @param l1 the edge connecting v and v1 (int)
//' @param d the actual distance traveled before the recursion
//' @param alpha the actual alpha value before the recursion
//' @param bw the kernel bandwidth
//' @param line_weights a vector with the length of the edges
//' @param samples_edgeid a vector associating each sample to an edge
//' @param samples_x a vector with x coordinates of each sample
//' @param samples_y a vector with y coordinates of each sample
//' @param nodes_x a vector with x coordinates of each node
//' @param nodes_y a vector with y coordinates of each node
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a vector with the kernel values calculated for each samples from
//' the first node given
arma::vec esc_kernel_rcpp_arma(fptr kernel_func, arma::vec samples_k, List neighbour_list, IntegerMatrix edge_mat, int v, int v1, int l1, double d,double alpha, double bw, NumericVector line_weights, arma::vec samples_edgeid, arma::vec samples_x, arma::vec samples_y, arma::vec nodes_x, arma::vec nodes_y , int depth, int max_depth){
  //step1 find the index of the right samples
  arma::uvec test = arma::find(samples_edgeid==l1);
  arma::vec sampling_x = samples_x.elem(test);
  arma::vec sampling_y = samples_y.elem(test);
  //extracting the X and Y coordinates of the starting node
  int v_m1 = v-1;
  double node_x = nodes_x[v_m1];
  double node_y = nodes_y[v_m1];
  //step2 calculating the distances for each sample
  arma::vec x_dists = arma::sqrt(arma::pow((sampling_x - node_x),2) + arma::pow((sampling_y - node_y),2)) + d;
  //step3 calculating the values of the new kernel
  arma::vec new_k = kernel_func(x_dists,bw)*alpha;
  //note:here we seems to always have 0 values
  samples_k.elem(test) += new_k;
  //mettre a jour d
  double d2 = line_weights[l1-1] + d;
  if((bw>d2) && (depth < max_depth)){
    //on veut trouver toutes les lignes emannant de v (Lv)
    IntegerVector v_neighbours = neighbour_list[v1-1];
    int n = v_neighbours.length();
    int new_depth;
    //updating depth
    if(n>2){
      new_depth = depth+1;
    }else{
      new_depth = depth+0;
    }
    if(n>1){
      for(int j=0; j < n; ++j){
        //int li = Lv[j];
        int vi = v_neighbours[j];
        //int li = edge_dict[v1][vi];
        int li = edge_mat(v1,vi);
        if(li==l1){
          //we must backfire only if we have a true intersection
          if(n>2){
            double p2 = (2.0-n)/n;
            double n_alpha = alpha * p2;
            samples_k = esc_kernel_rcpp_arma(kernel_func, samples_k, neighbour_list, edge_mat, v1,vi, li, d2, n_alpha, bw, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, new_depth, max_depth);
          }
        }else{
          double n_alpha = alpha * (2.0/n);
          samples_k = esc_kernel_rcpp_arma(kernel_func, samples_k, neighbour_list, edge_mat, v1,vi, li, d2, n_alpha, bw, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, new_depth, max_depth);
        };
      };
    };
  };
  return samples_k;
}




//' @title The main function to calculate continuous NKDE (with ARMADILO and sparse matrix)
//' @name continuousfunction
//' @param neighbour_list a list of the neighbours of each node
//' @param events a numeric vector of the node id of each event
//' @param weights a numeric vector of the weight of each event
//' @param samples a DataFrame of the samples (with spatial coordinates and belonging edge)
//' @param bws the kernel bandwidths for each event
//' @param kernel_name the name of the kernel to use
//' @param nodes a DataFrame representing the nodes of the graph (with spatial coordinates)
//' @param line_list a DataFrame representing the lines of the graph
//' @param max_depth the maximum recursion depth (after which recursion is stopped)
//' @param verbose a boolean indicating if the function must print its progress
//' @return a DataFrame with two columns : the kernel values (sum_k) and the number of events for each sample (n)
//' @export
//'
// [[Rcpp::export]]
DataFrame continuous_nkde_cpp_arma_sparse(List neighbour_list, NumericVector events, NumericVector weights, DataFrame samples, NumericVector bws, std::string kernel_name, DataFrame nodes, DataFrame line_list, int max_depth, bool verbose){

  //continuous_nkde_cpp_arma_sparse(neighbour_list,events$vertex_id, events$weight, samples@data, bws, kernel_name, nodes@data, graph_result$linelist, max_depth, tol, verbose)
  //selecting the kernel function
  fptr kernel_func = select_kernel(kernel_name);

  //step0 extract the columns of the dataframe
  NumericVector line_weights = line_list["weight"];
  arma::vec samples_edgeid = as<arma::vec>(samples["edge_id"]);
  arma::vec samples_x =  as<arma::vec>(samples["X_coords"]);
  arma::vec samples_y =  as<arma::vec>(samples["Y_coords"]);
  arma::vec nodes_x = as<arma::vec>(nodes["X_coords"]);
  arma::vec nodes_y = as<arma::vec>(nodes["Y_coords"]);

  //step 1 : mettre toutes les valeurs a 0
  arma::vec base_k(samples.nrow());
  arma::vec base_count(samples.nrow());
  arma::vec samples_k(samples.nrow());
  arma::vec count(samples.nrow());
  base_k.fill(0.0);
  base_count.fill(0.0);
  count.fill(0.0);
  samples_k.fill(0.0);

  //calculer le dictionnaire des lignes
  //IntegerMatrix edge_mat = make_matrix(line_list,neighbour_list);
  arma::sp_mat edge_mat = make_matrix_sparse(line_list,neighbour_list);
  //step2 : iterer sur chaque event
  int cnt_e = events.length()-1;
  Progress p(cnt_e, verbose);
  for(int i=0; i <= cnt_e; ++i){
    p.increment(); // update progress
    //preparer les differentes valeurs de departs pour l'event y
    int y = events[i];
    double w = weights[i];
    double bw = bws[i];
    //on veut trouver toutes les voisins emannant de y
    IntegerVector y_neighbours = neighbour_list[y-1];
    int cnt_y = y_neighbours.length();
    double alpha = 2.0/cnt_y;
    for(int j=0; j < cnt_y; ++j){
      //preparing all values for this loop
      int vi = y_neighbours[j];
      //int li = edge_dict[y][vi];
      int li = edge_mat(y,vi);
      double d = 0.0 ;
      int depth = 0 ;
      // launching recursion
      samples_k.fill(0.0);
      samples_k = esc_kernel_rcpp_arma_sparse(kernel_func, samples_k, neighbour_list, edge_mat ,y,vi,li,d,alpha,bw, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, depth,max_depth);
      base_k += samples_k*w;
      // calculating the new value
      count.fill(0.0);
      arma::uvec ids = arma::find(samples_k>0);
      count.elem(ids).fill(1.0);
      base_count += count*w;
    };

  };
  NumericVector v1 = Rcpp::NumericVector(base_k.begin(), base_k.end());
  NumericVector v2 = Rcpp::NumericVector(base_count.begin(), base_count.end());
  DataFrame df =  DataFrame::create( Named("sum_k") = v1 ,Named("n") = v2);
  return df;
}


//' @title The main function to calculate continuous NKDE (with ARMADILO and Integer matrix)
//' @name continuousfunction
//' @param neighbour_list a list of the neighbours of each node
//' @param events a numeric vector of the node id of each event
//' @param weights a numeric vector of the weight of each event
//' @param samples a DataFrame of the samples (with spatial coordinates and belonging edge)
//' @param bws the kernel bandwidths for each event
//' @param kernel_name the name of the kernel to use
//' @param nodes a DataFrame representing the nodes of the graph (with spatial coordinates)
//' @param line_list a DataFrame representing the lines of the graph
//' @param max_depth the maximum recursion depth (after which recursion is stopped)
//' @param verbose a boolean indicating if the function must print its progress
//' @return a DataFrame with two columns : the kernel values (sum_k) and the number of events for each sample (n)
//' @export
//'
// [[Rcpp::export]]
DataFrame continuous_nkde_cpp_arma(List neighbour_list, NumericVector events, NumericVector weights, DataFrame samples, NumericVector bws, std::string kernel_name, DataFrame nodes, DataFrame line_list, int max_depth, bool verbose){

  //selecting the kernel function
  fptr kernel_func = select_kernel(kernel_name);

  //step0 extract the columns of the dataframe
  NumericVector line_weights = line_list["weight"];
  arma::vec samples_edgeid = as<arma::vec>(samples["edge_id"]);
  arma::vec samples_x =  as<arma::vec>(samples["X_coords"]);
  arma::vec samples_y =  as<arma::vec>(samples["Y_coords"]);
  arma::vec nodes_x = as<arma::vec>(nodes["X_coords"]);
  arma::vec nodes_y = as<arma::vec>(nodes["Y_coords"]);

  //step 1 : mettre toutes les valeurs a 0
  arma::vec base_k(samples.nrow());
  arma::vec base_count(samples.nrow());
  arma::vec samples_k(samples.nrow());
  arma::vec count(samples.nrow());
  base_k.fill(0.0);
  base_count.fill(0.0);
  count.fill(0.0);
  samples_k.fill(0.0);

  //calculer le dictionnaire des lignes
  IntegerMatrix edge_mat = make_matrix(line_list,neighbour_list);
  //step2 : iterer sur chaque event
  int cnt_e = events.length()-1;
  Progress p(cnt_e, verbose);
  for(int i=0; i <= cnt_e; ++i){
    p.increment(); // update progress
    //preparer les differentes valeurs de departs pour l'event y
    int y = events[i];
    double w = weights[i];
    double bw = bws[i];
    //on veut trouver toutes les voisins emannant de y
    IntegerVector y_neighbours = neighbour_list[y-1];
    int cnt_y = y_neighbours.length();
    double alpha = 2.0/cnt_y;
    for(int j=0; j < cnt_y; ++j){
      //preparing all values for this loop
      int vi = y_neighbours[j];
      //int li = edge_dict[y][vi];
      int li = edge_mat(y,vi);
      double d = 0.0 ;
      int depth = 0 ;
      // launching recursion
      samples_k.fill(0.0);
      samples_k = esc_kernel_rcpp_arma(kernel_func, samples_k, neighbour_list, edge_mat ,y,vi,li,d,alpha,bw, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, depth,max_depth);
      base_k += samples_k*w;
      // calculating the new value
      count.fill(0.0);
      arma::uvec ids = arma::find(samples_k>0);
      count.elem(ids).fill(1.0);
      base_count += count*w;
    };

  };
  NumericVector v1 = Rcpp::NumericVector(base_k.begin(), base_k.end());
  NumericVector v2 = Rcpp::NumericVector(base_count.begin(), base_count.end());
  DataFrame df =  DataFrame::create( Named("sum_k") = v1 ,Named("n") = v2);
  return df;
}




//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### NKDE discontinuous functions ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//' @title The worker function to calculate discontinuous NKDE (with ARMADILLO and sparse matrix)
//' @name discontinuousWorker_sparse
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param edge_mat matrix, to find the id of each edge given two neighbours
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param v the actual node to consider for the recursion (int)
//' @param bw the kernel bandiwdth
//' @param line_weights a vector with the length of the edges
//' @param samples_edgeid a vector associating each sample to an edge
//' @param samples_x a vector with x coordinates of each sample
//' @param samples_ya vector with y coordinates of each sample
//' @param nodes_x a vector with x coordinates of each node
//' @param nodes_y a vector with y coordinates of each node
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a vector with the kernel values calculated for each samples from
//' the first node given
arma::vec esd_kernel_rcpp_arma_sparse(fptr kernel_func, arma::sp_mat edge_mat,
                                      List neighbour_list ,int v, double bw,
                                      arma::vec line_weights, arma::vec samples_edgeid,
                                      arma::vec samples_x, arma::vec samples_y,
                                      arma::vec nodes_x, arma::vec nodes_y, int depth, int max_depth){
  //step0 : generate the queue
  queue <List> data_holder;
  arma::vec samples_k(samples_x.n_elem);
  samples_k.fill(0.0);
  //step1 : generate the first case
  List cas1 = List::create(Named("d")=0.0,
                           Named("alpha")=1.0,
                           Named("v") = v,
                           Named("prev_node") = -999,
                           Named("depth") = 0
  );
  data_holder.push(cas1);
  //lancement des iterations
  while(data_holder.empty()==FALSE){
    //unpacking (imagine some loop unrolling here with a function to deal with.)
    List cas = data_holder.front();
    data_holder.pop();
    int v = cas["v"];
    int depth = cas["depth"];
    double d = cas["d"];
    double alpha = cas["alpha"];
    int prev_node = cas["prev_node"];

    //step1 : find all the neighbours
    IntegerVector neighbours = neighbour_list[v-1];

    //step2 : iterate over the neighbours
    int cnt_n = neighbours.length();
    int new_depth;
    if(cnt_n>2){
      new_depth = depth+1;
    }else{
      new_depth = depth;
    }

    double new_alpha;
    if((prev_node < 0)  && (cnt_n > 2)){
      new_alpha = 2.0/(cnt_n);
    }else{
      new_alpha = alpha * (1.0/(cnt_n-1.0));
    }

    //if we have only one neighbour, we must stop
    if(cnt_n>1 or prev_node<=0){
      for(int i=0; i < cnt_n; ++i){
        int v2 = neighbours[i];
        //on ne veut pas revenir en arriere !
        if(v2!=prev_node){
          //find the edge between the two nodes
          int edge_id = edge_mat(v,v2);
          //find the samples on that edge
          arma::uvec test = arma::find(samples_edgeid==edge_id);
          arma::vec sampling_x = samples_x.elem(test);
          arma::vec sampling_y = samples_y.elem(test);

          //extracting the X and Y coordinates of the starting node
          int v_m1 = v-1;
          double node_x = nodes_x[v_m1];
          double node_y = nodes_y[v_m1];

          //calculating the distances
          arma::vec x_dists = arma::sqrt(arma::pow((sampling_x - node_x),2) + arma::pow((sampling_y - node_y),2)) + d;

          //step3 calculating the values of the new kernel
          arma::vec new_k = kernel_func(x_dists,bw)*new_alpha;
          samples_k.elem(test) += new_k;

          //evaluating for the next move
          double d2 = line_weights[edge_id-1] + d;

          if (d2<bw and new_depth<max_depth){
            List new_cas = List::create(Named("d")=d2,
                                        Named("alpha")=new_alpha,
                                        Named("v") = v2,
                                        Named("prev_node") = v,
                                        Named("depth") = new_depth
            );
            data_holder.push(new_cas);
          }
        }
      }
    }
  }

  return samples_k;

}

//' @title The worker function to calculate discontinuous NKDE (with ARMADILLO and Integer matrix)
//' @name discontinuousWorker_int
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param edge_mat matrix, to find the id of each edge given two neighbours
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param v the actual node to consider for the recursion (int)
//' @param bw the kernel bandiwdth
//' @param line_weights a vector with the length of the edges
//' @param samples_edgeid a vector associating each sample to an edge
//' @param samples_x a vector with x coordinates of each sample
//' @param samples_ya vector with y coordinates of each sample
//' @param nodes_x a vector with x coordinates of each node
//' @param nodes_y a vector with y coordinates of each node
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a vector with the kernel values calculated for each samples from
//' the first node given
arma::vec esd_kernel_rcpp_arma(fptr kernel_func, IntegerMatrix edge_mat,
                               List neighbour_list ,int v, double bw,
                               arma::vec line_weights, arma::vec samples_edgeid,
                               arma::vec samples_x, arma::vec samples_y,
                               arma::vec nodes_x, arma::vec nodes_y, int depth, int max_depth){
  //step0 : generate the queue
  queue <List> data_holder;
  arma::vec samples_k(samples_x.n_elem);
  samples_k.fill(0.0);
  //step1 : generate the first case
  List cas1 = List::create(Named("d")=0.0,
                           Named("alpha")=1.0,
                           Named("v") = v,
                           Named("prev_node") = -999,
                           Named("depth") = 0
  );
  data_holder.push(cas1);
  //lancement des iterations
  while(data_holder.empty()==FALSE){
    //unpacking (imagine some loop unrolling here with a function to deal with.)
    List cas = data_holder.front();
    data_holder.pop();
    int v = cas["v"];
    int depth = cas["depth"];
    double d = cas["d"];
    double alpha = cas["alpha"];
    int prev_node = cas["prev_node"];

    //step1 : find all the neighbours
    IntegerVector neighbours = neighbour_list[v-1];

    //step2 : iterate over the neighbours
    int cnt_n = neighbours.length();
    int new_depth;
    if(cnt_n>2){
      new_depth = depth+1;
    }else{
      new_depth = depth;
    }

    double new_alpha;
    if((prev_node < 0)  && (cnt_n > 2)){
      new_alpha = 2.0/(cnt_n);
    }else{
      new_alpha = alpha * (1.0/(cnt_n-1.0));
    }

    //if we have only one neighbour, we must stop
    if(cnt_n>1 or prev_node<=0){
      for(int i=0; i < cnt_n; ++i){
        int v2 = neighbours[i];
        //on ne veut pas revenir en arriere !
        if(v2!=prev_node){
          //find the edge between the two nodes
          int edge_id = edge_mat(v,v2);
          //find the samples on that edge
          arma::uvec test = arma::find(samples_edgeid==edge_id);
          arma::vec sampling_x = samples_x.elem(test);
          arma::vec sampling_y = samples_y.elem(test);

          //extracting the X and Y coordinates of the starting node
          int v_m1 = v-1;
          double node_x = nodes_x[v_m1];
          double node_y = nodes_y[v_m1];

          //calculating the distances
          arma::vec x_dists = arma::sqrt(arma::pow((sampling_x - node_x),2) + arma::pow((sampling_y - node_y),2)) + d;

          //step3 calculating the values of the new kernel
          arma::vec new_k = kernel_func(x_dists,bw)*new_alpha;
          samples_k.elem(test) += new_k;

          //evaluating for the next move
          double d2 = line_weights[edge_id-1] + d;

          if (d2<bw and new_depth<max_depth){
            List new_cas = List::create(Named("d")=d2,
                                        Named("alpha")=new_alpha,
                                        Named("v") = v2,
                                        Named("prev_node") = v,
                                        Named("depth") = new_depth
            );
            data_holder.push(new_cas);
          }
        }
      }
    }
  }

  return samples_k;

}




//' @title The main function to calculate discontinuous NKDE (ARMA and sparse matrix)
//' @name discontinuousfunction
//' @param neighbour_list a list of the neighbours of each node
//' @param events a numeric vector of the node id of each event
//' @param weights a numeric vector of the weight of each event
//' @param samples a DataFrame of the samples (with spatial coordinates and belonging edge)
//' @param bws the kernel bandwidths for each event
//' @param kernel_name the name of the kernel function to use
//' @param nodes a DataFrame representing the nodes of the graph (with spatial coordinates)
//' @param line_list a DataFrame representing the lines of the graph
//' @param max_depth the maximum recursion depth (after which recursion is stopped)
//' @param verbose a boolean indicating if the function must print its progress
//' @return a DataFrame with two columns : the kernel values (sum_k) and the number of events for each sample (n)
//' @export
//'
// [[Rcpp::export]]
DataFrame discontinuous_nkde_cpp_arma_sparse(List neighbour_list, NumericVector events, NumericVector weights, DataFrame samples, NumericVector bws, std::string kernel_name, DataFrame nodes, DataFrame line_list, int max_depth, bool verbose){

  //selecting the kernel function
  fptr kernel_func = select_kernel(kernel_name);

  //step0 extract the columns of the dataframe
  arma::vec line_weights =  as<arma::vec>(line_list["weight"]);
  arma::vec samples_edgeid = as<arma::vec>(samples["edge_id"]);
  arma::vec samples_x = as<arma::vec>(samples["X_coords"]);
  arma::vec samples_y = as<arma::vec>(samples["Y_coords"]);
  arma::vec nodes_x = as<arma::vec>(nodes["X_coords"]);
  arma::vec nodes_y = as<arma::vec>(nodes["Y_coords"]);
  //step 1 : mettre toutes les valeurs a 0
  arma::vec base_k(samples_x.n_elem);
  arma::vec base_count(samples_x.n_elem);
  arma::vec count(samples_x.n_elem);
  base_k.fill(0.0);
  base_count.fill(0.0);
  count.fill(0.0);


  //calculer la matrice des lignes
  //IntegerMatrix edge_mat = make_matrix(line_list,neighbour_list);
  arma::sp_mat edge_mat = make_matrix_sparse(line_list,neighbour_list);
  int depth = 0;

  //step2 : iterer sur chaque event
  int cnt_e = events.length()-1;
  Progress p(cnt_e, verbose);
  for(int i=0; i <= cnt_e; ++i){
    p.increment(); // update progress
    //preparer les differentes valeurs de departs pour l'event y
    int y = events[i];
    double w = weights[i];
    double bw = bws[i];
    // launching recursion
    //NumericVector k = esd_kernel_rcpp(samples_k, edge_dict, neighbour_list ,y,prev_node,d,alpha,bw,kernel_func, line_weights, samples_edgeid, samples_x, samples_y, samples_oid, nodes_x, nodes_y, depth,max_depth);
    arma::vec k = esd_kernel_rcpp_arma_sparse(kernel_func,edge_mat, neighbour_list ,y,bw, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, depth,max_depth);
    // getting the actual base_k values (summed at each iteration)
    base_k += (k*w);
    // calculating the new value
    count.fill(0.0);
    arma::uvec ids = arma::find(k>0);
    count.elem(ids).fill(1.0);
    //NumericVector new_base_count = base_count + (count*w);
    base_count += count*w;
  };
  NumericVector v1 = Rcpp::NumericVector(base_k.begin(), base_k.end());
  NumericVector v2 = Rcpp::NumericVector(base_count.begin(), base_count.end());
  DataFrame df =  DataFrame::create( Named("sum_k") = v1 ,Named("n") = v2);
  return df;
}


//' @title The main function to calculate discontinuous NKDE (ARMA and Integer matrix)
//' @name discontinuousfunction
//' @param neighbour_list a list of the neighbours of each node
//' @param events a numeric vector of the node id of each event
//' @param weights a numeric vector of the weight of each event
//' @param samples a DataFrame of the samples (with spatial coordinates and belonging edge)
//' @param bws the kernel bandwidth for each event
//' @param kernel_name the name of the kernel function to use
//' @param nodes a DataFrame representing the nodes of the graph (with spatial coordinates)
//' @param line_list a DataFrame representing the lines of the graph
//' @param max_depth the maximum recursion depth (after which recursion is stopped)
//' @param verbose a boolean indicating if the function must print its progress
//' @return a DataFrame with two columns : the kernel values (sum_k) and the number of events for each sample (n)
//' @export
//'
// [[Rcpp::export]]
DataFrame discontinuous_nkde_cpp_arma(List neighbour_list, NumericVector events, NumericVector weights, DataFrame samples, NumericVector bws, std::string kernel_name, DataFrame nodes, DataFrame line_list, int max_depth, bool verbose){

  //selecting the kernel function
  fptr kernel_func = select_kernel(kernel_name);

  //step0 extract the columns of the dataframe
  arma::vec line_weights =  as<arma::vec>(line_list["weight"]);
  arma::vec samples_edgeid = as<arma::vec>(samples["edge_id"]);
  arma::vec samples_x = as<arma::vec>(samples["X_coords"]);
  arma::vec samples_y = as<arma::vec>(samples["Y_coords"]);
  arma::vec nodes_x = as<arma::vec>(nodes["X_coords"]);
  arma::vec nodes_y = as<arma::vec>(nodes["Y_coords"]);
  //step 1 : mettre toutes les valeurs a 0
  arma::vec base_k(samples_x.n_elem);
  arma::vec base_count(samples_x.n_elem);
  arma::vec count(samples_x.n_elem);
  base_k.fill(0.0);
  base_count.fill(0.0);
  count.fill(0.0);


  //calculer la matrice des lignes
  IntegerMatrix edge_mat = make_matrix(line_list,neighbour_list);
  int depth = 0;

  //step2 : iterer sur chaque event
  int cnt_e = events.length()-1;
  Progress p(cnt_e, verbose);
  for(int i=0; i <= cnt_e; ++i){
    p.increment(); // update progress
    //preparer les differentes valeurs de departs pour l'event y
    int y = events[i];
    double w = weights[i];
    double bw = bws[i];
    // launching recursion
    //NumericVector k = esd_kernel_rcpp(samples_k, edge_dict, neighbour_list ,y,prev_node,d,alpha,bw,kernel_func, line_weights, samples_edgeid, samples_x, samples_y, samples_oid, nodes_x, nodes_y, depth,max_depth);
    arma::vec k = esd_kernel_rcpp_arma(kernel_func,edge_mat, neighbour_list ,y,bw, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, depth,max_depth);
    // getting the actual base_k values (summed at each iteration)
    base_k += (k*w);
    // calculating the new value
    count.fill(0.0);
    arma::uvec ids = arma::find(k>0);
    count.elem(ids).fill(1.0);
    //NumericVector new_base_count = base_count + (count*w);
    base_count += count*w;
  };
  NumericVector v1 = Rcpp::NumericVector(base_k.begin(), base_k.end());
  NumericVector v2 = Rcpp::NumericVector(base_count.begin(), base_count.end());
  DataFrame df =  DataFrame::create( Named("sum_k") = v1 ,Named("n") = v2);
  return df;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### Correction factor functions ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// DISCONTINUOUS VERSION

//' @title A function to calculate the necessary information to apply the
//' Diggle correction factor with a discontinuous method (sparse)
//' @name corrfactor_discontinuous_sparse
//' @param neighbour_list a list of the neighbours of each node
//' @param events a numeric vector of the node id of each event
//' @param line_list a DataFrame representing the lines of the graph
//' @param bws the kernel bandwidth for each event
//' @param max_depth the maximum recursion depth (after which recursion is stopped)
//' @return a list of dataframes, used to calculate the Diggel correction factor
//' @export
//'
// [[Rcpp::export]]
List corrfactor_discontinuous_sparse(List neighbour_list, NumericVector events, DataFrame line_list, NumericVector bws, int max_depth){

  //extraire le poids des lignes
  arma::vec line_weights =  as<arma::vec>(line_list["weight"]);
  //preparer la matrice de voisinage
  arma::sp_mat edge_mat = make_matrix_sparse(line_list,neighbour_list);
  //preparer la list de DataFrame en sortie
  List list_df;

  //debut des iterations sur les evenements
  int cnt_e = events.length();
  for(int i=0; i < cnt_e; ++i){
    //preparer des vecteurs qui recevront les valeurs de alpha et de distance
    NumericVector edge_ids;
    NumericVector alphas;
    NumericVector dists;
    NumericVector size;
    double bw = bws[i];
    int y = events[i];
    queue<List> data_holder;
    List start = List::create(Named("v") = y ,
                              Named("d") = 0.0,
                              Named("depth") = 0,
                              Named("prev_node") = -1,
                              Named("alpha") = 1.0);
    data_holder.push(start);

    while(data_holder.empty()==FALSE){

      //extraire le 1er cas
      List cas = data_holder.front();
      data_holder.pop();
      int v = cas["v"];
      double alpha = cas["alpha"];
      double d = cas["d"];
      int prev_node = cas["prev_node"];
      int depth = cas["depth"];

      //step1 : find all the neighbours
      IntegerVector neighbours = neighbour_list[v-1];

      //step2 : iterate over the neighbours
      int cnt_n = neighbours.length();
      int new_depth;
      if(cnt_n>2){
        new_depth = depth+1;
      }else{
        new_depth = depth;
      }
      double new_alpha;
      if((prev_node == -1) && (cnt_n > 2)){
        new_alpha = 2.0/(cnt_n);
      }else{
        new_alpha = alpha * (1.0/(cnt_n-1.0));
      }

      //if we have only one neighbour, we must stop
      if(cnt_n>1 or prev_node<=0){
        for(int i=0; i < cnt_n; ++i){
          int v2 = neighbours[i];
          //on ne veut pas revenir en arriere !
          if(v2!=prev_node){
            //find the edge between the two nodes
            int edge_id = edge_mat(v,v2);
            double w = line_weights[edge_id-1];
            //adding the values to the vectors
            edge_ids.push_back(edge_id);
            alphas.push_back(new_alpha);
            dists.push_back(d);
            size.push_back(w);

            //evaluating for the next move
            double d2 = w + d;

            if (d2<bw and new_depth<max_depth){
              List new_cas = List::create(Named("d")=d2,
                                          Named("alpha")=new_alpha,
                                          Named("v") = v2,
                                          Named("prev_node") = v,
                                          Named("depth") = new_depth
              );
              data_holder.push(new_cas);
            }
          }
        }
      }
    }
    DataFrame df = DataFrame::create( Named("edge_id") = edge_ids,
                                      Named("alpha") = alphas,
                                      Named("distances") = dists,
                                      Named("edge_size") = size
    );
    list_df.push_back(df);
  }
  return list_df;
}

//' @title A function to calculate the necessary informations to apply the
//' Diggle correction factor with a discontinuous method
//' @name corrfactor_discontinuous
//' @param neighbour_list a list of the neighbours of each node
//' @param events a numeric vector of the node id of each event
//' @param line_list a DataFrame representing the lines of the graph
//' @param bws the kernel bandwidth for each event
//' @param max_depth the maximum recursion depth (after which recursion is stopped)
//' @return a list of dataframes, used to calculate the Diggel correction factor
//' @export
//
// [[Rcpp::export]]
List corrfactor_discontinuous(List neighbour_list, NumericVector events, DataFrame line_list, NumericVector bws, int max_depth){

  //extraire le poids des lignes
  arma::vec line_weights =  as<arma::vec>(line_list["weight"]);
  //preparer la matrice de voisinage
  IntegerMatrix edge_mat = make_matrix(line_list,neighbour_list);
  //preparer la list de DataFrame en sortie
  List list_df;

  //debut des iterations sur les evenements
  int cnt_e = events.length();
  for(int i=0; i < cnt_e; ++i){
    //preparer des vecteurs qui recevront les valeurs de alpha et de distance
    NumericVector edge_ids;
    NumericVector alphas;
    NumericVector dists;
    NumericVector size;
    double bw = bws[i];
    int y = events[i];
    queue<List> data_holder;
    List start = List::create(Named("v") = y ,
                              Named("d") = 0.0,
                              Named("depth") = 0,
                              Named("prev_node") = -1,
                              Named("alpha") = 1.0);
    data_holder.push(start);

    while(data_holder.empty()==FALSE){

      //extraire le 1er cas
      List cas = data_holder.front();
      data_holder.pop();
      int v = cas["v"];
      double alpha = cas["alpha"];
      double d = cas["d"];
      int prev_node = cas["prev_node"];
      int depth = cas["depth"];

      //step1 : find all the neighbours
      IntegerVector neighbours = neighbour_list[v-1];

      //step2 : iterate over the neighbours
      int cnt_n = neighbours.length();
      int new_depth;
      if(cnt_n>2){
        new_depth = depth+1;
      }else{
        new_depth = depth;
      }
      double new_alpha;
      if((prev_node == -1) && (cnt_n > 2)){
        new_alpha = 2.0/(cnt_n);
      }else{
        new_alpha = alpha * (1.0/(cnt_n-1.0));
      }

      //if we have only one neighbour, we must stop
      if(cnt_n>1 or prev_node<=0){
        for(int i=0; i < cnt_n; ++i){
          int v2 = neighbours[i];
          //on ne veut pas revenir en arriere !
          if(v2!=prev_node){
            //find the edge between the two nodes
            int edge_id = edge_mat(v,v2);
            double w = line_weights[edge_id-1];
            //adding the values to the vectors
            edge_ids.push_back(edge_id);
            alphas.push_back(new_alpha);
            dists.push_back(d);
            size.push_back(w);

            //evaluating for the next move
            double d2 = w + d;

            if (d2<bw and new_depth<max_depth){
              List new_cas = List::create(Named("d")=d2,
                                          Named("alpha")=new_alpha,
                                          Named("v") = v2,
                                          Named("prev_node") = v,
                                          Named("depth") = new_depth
              );
              data_holder.push(new_cas);
            }
          }
        }
      }
    }
    DataFrame df = DataFrame::create( Named("edge_id") = edge_ids,
                                      Named("alpha") = alphas,
                                      Named("distances") = dists,
                                      Named("edge_size") = size
    );
    list_df.push_back(df);
  }
  return list_df;
}



// CONTINUOUS VERSION

//' @title A function to calculate the necessary information to apply the
//' Diggle correction factor with a continuous method (sparse)
//' @name corrfactor_continuous_sparse
//' @param neighbour_list a list of the neighbours of each node
//' @param events a numeric vector of the node id of each event
//' @param line_list a DataFrame representing the lines of the graph
//' @param bws the kernel bandwidth for each event
//' @param max_depth the maximum recursion depth (after which recursion is stopped)
//' @return a list of dataframes, used to calculate the Diggel correction factor
//' @export
//'
// [[Rcpp::export]]
List corrfactor_continuous_sparse(List neighbour_list, NumericVector events, DataFrame line_list, NumericVector bws, int max_depth){

  //extraire le poids des lignes
  arma::vec line_weights =  as<arma::vec>(line_list["weight"]);
  //preparer la matrice de voisinage
  arma::sp_mat edge_mat = make_matrix_sparse(line_list,neighbour_list);
  //preparer la list de DataFrame en sortie
  List list_df;

  //debut des iterations sur les evenements
  int cnt_e = events.length();
  for(int i=0; i < cnt_e; ++i){
    //preparer des vecteurs qui recevront les valeurs de alpha et de distance
    NumericVector edge_ids;
    NumericVector alphas;
    NumericVector dists;
    NumericVector size;
    double bw = bws[i];
    int y = events[i];
    queue<List> data_holder;
    List start = List::create(Named("v") = y ,
                              Named("d") = 0.0,
                              Named("depth") = 0,
                              Named("prev_node") = -1,
                              Named("alpha") = 1.0);
    data_holder.push(start);

    while(data_holder.empty()==FALSE){

      //extraire le 1er cas
      List cas = data_holder.front();
      data_holder.pop();
      int v = cas["v"];
      double alpha = cas["alpha"];
      double d = cas["d"];
      int prev_node = cas["prev_node"];
      int depth = cas["depth"];

      //step1 : find all the neighbours
      IntegerVector neighbours = neighbour_list[v-1];
      int cnt_n = neighbours.length();

      int new_depth;
      if(cnt_n>2){
        new_depth = depth+1;
      }else{
        new_depth = depth;
      }
      //if we have only one neighbour, we must stop
      if(cnt_n>1 or prev_node<=0){
        for(int i=0; i < cnt_n; ++i){
          int v2 = neighbours[i];
          //find the edge between the two nodes
          int edge_id = edge_mat(v,v2);
          double w = line_weights[edge_id-1];
          double new_alpha;


          if(v2!=prev_node){
            //----- cas simple, on ne revient pas en arriere !
            new_alpha = alpha * (2.0/cnt_n);
            //adding the values in the vectors !
            edge_ids.push_back(edge_id);
            alphas.push_back(new_alpha);
            dists.push_back(d);
            size.push_back(w);
            //and preparing the next step !
            double d2 = w + d;
            if (d2<bw and new_depth<max_depth){
              List new_cas = List::create(Named("d")=d2,
                                          Named("alpha")=new_alpha,
                                          Named("v") = v2,
                                          Named("prev_node") = v,
                                          Named("depth") = new_depth
              );
              data_holder.push(new_cas);

            }
          }// cas special pour le back fire !
          else{
            //on ne peut faire un backfire qu'a une vraie intersection
            if(cnt_n>2){
              IntegerVector l_neighbours = neighbour_list[v2-1];
              int n2 = l_neighbours.length();
              double p2 = (n2-2.0)/n2;
              new_alpha = -1.0 * alpha * p2;
              //adding the values in the vectors !
              edge_ids.push_back(edge_id);
              alphas.push_back(new_alpha);
              dists.push_back(d);
              size.push_back(w);
              //and preparing the next step !
              double d2 = w + d;
              if (d2<bw and new_depth<max_depth){
                List new_cas = List::create(Named("d")=d2,
                                            Named("alpha")=new_alpha,
                                            Named("v") = v2,
                                            Named("prev_node") = v,
                                            Named("depth") = new_depth
                );
                data_holder.push(new_cas);
              }
            }
          }
        }
      }
    }
    DataFrame df = DataFrame::create( Named("edge_id") = edge_ids,
                                      Named("alpha") = alphas,
                                      Named("distances") = dists,
                                      Named("edge_size") = size);
    list_df.push_back(df);
  }
  return list_df;
}

//' @title A function to calculate the necessary information to apply the
//' Diggle correction factor with a continuous method
//' @name corrfactor_continuous
//' @param neighbour_list a list of the neighbours of each node
//' @param events a numeric vector of the node id of each event
//' @param line_list a DataFrame representing the lines of the graph
//' @param bws the kernel bandwidth for each event
//' @param max_depth the maximum recursion depth (after which recursion is stopped)
//' @return a list of dataframes, used to calculate the Diggel correction factor
//' @export
//'
// [[Rcpp::export]]
List corrfactor_continuous(List neighbour_list, NumericVector events, DataFrame line_list, NumericVector bws, int max_depth){

  //extraire le poids des lignes
  arma::vec line_weights =  as<arma::vec>(line_list["weight"]);
  //preparer la matrice de voisinage
  IntegerMatrix edge_mat = make_matrix(line_list,neighbour_list);
  //preparer la list de DataFrame en sortie
  List list_df;

  //debut des iterations sur les evenements
  int cnt_e = events.length();
  for(int i=0; i < cnt_e; ++i){
    //preparer des vecteurs qui recevront les valeurs de alpha et de distance
    NumericVector edge_ids;
    NumericVector alphas;
    NumericVector dists;
    NumericVector size;
    double bw = bws[i];
    int y = events[i];
    queue<List> data_holder;
    List start = List::create(Named("v") = y ,
                              Named("d") = 0.0,
                              Named("depth") = 0,
                              Named("prev_node") = -1,
                              Named("alpha") = 1.0);
    data_holder.push(start);

    while(data_holder.empty()==FALSE){

      //extraire le 1er cas
      List cas = data_holder.front();
      data_holder.pop();
      int v = cas["v"];
      double alpha = cas["alpha"];
      double d = cas["d"];
      int prev_node = cas["prev_node"];
      int depth = cas["depth"];

      //step1 : find all the neighbours
      IntegerVector neighbours = neighbour_list[v-1];
      int cnt_n = neighbours.length();

      int new_depth;
      if(cnt_n>2){
        new_depth = depth+1;
      }else{
        new_depth = depth;
      }
      //if we have only one neighbour, we must stop
      if(cnt_n>1 or prev_node<=0){
        for(int i=0; i < cnt_n; ++i){
          int v2 = neighbours[i];
          //find the edge between the two nodes
          int edge_id = edge_mat(v,v2);
          double w = line_weights[edge_id-1];
          double new_alpha;


          if(v2!=prev_node){
            //----- cas simple, on ne revient pas en arriere !
            new_alpha = alpha * (2.0/cnt_n);
            //adding the values in the vectors !
            edge_ids.push_back(edge_id);
            alphas.push_back(new_alpha);
            dists.push_back(d);
            size.push_back(w);
            //and preparing the next step !
            double d2 = w + d;
            if (d2<bw and new_depth<max_depth){
              List new_cas = List::create(Named("d")=d2,
                                          Named("alpha")=new_alpha,
                                          Named("v") = v2,
                                          Named("prev_node") = v,
                                          Named("depth") = new_depth
              );
              data_holder.push(new_cas);

            }
          }// cas special pour le back fire !
          else{
            //on ne peut faire un backfire qu'a une vraie intersection
            if(cnt_n>2){
              IntegerVector l_neighbours = neighbour_list[v2-1];
              int n2 = l_neighbours.length();
              double p2 = (n2-2.0)/n2;
              new_alpha = -1.0 * alpha * p2;
              //adding the values in the vectors !
              edge_ids.push_back(edge_id);
              alphas.push_back(new_alpha);
              dists.push_back(d);
              size.push_back(w);
              //and preparing the next step !
              double d2 = w + d;
              if (d2<bw and new_depth<max_depth){
                List new_cas = List::create(Named("d")=d2,
                                            Named("alpha")=new_alpha,
                                            Named("v") = v2,
                                            Named("prev_node") = v,
                                            Named("depth") = new_depth
                );
                data_holder.push(new_cas);
              }
            }
          }
        }
      }
    }
    DataFrame df = DataFrame::create( Named("edge_id") = edge_ids,
                                      Named("alpha") = alphas,
                                      Named("distances") = dists,
                                      Named("edge_size") = size);
    list_df.push_back(df);
  }
  return list_df;
}



//########################################################################
//#### bandwidth selection by likelihood loo cross validation #####
//########################################################################

// a simple function to find the index of the first occurence of value in a numeric vector
int get_first_index(NumericVector v1, double x){
  int i;
  for( i = 0; i < v1.size(); ++i) {
    if(v1[i] == x){
      return i;
    }
  }
  return -1;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Here we have the working functions for loo with continuous kernel
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title The recursive function to calculate continuous NKDE likelihood cv
//' @name esc_kernel_loo
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param samples_k a numeric vector of the actual kernel values, updated at
//' each recursion
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param v the actual node to consider for the recursion (int)
//' @param v1 the connected node to consider for the recursion (int)
//' @param l1 the edge connecting v and v1 (int)
//' @param d the actual distance traveled before the recursion
//' @param alpha the actual alpha value before the recursion
//' @param bw the kernel bandwidth
//' @param line_weights a vector with the length of the edges
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a vector with the kernel values calculated for each samples from
//' the first node given
NumericVector esc_kernel_loo(fptros kernel_func, NumericVector samples_k, List neighbour_list, arma::sp_mat edge_mat, int v, int v1, int l1, double d,double alpha, double bw, NumericVector line_weights, NumericVector events, int depth, int max_depth){

  //mettre a jour d
  double d2 = line_weights[l1-1] + d;
  if((bw>d2) && (depth < max_depth)){

    // on veut mettre a jour les valeurs dans samples_k
    int index = get_first_index(events,v1);
    if( index > -1){
      samples_k[index] = samples_k[index] + kernel_func(d2,bw)*alpha;

    }

    //on veut trouver toutes les lignes emannant de v (Lv)
    IntegerVector v_neighbours = neighbour_list[v1-1];
    int n = v_neighbours.length();
    int new_depth;
    //updating depth
    if(n>2){
      new_depth = depth+1;
    }else{
      new_depth = depth+0;
    }
    if(n>1){
      for(int j=0; j < n; ++j){
        //int li = Lv[j];
        int vi = v_neighbours[j];
        int li = edge_mat(v1,vi);
        if(li==l1){
          //we must backfire only if we have a true intersection
          if(n>2){
            double p2 = (2.0-n)/n;
            double n_alpha = alpha * p2;
            samples_k = esc_kernel_loo(kernel_func, samples_k, neighbour_list, edge_mat, v1,vi, li, d2, n_alpha, bw, line_weights, events,new_depth, max_depth);
          }
        }else{
          double n_alpha = alpha * (2.0/n);
          samples_k = esc_kernel_loo(kernel_func, samples_k, neighbour_list, edge_mat, v1,vi, li, d2, n_alpha, bw, line_weights, events, new_depth, max_depth);
        };
      };
    };
  };
  return samples_k;
}


//' @title The exported function to calculate continuous NKDE likelihood cv
//' @name get_loo_values_continuous
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param samples a vector with the node index of the samples
//' @param sweights a vector with the weight of each sample
//' @param events a vector with the node index of the events
//' @param weights a vector with the weights of the events
//' @param bws a vector with the bandwidths
//' @param kernel_name the name of the kernel function to use
//' @param line_list a dataframe representing the edges of the network
//' @param max_depth the maximum recursion depth
//' @return a vector with the kernel values calculated for each event
//' @export
//'
// [[Rcpp::export]]
DataFrame get_loo_values_continuous(List neighbour_list,
                                    NumericVector samples, NumericVector sweights,
                                    NumericVector events, NumericVector weights,
                                    NumericVector bws, std::string kernel_name, DataFrame line_list,
                                    int max_depth){

  //selecting the kernel function
  fptros kernel_func = select_kernelos(kernel_name);

  //step0 extract the columns of the dataframe
  NumericVector line_weights = line_list["weight"];

  //step 1 : mettre toutes les valeurs a 0
  //idee generale : calculer notre base_k comme toujours
  //mais a la fin du calcul sur un event, recuperer la valeur sur cet
  //event pour le loo. + ne calculer les densites que sur des evenements
  int Ne = samples.length();
  double startval = kernel_func(0,bws[0]);
  NumericVector base_k = sweights * startval;
  NumericVector loo_k = sweights * startval;
  NumericVector samples_k(Ne);

  //calculer le dictionnaire des lignes
  arma::sp_mat edge_mat = make_matrix_sparse(line_list,neighbour_list);
  //step2 : iterer sur chaque event
  int cnt_e = events.length()-1;
  for(int i=0; i <= cnt_e; ++i){
    //preparer les differentes valeurs de departs pour l'event y
    int y = events[i]; // NOTE y est le code du noeud correspondant a l'evenement a la place y dans son vecteur
    double w = weights[i];
    double bw = bws[i];
    //on veut trouver toutes les voisins emannant de y
    IntegerVector y_neighbours = neighbour_list[y-1];
    int cnt_y = y_neighbours.length();
    double alpha = 2.0/cnt_y;
    for(int j=0; j < cnt_y; ++j){
      //preparing all values for this loop
      int vi = y_neighbours[j];
      int li = edge_mat(y,vi);
      double d = 0.0 ;
      int depth = 0 ;
      // launching recursion
      samples_k.fill(0.0); // samples_k here will be the effect of the event on all events
      samples_k = esc_kernel_loo(kernel_func, samples_k, neighbour_list, edge_mat ,y,vi,li,d,alpha,bw, line_weights, samples, depth,max_depth);
      base_k = base_k + samples_k*w;
      // adding to loo_k if the event is in samples !
      int idx = get_first_index(samples,y);
      if(idx >= 0){
        loo_k[idx] = loo_k[idx] + samples_k[idx]*w;
      };

    };
  };

  DataFrame df =  DataFrame::create( Named("sum_k") = base_k , Named("loo") = loo_k);
  return df;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Here we have the working functions for loo with discontinuous kernel
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title The function to calculate discontinuous NKDE likelihood cv
//' @name esd_kernel_loo
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param samples_k a numeric vector of the actual kernel values, updated at
//' each recursion
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param v the actual node to consider for the recursion (int)
//' @param bw the kernel bandwidth
//' @param line_weights a vector with the length of the edges
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a vector with the kernel values calculated for each samples from
//' the first node given
NumericVector esd_kernel_loo(fptros kernel_func, arma::sp_mat edge_mat, NumericVector events,
                               List neighbour_list ,int v, double bw,
                               NumericVector line_weights, int depth, int max_depth){
  //step0 : generate the queue
  queue <List> data_holder;
  NumericVector kvalues(events.length());
  kvalues.fill(0.0);
  //step1 : generate the first case
  List cas1 = List::create(Named("d")=0.0,
                           Named("alpha")=1.0,
                           Named("v") = v,
                           Named("prev_node") = -999,
                           Named("depth") = 0
  );
  data_holder.push(cas1);
  //lancement des iterations
  while(data_holder.empty()==FALSE){
    //unpacking (imagine some loop unrolling here with a function to deal with.)
    List cas = data_holder.front();
    data_holder.pop();
    int v = cas["v"];
    int depth = cas["depth"];
    double d = cas["d"];
    double alpha = cas["alpha"];
    int prev_node = cas["prev_node"];

    //step1 : find all the neighbours
    IntegerVector neighbours = neighbour_list[v-1];

    //step2 : iterate over the neighbours
    int cnt_n = neighbours.length();
    int new_depth;
    if(cnt_n>2){
      new_depth = depth+1;
    }else{
      new_depth = depth;
    }

    double new_alpha;
    if((prev_node < 0)  && (cnt_n > 2)){
      new_alpha = 2.0/(cnt_n);
    }else{
      new_alpha = alpha * (1.0/(cnt_n-1.0));
    }

    //if we have only one neighbour, we must stop
    if(cnt_n>1 or prev_node<=0){
      for(int i=0; i < cnt_n; ++i){
        int v2 = neighbours[i];
        //on ne veut pas revenir en arriere !
        if(v2!=prev_node){
          //find the edge between the two nodes
          int edge_id = edge_mat(v,v2);
          double d2 = line_weights[edge_id-1] + d;

          //est ce que v2 est un evenement ?
          int index = get_first_index(events,v2);
          if(index >=0){
            kvalues[index] = kvalues[index] + kernel_func(d2,bw)*new_alpha;
          }
          //evaluating for the next move


          if (d2<bw and new_depth<max_depth){
            List new_cas = List::create(Named("d")=d2,
                                        Named("alpha")=new_alpha,
                                        Named("v") = v2,
                                        Named("prev_node") = v,
                                        Named("depth") = new_depth
            );
            data_holder.push(new_cas);
          }
        }
      }
    }
  }

  return kvalues;

}




//' @title The exported function to calculate discontinuous NKDE likelihood cv
//' @name get_loo_values_discontinuous
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param samples a vector with the node index of the samples
//' @param sweights a vector with the weight of each sample
//' @param events a vector with the node index of the events
//' @param weights a vector with the weights of the events
//' @param bws a vector with the bandwidths
//' @param kernel_name the name of the kernel function to use
//' @param line_list a dataframe representing the edges of the network
//' @param max_depth the maximum recursion depth
//' @return a vector with the kernel values calculated for each event
//' @export
//'
// [[Rcpp::export]]
NumericVector get_loo_values_discontinuous(List neighbour_list, NumericVector samples, NumericVector sweights,
                                           NumericVector events, NumericVector weights,
                                           NumericVector bws, std::string kernel_name, DataFrame line_list,
                                           int max_depth){

  //selecting the kernel function
  fptros kernel_func = select_kernelos(kernel_name);

  //step0 extract the columns of the dataframe
  NumericVector line_weights = line_list["weight"];

  //step 1 : mettre toutes les valeurs a 0
  double startValue = kernel_func(0,bws[0]);
  NumericVector base_k = sweights * startValue;


  //calculer la matrice des lignes
  //IntegerMatrix edge_mat = make_matrix(line_list,neighbour_list);
  arma::sp_mat edge_mat = make_matrix_sparse(line_list,neighbour_list);
  int depth = 0;

  //step2 : iterer sur chaque event
  int cnt_e = events.length()-1;
  for(int i=0; i <= cnt_e; ++i){
    //preparer les differentes valeurs de departs pour l'event y
    int y = events[i];
    double w = weights[i];
    double bw = bws[i];
    // launching recursion
    NumericVector k = esd_kernel_loo(kernel_func,edge_mat, samples, neighbour_list ,y,bw, line_weights, depth,max_depth);
    // getting the actual base_k values (summed at each iteration)
    base_k += (k*w);
  };
  return base_k;
}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Here we have the functions in development for loo with simple kernel
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title The function to calculate simple NKDE likelihood cv
//' @name esd_kernel_loo
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param samples_k a numeric vector of the actual kernel values, updated at
//' each recursion
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param v the actual node to consider for the recursion (int)
//' @param bw the kernel bandwidth
//' @param line_weights a vector with the length of the edges
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a vector with the kernel values calculated for each samples from
//' the first node given
NumericVector ess_kernel_loo(fptros kernel_func, arma::sp_mat edge_mat, NumericVector events,
                             List neighbour_list ,int v, double bw,
                             NumericVector line_weights, int depth, int max_depth){
  //step0 : generate the queue
  queue <List> data_holder;
  NumericVector kvalues(events.length());
  kvalues.fill(0.0);
  //step1 : generate the first case
  List cas1 = List::create(Named("d")=0.0,
                           Named("v") = v,
                           Named("prev_node") = -999,
                           Named("depth") = 0
  );
  data_holder.push(cas1);
  //lancement des iterations
  while(data_holder.empty()==FALSE){
    //unpacking (imagine some loop unrolling here with a function to deal with.)
    List cas = data_holder.front();
    data_holder.pop();
    int v = cas["v"];
    int depth = cas["depth"];
    double d = cas["d"];
    int prev_node = cas["prev_node"];

    //step1 : find all the neighbours
    IntegerVector neighbours = neighbour_list[v-1];

    //step2 : iterate over the neighbours
    int cnt_n = neighbours.length();
    int new_depth;
    if(cnt_n>2){
      new_depth = depth+1;
    }else{
      new_depth = depth;
    }

    //if we have only one neighbour, we must stop
    if(cnt_n>1 or prev_node<=0){
      for(int i=0; i < cnt_n; ++i){
        int v2 = neighbours[i];
        //on ne veut pas revenir en arriere !
        if(v2!=prev_node){
          //find the edge between the two nodes
          int edge_id = edge_mat(v,v2);
          double d2 = line_weights[edge_id-1] + d;

          //est ce que v2 est un evenement pour lequel on doit garder la valeur ?
          int index = get_first_index(events,v2);
          if(index >=0){
            kvalues[index] = kvalues[index] + kernel_func(d2,bw);
          }
          //evaluating for the next move


          if (d2<bw and new_depth<max_depth){
            List new_cas = List::create(Named("d")=d2,
                                        Named("v") = v2,
                                        Named("prev_node") = v,
                                        Named("depth") = new_depth
            );
            data_holder.push(new_cas);
          }
        }
      }
    }
  }

  return kvalues;

}




//' @title The exported function to calculate simple NKDE likelihood cv
//' @name get_loo_values_simple
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param samples a vector with the node index of the samples
//' @param sweights a vector with the weight of each sample
//' @param events a vector with the node index of the events
//' @param weights a vector with the weights of the events
//' @param bws a vector with the bandwidths
//' @param kernel_name the name of the kernel function to use
//' @param line_list a dataframe representing the edges of the network
//' @param max_depth the maximum recursion depth
//' @return a vector with the kernel values calculated for each event
//' @export
//'
// [[Rcpp::export]]
NumericVector get_loo_values_simple(List neighbour_list, NumericVector samples, NumericVector sweights,
                                    NumericVector events, NumericVector weights,
                                    NumericVector bws, std::string kernel_name,
                                    DataFrame line_list, int max_depth){

  //selecting the kernel function
  fptros kernel_func = select_kernelos(kernel_name);

  //step0 extract the columns of the dataframe
  NumericVector line_weights = line_list["weight"];

  //step 1 : mettre toutes les valeurs a 0

  double startValue = kernel_func(0,bws[0]);
  NumericVector base_k = sweights * startValue;


  //calculer la matrice des lignes
  //IntegerMatrix edge_mat = make_matrix(line_list,neighbour_list);
  arma::sp_mat edge_mat = make_matrix_sparse(line_list,neighbour_list);
  int depth = 0;

  //step2 : iterer sur chaque event
  int cnt_e = events.length()-1;
  for(int i=0; i <= cnt_e; ++i){
    //preparer les differentes valeurs de departs pour l'event y
    int y = events[i];
    double w = weights[i];
    double bw = bws[i];
    // launching recursion
    NumericVector k = ess_kernel_loo(kernel_func,edge_mat, samples, neighbour_list ,y,bw, line_weights, depth,max_depth);
    // getting the actual base_k values (summed at each iteration)
    base_k += (k*w);
  };
  return base_k;
}
