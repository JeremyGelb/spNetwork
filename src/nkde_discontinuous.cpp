#include "spNetwork.h"
#include "base_kernel_funtions.h"
#include "matrices_functions.h"


//#####################################################################################
// ######################  THE WORKER FUNCTIONS  ######################################
//#####################################################################################


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
//' @param samples_coords a matrix with the X and Y coordinates of the samples
//' @param nodes_coords a matrix with the X and Y coordinates of the nodes
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a vector with the kernel values calculated for each samples from
//' the first node given
//' @keywords internal
arma::vec esd_kernel_rcpp_arma_sparse(fptr kernel_func, arma::sp_imat &edge_mat,
                                      List &neighbour_list ,int v, double bw,
                                      arma::vec &line_weights,
                                      arma::ivec &samples_edgeid,
                                      arma::mat &samples_coords,
                                      arma::mat &nodes_coords,
                                      int depth, int max_depth){
  //step0 : generate the queue
  //queue <List> data_holder;

  struct acase {
    double d;
    double alpha;
    int v;
    int prev_node;
    int depth;
  };

  std::vector<acase> data_holder;

  arma::vec samples_k(samples_edgeid.n_elem);
  samples_k.fill(0.0);
  //step1 : generate the first case
  acase cas1 = {0.0,1.0,v,-999,0};
  data_holder.push_back(cas1);

  int new_depth;
  int cnt_n;
  double new_alpha;

  //lancement des iterations
  while(data_holder.empty()==FALSE){
    //unpacking (imagine some loop unrolling here with a function to deal with.)
    acase cas = data_holder.back();
    data_holder.pop_back();

    int v = cas.v;

    //step1 : find all the neighbours
    IntegerVector neighbours = neighbour_list[v-1];

    //step2 : iterate over the neighbours
    cnt_n = neighbours.length();
    if(cnt_n>2){
      new_depth = cas.depth+1;
    }else{
      new_depth = cas.depth;
    }

    if((cas.prev_node < 0)  && (cnt_n > 2)){
      new_alpha = 2.0/(cnt_n);
    }else if((cas.prev_node < 0)  && (cnt_n == 1)){
      new_alpha = 1;
    }else{
      new_alpha = cas.alpha * (1.0/(cnt_n-1.0));
    }

    //if we have only one neighbour, we must stop
    if(cnt_n>1 or cas.prev_node<=0){
      for(int i=0; i < cnt_n; ++i){
        int v2 = neighbours[i];
        //on ne veut pas revenir en arriere !
        if(v2!=cas.prev_node){
          //find the edge between the two nodes
          int edge_id = edge_mat(v,v2);
          //find the samples on that edge
          arma::uvec test = arma::find(samples_edgeid==edge_id);

          arma::mat sampling_coords = samples_coords.rows(test);


          //extracting the X and Y coordinates of the starting node

          arma::colvec x_dists =  arma::sqrt(arma::sum(arma::pow(sampling_coords.each_row() - nodes_coords.row((v - 1)),2),1)) + cas.d;

          //step3 calculating the values of the new kernel
          arma::vec new_k = kernel_func(x_dists,bw)*new_alpha;
          samples_k.elem(test) += new_k;

          //evaluating for the next move
          double d2 = line_weights[edge_id-1] + cas.d;

          if (d2<bw and new_depth<max_depth){
            acase new_cas = {d2,new_alpha,v2,v,new_depth};
            data_holder.push_back(new_cas);
            //data_holder.push_back((struct acase){d2,new_alpha,v2,cas.v,new_depth});
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
//' @keywords internal
arma::vec esd_kernel_rcpp_arma(fptr kernel_func, IntegerMatrix& edge_mat,
                               List& neighbour_list ,int v, double bw,
                               arma::vec& line_weights, arma::vec& samples_edgeid,
                               arma::vec& samples_x, arma::vec& samples_y,
                               arma::vec& nodes_x, arma::vec& nodes_y, int depth, int max_depth){
  //step0 : generate the queue
  //queue <List> data_holder;
  std::vector <List> data_holder;
  arma::vec samples_k(samples_x.n_elem);
  //samples_k.fill(0.0);
  //step1 : generate the first case
  List cas1 = List::create(Named("d")=0.0,
                           Named("alpha")=1.0,
                           Named("v") = v,
                           Named("prev_node") = -999,
                           Named("depth") = 0
  );
  data_holder.push_back(cas1);
  //lancement des iterations
  while(data_holder.empty()==FALSE){
    //unpacking (imagine some loop unrolling here with a function to deal with.)
    List cas = data_holder.back();
    data_holder.pop_back();
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
    }else if((prev_node < 0)  && (cnt_n == 1)){
      new_alpha = 1;
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
            data_holder.push_back(new_cas);
          }
        }
      }
    }
  }

  return samples_k;

}


//#####################################################################################
// #####################  THE EXECUTION FUNCTIONS  ####################################
//#####################################################################################

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
//' @param div The divisor to use for the kernel. Must be "n" (the number of events within the radius around each sampling point), "bw" (the bandwidth) "none" (the simple sum).
//' @return a DataFrame with two columns : the kernel values (sum_k) and the number of events for each sample (n)
//' @export
//' @keywords internal
// [[Rcpp::export]]
DataFrame discontinuous_nkde_cpp_arma_sparse(List neighbour_list, NumericVector events,
                                             NumericVector weights, DataFrame samples,
                                             NumericVector bws, std::string kernel_name,
                                             DataFrame nodes, DataFrame line_list,
                                             int max_depth, bool verbose, std::string div = "bw"){

  //selecting the kernel function
  fptr kernel_func = select_kernel(kernel_name);

  //step0 extract the columns of the dataframe
  arma::vec line_weights =  as<arma::vec>(line_list["weight"]);
  arma::ivec samples_edgeid = as<arma::ivec>(samples["edge_id"]);


  arma::mat samples_coords(samples.nrows(),2);
  samples_coords.col(0) = as<arma::vec>(samples["X_coords"]);
  samples_coords.col(1) = as<arma::vec>(samples["Y_coords"]);

  arma::mat nodes_coords(nodes.nrows(),2);
  nodes_coords.col(0) = as<arma::vec>(nodes["X_coords"]);
  nodes_coords.col(1) = as<arma::vec>(nodes["Y_coords"]);

  //step 1 : mettre toutes les valeurs a 0
  arma::vec base_k(samples.nrow());
  arma::vec base_count(samples.nrow());
  arma::vec count(samples.nrow());

  // base_k.fill(0.0);
  // base_count.fill(0.0);
  // count.fill(0.0);


  //calculer la matrice des lignes
  //arma::sp_mat edge_mat = make_matrix_sparse(line_list,neighbour_list);
  arma::sp_imat edge_mat = make_imatrix_sparse(line_list, neighbour_list);
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
    arma::vec k = esd_kernel_rcpp_arma_sparse(kernel_func,edge_mat, neighbour_list ,y, bw, line_weights, samples_edgeid, samples_coords, nodes_coords, depth,max_depth);
    // getting the actual base_k values (summed at each iteration)
    if(div == "bw"){
      base_k += ((k*w))/bw;
    }else{
      base_k += (k*w);
    }
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
//' @param div The divisor to use for the kernel. Must be "n" (the number of events within the radius around each sampling point), "bw" (the bandwidth) "none" (the simple sum).
//' @return a DataFrame with two columns : the kernel values (sum_k) and the number of events for each sample (n)
//' @export
//' @keywords internal
// [[Rcpp::export]]
DataFrame discontinuous_nkde_cpp_arma(List neighbour_list, NumericVector events,
                                      NumericVector weights, DataFrame samples,
                                      NumericVector bws, std::string kernel_name,
                                      DataFrame nodes, DataFrame line_list,
                                      int max_depth, bool verbose, std::string div = "bw"){

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
  // base_k.fill(0.0);
  // base_count.fill(0.0);
  // count.fill(0.0);


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
    if(div == "bw"){
      base_k += ((k*w))/bw;
    }else{
      base_k += (k*w);
    }
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


//#####################################################################################
// ##################  THE EXECUTION FUNCTIONS FOT TNKDE  #############################
//#####################################################################################

//' @title The main function to calculate discontinuous NKDE (ARMA and sparse matrix)
//' @name tnkdediscontinuousfunctionsparse
//' @param neighbour_list a list of the neighbours of each node
//' @param events a numeric vector of the node id of each event
//' @param events_time a numeric vector with the time for the events
//' @param weights a numeric vector of the weight of each event
//' @param samples a DataFrame of the samples (with spatial coordinates and belonging edge)
//' @param samples_time a NumericVector indicating when to do the samples
//' @param bws_net the network kernel bandwidths for each event
//' @param kernel_name the name of the kernel function to use
//' @param nodes a DataFrame representing the nodes of the graph (with spatial coordinates)
//' @param line_list a DataFrame representing the lines of the graph
//' @param max_depth the maximum recursion depth (after which recursion is stopped)
//' @param verbose a boolean indicating if the function must print its progress
//' @param div a string indicating how to standardize the kernel values
//' @return a List with two matrices: the kernel values (sum_k) and the number of events for each sample (n)
//' @export
//' @keywords internal
// [[Rcpp::export]]
List discontinuous_tnkde_cpp_arma_sparse(List neighbour_list,
                                         IntegerVector events, NumericVector weights, NumericVector events_time,
                                         DataFrame samples, arma::vec samples_time,
                                         NumericVector bws_net,
                                         NumericVector bws_time,
                                         std::string kernel_name, DataFrame nodes, DataFrame line_list,
                                         int max_depth, bool verbose, std::string div = "bw"){

  //selecting the kernel function
  fptr kernel_func = select_kernel(kernel_name);

  //step0 extract the columns of the dataframe
  arma::vec line_weights =  as<arma::vec>(line_list["weight"]);
  arma::ivec samples_edgeid = as<arma::ivec>(samples["edge_id"]);

  arma::mat samples_coords(samples.nrows(),2);
  samples_coords.col(0) = as<arma::vec>(samples["X_coords"]);
  samples_coords.col(1) = as<arma::vec>(samples["Y_coords"]);

  arma::mat nodes_coords(nodes.nrows(),2);
  nodes_coords.col(0) = as<arma::vec>(nodes["X_coords"]);
  nodes_coords.col(1) = as<arma::vec>(nodes["Y_coords"]);

  //step 1 : set all values to 0
  // NOTE : we will produce matrices here (tnkde)
  arma::mat base_k(samples.nrows(), samples_time.n_elem);
  arma::mat base_count(samples.nrows(), samples_time.n_elem);
  arma::mat count(samples.nrows(), samples_time.n_elem);

  // NOTE It is possible that we must recalculate the density for the same vertex
  // we will store them instead
  std::map<int, int> event_counter = count_values_intvec(events);
  std::map<int, arma::vec> saved_values;

  //calculer la matrice des lignes
  //arma::sp_mat edge_mat = make_matrix_sparse(line_list,neighbour_list);
  arma::sp_imat edge_mat = make_imatrix_sparse(line_list, neighbour_list);
  int depth = 0;
  //step2 : iterer sur chaque event
  int cnt_e = events.length()-1;
  Progress p(cnt_e, verbose);
  for(int i=0; i <= cnt_e; ++i){
    p.increment(); // update progress
    //preparer les differentes valeurs de departs pour l'event y
    int y = events[i];
    double w = weights[i];
    double bw_net = bws_net[i];
    double bw_time = bws_time[i];
    double et = events_time[i];
    arma::vec k;
    // calculating the density value on the network
    if(event_counter[y] > 1){
      // here we have a dupplicated event, we must either use the saved nkde or calculate it and save it latter
      if(map_contains_key(saved_values,y)){
        // we already have the value !
        k = saved_values[y];
      }else{
        // we have not yet the value
        k = esd_kernel_rcpp_arma_sparse(kernel_func,edge_mat, neighbour_list ,y,bw_net, line_weights, samples_edgeid, samples_coords, nodes_coords, depth,max_depth);
        saved_values[y] = k;
      }
    }else{
      // this is not a duplicate event, no need to store the value in memory
      k = esd_kernel_rcpp_arma_sparse(kernel_func,edge_mat, neighbour_list ,y, bw_net, line_weights, samples_edgeid, samples_coords, nodes_coords, depth,max_depth);
    }
    // ok, now we must calculate the temporal densities
    arma::vec temporal_density = kernel_func(arma::abs(samples_time - et), bw_time);

    //applying the scaling here if required
    if(div == "bw"){
      k = k / bw_net;
      temporal_density = temporal_density / bw_time;
    }

    // and create an awesome spatio-temporal matrix
    arma::mat k_mat(samples.nrows(), samples_time.n_elem);

    int mj = temporal_density.n_elem;
    for(int j = 0; j < mj; j++){
      k_mat.col(j) = k * temporal_density[j];
    }
    // getting the actual base_k values (summed at each iteration)
    base_k += (k_mat*w);
    // calculating the new value
    count.fill(0.0);
    arma::umat ids = arma::find(k_mat>0);
    count.elem(ids).fill(1.0);
    base_count += count*w;
  };

  // standardise the result if div = n
  if(div == "n"){
    base_k = base_k / base_count;
  }

  List results = List::create(
    Named("sum_k") = base_k,
    Named("n") = base_count );
  return results;
}



//' @title The main function to calculate discontinuous NKDE (ARMA and Integer matrix)
//' @name tnkdediscontinuousfunction
//' @param neighbour_list a list of the neighbours of each node
//' @param events a numeric vector of the node id of each event
//' @param events_time a numeric vector with the time for the events
//' @param weights a numeric vector of the weight of each event
//' @param samples a DataFrame of the samples (with spatial coordinates and belonging edge)
//' @param samples_time a NumericVector indicating when to do the samples
//' @param bws_net the network kernel bandwidths for each event
//' @param kernel_name the name of the kernel function to use
//' @param nodes a DataFrame representing the nodes of the graph (with spatial coordinates)
//' @param line_list a DataFrame representing the lines of the graph
//' @param max_depth the maximum recursion depth (after which recursion is stopped)
//' @param verbose a boolean indicating if the function must print its progress
//' @param div a string indicating how to standardize the kernel values
//' @return a List with two matrices: the kernel values (sum_k) and the number of events for each sample (n)
//' @export
//' @keywords internal
// [[Rcpp::export]]
List discontinuous_tnkde_cpp_arma(List neighbour_list,
                                  IntegerVector events, NumericVector weights, NumericVector events_time,
                                  DataFrame samples, arma::vec samples_time,
                                  NumericVector bws_net,
                                  NumericVector bws_time,
                                  std::string kernel_name, DataFrame nodes, DataFrame line_list,
                                  int max_depth, bool verbose, std::string div = "bw"){

  //selecting the kernel function
  fptr kernel_func = select_kernel(kernel_name);

  //step0 extract the columns of the dataframe
  arma::vec line_weights =  as<arma::vec>(line_list["weight"]);
  arma::vec samples_edgeid = as<arma::vec>(samples["edge_id"]);
  arma::vec samples_x = as<arma::vec>(samples["X_coords"]);
  arma::vec samples_y = as<arma::vec>(samples["Y_coords"]);
  arma::vec nodes_x = as<arma::vec>(nodes["X_coords"]);
  arma::vec nodes_y = as<arma::vec>(nodes["Y_coords"]);

  //step 1 : set all values to 0
  // NOTE : we will produce matrices here (tnkde)
  arma::mat base_k(samples_x.n_elem, samples_time.n_elem);
  arma::mat base_count(samples_x.n_elem, samples_time.n_elem);
  arma::mat count(samples_x.n_elem, samples_time.n_elem);

  // NOTE It is possible that we must recalculate the density for the same vertex
  // we will store them instead
  std::map<int, int> event_counter = count_values_intvec(events);
  std::map<int, arma::vec> saved_values;

  //calculer la matrice des lignes
  //IntegerMatrix edge_mat = make_matrix(line_list,neighbour_list);
  IntegerMatrix edge_mat = make_matrix(line_list, neighbour_list);
  int depth = 0;


  //step2 : iterer sur chaque event
  int cnt_e = events.length()-1;
  Progress p(cnt_e, verbose);
  for(int i=0; i <= cnt_e; ++i){
    p.increment(); // update progress
    //preparer les differentes valeurs de departs pour l'event y
    int y = events[i];
    double w = weights[i];
    double bw_net = bws_net[i];
    double bw_time = bws_time[i];
    double et = events_time[i];
    arma::vec k;
    // calculating the density value on the network
    if(event_counter[y] > 1){
      // here we have a dupplicated event, we must either use the savec nkde or calculate it and save it latter
      if(map_contains_key(saved_values,y)){
        // we already have the value !
        k = saved_values[y];
      }else{
        // we have not yet the value
        k = esd_kernel_rcpp_arma(kernel_func,edge_mat, neighbour_list ,y,bw_net, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, depth,max_depth);
        saved_values[y] = k;
      }
    }else{
      // this is not a duplicate event, no need to store the value in memory
      k = esd_kernel_rcpp_arma(kernel_func,edge_mat, neighbour_list ,y, bw_net, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, depth,max_depth);
    }

    // ok, now we must calculate the temporal densities
    arma::vec temporal_density = kernel_func(arma::abs(samples_time - et), bw_time);

    //applying the scaling here if required
    if(div == "bw"){
      k = k / bw_net;
      temporal_density = temporal_density / bw_time;
    }

    // and create an awesome spatio-temporal matrix
    arma::mat k_mat(samples_x.n_elem, samples_time.n_elem);

    int mj = temporal_density.n_elem;
    for(int j = 0; j < mj; j++){
      k_mat.col(j) = k * temporal_density[j];
    }


    // getting the actual base_k values (summed at each iteration)
    base_k += (k_mat*w);
    // calculating the new value
    count.fill(0.0);
    arma::umat ids = arma::find(k_mat>0);
    count.elem(ids).fill(1.0);
    base_count += count*w;
  };

  // standardise the result if div = n
  if(div == "n"){
    base_k = base_k / base_count;
  }

  List results = List::create(
    Named("sum_k") = base_k,
    Named("n") = base_count );
  return results;
}

