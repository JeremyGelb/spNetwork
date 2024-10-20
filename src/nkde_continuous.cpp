#include "spNetwork.h"
#include "base_kernel_funtions.h"
#include "matrices_functions.h"


//#####################################################################################
// ######################  THE WORKER FUNCTIONS  ######################################
//#####################################################################################


//' @title The worker function to calculate continuous NKDE (with ARMADILLO and sparse matrix)
//' @name continuousWorker_sparse
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param samples_k a numeric vector of the actual kernel values, updates at
//' each recursion
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param v the actual node to consider for the recursion (int)
//' @param bw the kernel bandwidth
//' @param line_weights a vector with the length of the edges
//' @param samples_edgeid a vector associating each sample to an edge
//' @param samples_coords a matrix with the X and Y coordinates of the samples
//' @param nodes_coords a matrix with the X and Y coordinates of the nodes
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a vector with the kernel values calculated for each samples from
//' the first node given
//' @keywords internal
arma::vec esc_kernel_rcpp_arma_sparse(fptr kernel_func, List &neighbour_list, arma::sp_imat &edge_mat, int v, double bw, NumericVector &line_weights, arma::ivec &samples_edgeid, arma::mat &samples_coords, arma::mat &nodes_coords, int max_depth){

  // let me create the vector of kernel values
  arma::vec samples_k(samples_edgeid.n_elem);
  //samples_k.fill(0.0);

  // instead of a recursion, we will iterate over cases stored in a queue
  // a case will be an edge for which we want to recalculate the density
  // of events on its
  struct acase{
    int v1;
    int v2;
    int l;
    int depth;
    double d;
    double alpha;
  };

  //queue<acase> data_holder;
  std::vector<acase> data_holder;

  // let us prepare the first cases
  IntegerVector v_neighbours = neighbour_list[v-1];
  int vl = v_neighbours.length();

  double alpha = 2.0/vl;
  for(int j = 0; j < vl ; j++){
    int v2 = v_neighbours[j];
    int l = edge_mat(v,v2);
    acase el = {v,v2,l,0,0.0,alpha};
    //Rcout << "Adding this cas : d="<<0.0<<", v1="<<v<<", v2="<<v2<<", l="<<l<<"\n";
    data_holder.push_back(el);
  }

  // let us start the iterations
  while(data_holder.empty()==FALSE){
    //acase cas = data_holder.front();
    //data_holder.pop();
    acase cas = data_holder.back();
    data_holder.pop_back();

    int v1 = cas.v1;
    int v2 = cas.v2;
    //Rcout << "Iterating on this cas : d="<<d<<", v1="<<v1<<", v2="<<v2<<", l="<<l<<" alpha="<<alpha<<"\n";


    // calculating the density values of the sample on that edge

    arma::uvec test = arma::find(samples_edgeid == cas.l);

    // we do the calculus only if we have some coords to calculate on
    if(test.n_elem  > 0){
      arma::mat sampling_coords = samples_coords.rows(test);

      // and calculate their distance to the start node
      // NOTE : I could use some memoization here ?

      //arma::colvec x_dists = calcEuclideanDistance3(sampling_coords, nodes_coords.row(cas.v1-1)) + cas.d;

      arma::colvec x_dists =  arma::sqrt(arma::sum(arma::pow(sampling_coords.each_row() - nodes_coords.row(v1-1),2),1)) + cas.d;

      //arma::vec x_dists = arma::sqrt(arma::pow((sampling_x - nodes_x[cas.v1-1]),2) + arma::pow((sampling_y - nodes_y[cas.v1-1]),2)) + cas.d;
      // NumericVector tempx = wrap(x_dists);
      //Rcout << tempx <<"\n";
      // and calculating the new kernel value
      arma::vec new_k = kernel_func(x_dists,bw)*cas.alpha;
      samples_k.elem(test) += new_k;
      //NumericVector tempk = wrap(samples_k);
      //Rcout << tempk <<"\n";
    }




    if(bw >= cas.d){

      // noice, now we can check the reachable nodes and create new cases
      IntegerVector v2_neighbours = neighbour_list[v2-1];
      int n = v2_neighbours.length();

      int new_depth;
      //updating depth
      if(n>2){
        new_depth = cas.depth+1;
      }else{
        new_depth = cas.depth;
      }

      if(n>1){
        // we must continue only if we are not too deep
        if(new_depth <= max_depth){
          // and iterate over the neighbours
          for(int j = 0; j < n ; j++){

            // ---- first, the simple case, we are not going back ----
            int v3 = v2_neighbours[j];
            // int l2 = edge_mat(cas.v2,v3);
            // double d2 = cas.d + line_weights[cas.l-1];
            // first case, we must back fire
            if(v3 == v1){
              if(n>2){
                //acase new_case = {cas.v2,v3,l2,new_depth,d2,(cas.alpha * ((2.0-n)/n))};
                //Rcout << "  adding this cas : d="<<d2<<", v2="<<v2<<", v3="<<v3<<", l2="<<l2<<" alpha="<<new_alpha<<"\n";
                //data_holder.push(new_case);
                acase new_case = {cas.v2,v3,edge_mat(v2,v3),new_depth,(cas.d + line_weights[cas.l-1]),(cas.alpha * ((2.0-n)/n))};
                data_holder.push_back(new_case);
                //data_holder.push_back((struct acase){cas.v2,v3,edge_mat(cas.v2,v3),new_depth,(cas.d + line_weights[cas.l-1]),(cas.alpha * ((2.0-n)/n))});
              }
            }else{
              //acase new_case = {cas.v2,v3,l2,new_depth,d2,(cas.alpha * (2.0/n))};
              //Rcout << "  adding this cas : d="<<d2<<", v2="<<v2<<", v3="<<v3<<", l2="<<l2<<" alpha="<<new_alpha<<"\n";
              //data_holder.push(new_case);
              acase new_case = {cas.v2,v3,edge_mat(cas.v2,v3),new_depth,(cas.d + line_weights[cas.l-1]),(cas.alpha * (2.0/n))};
              data_holder.push_back(new_case);
              //data_holder.push_back((struct acase){cas.v2,v3,edge_mat(cas.v2,v3),new_depth,(cas.d + line_weights[cas.l-1]),(cas.alpha * (2.0/n))});
            }
          }
        }
      }
    }
  }
  return samples_k;
}




// arma::vec esc_kernel_rcpp_arma_sparse_old(fptr kernel_func, List &neighbour_list, arma::sp_mat &edge_mat, int v, double bw, NumericVector &line_weights, arma::ivec &samples_edgeid, arma::vec &samples_x, arma::vec &samples_y, arma::vec &nodes_x, arma::vec &nodes_y, int max_depth){
//
//   // let me create the vector of kernel values
//   arma::vec samples_k(samples_edgeid.n_elem);
//   //samples_k.fill(0.0);
//
//   // instead of a recursion, we will iterate over cases stored in a queue
//   // a case will be an edge for which we want to recalculate the density
//   // of events on its
//   struct acase{
//     int v1;
//     int v2;
//     int l;
//     int depth;
//     double d;
//     double alpha;
//   };
//
//   //queue<acase> data_holder;
//   std::vector<acase> data_holder;
//
//   // let us prepare the first cases
//   IntegerVector v_neighbours = neighbour_list[v-1];
//   double alpha = 2.0/v_neighbours.length();
//   for(int j = 0; j < v_neighbours.length(); j++){
//     int v2 = v_neighbours[j];
//     int l = edge_mat(v,v2);
//     acase el = {v,v2,l,0,0.0,alpha};
//     //Rcout << "Adding this cas : d="<<0.0<<", v1="<<v<<", v2="<<v2<<", l="<<l<<"\n";
//     data_holder.push_back(el);
//   }
//
//   // let us start the iterations
//   while(data_holder.empty()==FALSE){
//     //acase cas = data_holder.front();
//     //data_holder.pop();
//     acase cas = data_holder.back();
//     data_holder.pop_back();
//
//     //Rcout << "Iterating on this cas : d="<<d<<", v1="<<v1<<", v2="<<v2<<", l="<<l<<" alpha="<<alpha<<"\n";
//
//
//     // calculating the density values of the sample on that edge
//
//     arma::uvec test = arma::find(samples_edgeid == cas.l);
//     arma::vec sampling_x = samples_x.elem(test);
//     arma::vec sampling_y = samples_y.elem(test);
//
//     // NOTE : I would like to use the precalcualte distance given here : mat_dist_samples
//
//     // and calculate their distance to the start node
//     // NOTE : I could use some memoization here ?
//     arma::vec x_dists = arma::sqrt(arma::pow((sampling_x - nodes_x[cas.v1-1]),2) + arma::pow((sampling_y - nodes_y[cas.v1-1]),2)) + cas.d;
//     NumericVector tempx = wrap(x_dists);
//     //Rcout << tempx <<"\n";
//     // and calculating the new kernel value
//     arma::vec new_k = kernel_func(x_dists,bw)*cas.alpha;
//     samples_k.elem(test) += new_k;
//     NumericVector tempk = wrap(samples_k);
//     //Rcout << tempk <<"\n";
//
//
//     if(bw >= cas.d){
//
//       // noice, now we can check the reachable nodes and create new cases
//       IntegerVector v2_neighbours = neighbour_list[cas.v2-1];
//       int n = v2_neighbours.length();
//
//       int new_depth;
//       //updating depth
//       if(n>2){
//         new_depth = cas.depth+1;
//       }else{
//         new_depth = cas.depth;
//       }
//
//       if(n>1){
//         // we must continue only if we are not too deep
//         if(new_depth <= max_depth){
//           // and iterate over the neighbours
//           for(int j = 0; j < n ; j++){
//
//             // ---- first, the simple case, we are not going back ----
//             int v3 = v2_neighbours[j];
//             int l2 = edge_mat(cas.v2,v3);
//             double d2 = cas.d + line_weights[cas.l-1];
//             // first case, we must back fire
//             if(v3 == cas.v1){
//               if(n>2){
//                 //acase new_case = {cas.v2,v3,l2,new_depth,d2,(cas.alpha * ((2.0-n)/n))};
//                 //Rcout << "  adding this cas : d="<<d2<<", v2="<<v2<<", v3="<<v3<<", l2="<<l2<<" alpha="<<new_alpha<<"\n";
//                 //data_holder.push(new_case);
//                 //data_holder.push_back(new_case);
//                 data_holder.push_back((struct acase){cas.v2,v3,l2,new_depth,d2,(cas.alpha * ((2.0-n)/n))});
//               }
//             }else{
//               //acase new_case = {cas.v2,v3,l2,new_depth,d2,(cas.alpha * (2.0/n))};
//               //Rcout << "  adding this cas : d="<<d2<<", v2="<<v2<<", v3="<<v3<<", l2="<<l2<<" alpha="<<new_alpha<<"\n";
//               //data_holder.push(new_case);
//               //data_holder.push_back(new_case);
//               data_holder.push_back((struct acase){cas.v2,v3,l2,new_depth,d2,(cas.alpha * (2.0/n))});
//             }
//           }
//         }
//       }
//     }
//   }
//   return samples_k;
// }


// arma::vec esc_kernel_rcpp_arma_sparse(fptr kernel_func, List &neighbour_list, arma::sp_mat &edge_mat, int v, double bw, NumericVector &line_weights, arma::vec &samples_edgeid, arma::vec &samples_x, arma::vec &samples_y, arma::vec &nodes_x, arma::vec &nodes_y, int max_depth){
//
//   // let me create the vector of kernel values
//   arma::vec samples_k(samples_edgeid.n_elem);
//   //samples_k.fill(0.0);
//
//   // instead of a recursion, we will iterate over cases stored in a queue
//   // a case will be an edge for which we want to recalculate the density
//   // of events on its
//   struct acase{
//     int v1;
//     int v2;
//     int l;
//     int depth;
//     double d;
//     double alpha;
//   };
//
//   //queue<acase> data_holder;
//   std::vector<acase> data_holder;
//
//   // let us prepare the first cases
//   IntegerVector v_neighbours = neighbour_list[v-1];
//   double alpha = 2.0/v_neighbours.length();
//   for(int j = 0; j < v_neighbours.length(); j++){
//     int v2 = v_neighbours[j];
//     int l = edge_mat(v,v2);
//     acase el = {v,v2,l,0,0.0,alpha};
//     //Rcout << "Adding this cas : d="<<0.0<<", v1="<<v<<", v2="<<v2<<", l="<<l<<"\n";
//     data_holder.push_back(el);
//   }
//
//   // let us start the iterations
//   while(data_holder.empty()==FALSE){
//     //acase cas = data_holder.front();
//     //data_holder.pop();
//     acase cas = data_holder.back();
//     data_holder.pop_back();
//     int v1 = cas.v1;
//     int v2 = cas.v2;
//     int depth = cas.depth;
//     double d = cas.d;
//     double alpha = cas.alpha;
//     int l = cas.l;
//
//     //Rcout << "Iterating on this cas : d="<<d<<", v1="<<v1<<", v2="<<v2<<", l="<<l<<" alpha="<<alpha<<"\n";
//
//     //extracting the X and Y coordinates of the starting node
//     int v_m1 = v1-1;
//     double node_x = nodes_x[v_m1];
//     double node_y = nodes_y[v_m1];
//
//     // calculating the density values of the sample on that edge
//     arma::uvec test = arma::find(samples_edgeid==l);
//     arma::vec sampling_x = samples_x.elem(test);
//     arma::vec sampling_y = samples_y.elem(test);
//
//     // and calculate their distance to the start node
//     arma::vec x_dists = arma::sqrt(arma::pow((sampling_x - node_x),2) + arma::pow((sampling_y - node_y),2)) + d;
//     NumericVector tempx = wrap(x_dists);
//     //Rcout << tempx <<"\n";
//     // and calculating the new kernel value
//     arma::vec new_k = kernel_func(x_dists,bw)*alpha;
//     samples_k.elem(test) += new_k;
//     NumericVector tempk = wrap(samples_k);
//     //Rcout << tempk <<"\n";
//
//
//     if(bw >= d){
//
//       // noice, now we can check the reachable nodes and create new cases
//       IntegerVector v2_neighbours = neighbour_list[v2-1];
//       int n = v2_neighbours.length();
//
//       int new_depth;
//       //updating depth
//       if(n>2){
//         new_depth = depth+1;
//       }else{
//         new_depth = depth;
//       }
//
//       if(n>1){
//         // we must continue only if we are not too deep
//         if(new_depth <= max_depth){
//           // and iterate over the neighbours
//           for(int j = 0; j < n ; j++){
//
//             // ---- first, the simple case, we are not going back ----
//             int v3 = v2_neighbours[j];
//             int l2 = edge_mat(v2,v3);
//             double d2 = d + line_weights[l-1];
//             // first case, we must back fire
//             if(v3 == v1){
//               if(n>2){
//                 double p2 = (2.0-n)/n;
//                 double new_alpha = alpha * p2;
//                 acase new_case = {v2,v3,l2,new_depth,d2,new_alpha};
//                 //Rcout << "  adding this cas : d="<<d2<<", v2="<<v2<<", v3="<<v3<<", l2="<<l2<<" alpha="<<new_alpha<<"\n";
//                 //data_holder.push(new_case);
//                 data_holder.push_back(new_case);
//               }
//             }else{
//               double new_alpha = alpha * (2.0/n);
//               acase new_case = {v2,v3,l2,new_depth,d2,new_alpha};
//               //Rcout << "  adding this cas : d="<<d2<<", v2="<<v2<<", v3="<<v3<<", l2="<<l2<<" alpha="<<new_alpha<<"\n";
//               //data_holder.push(new_case);
//               data_holder.push_back(new_case);
//             }
//           }
//         }
//       }
//     }
//   }
//   return samples_k;
// }

// THIS IS THE PREVIOUS VERSION USING RECURSION
// arma::vec esc_kernel_rcpp_arma_sparse(fptr kernel_func, arma::vec samples_k, List neighbour_list, arma::sp_mat edge_mat, int v, int v1, int l1, double d,double alpha, double bw, NumericVector line_weights, arma::vec samples_edgeid, arma::vec samples_x, arma::vec samples_y, arma::vec nodes_x, arma::vec nodes_y , int depth, int max_depth){
//
//   //step1 find the index of the right samples
//   arma::uvec test = arma::find(samples_edgeid==l1);
//   arma::vec sampling_x = samples_x.elem(test);
//   arma::vec sampling_y = samples_y.elem(test);
//   //extracting the X and Y coordinates of the starting node
//   int v_m1 = v-1;
//   double node_x = nodes_x[v_m1];
//   double node_y = nodes_y[v_m1];
//   //step2 calculating the distances for each sample
//   arma::vec x_dists = arma::sqrt(arma::pow((sampling_x - node_x),2) + arma::pow((sampling_y - node_y),2)) + d;
//   //step3 calculating the values of the new kernel
//   arma::vec new_k = kernel_func(x_dists,bw)*alpha;
//   //note:here we seems to always have 0 values
//   samples_k.elem(test) += new_k;
//   //mettre a jour d
//   double d2 = line_weights[l1-1] + d;
//   if((bw>d2) && (depth < max_depth)){
//     //on veut trouver toutes les lignes emannant de v (Lv)
//     IntegerVector v_neighbours = neighbour_list[v1-1];
//     int n = v_neighbours.length();
//     int new_depth;
//     //updating depth
//     if(n>2){
//       new_depth = depth+1;
//     }else{
//       new_depth = depth+0;
//     }
//     if(n>1){
//       for(int j=0; j < n; ++j){
//         //int li = Lv[j];
//         int vi = v_neighbours[j];
//         //int li = edge_dict[v1][vi];
//         int li = edge_mat(v1,vi);
//         if(li==l1){
//           //we must backfire only if we have a true intersection
//           if(n>2){
//             double p2 = (2.0-n)/n;
//             double n_alpha = alpha * p2;
//             samples_k = esc_kernel_rcpp_arma_sparse(kernel_func, samples_k, neighbour_list, edge_mat, v1,vi, li, d2, n_alpha, bw, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, new_depth, max_depth);
//           }
//         }else{
//           double n_alpha = alpha * (2.0/n);
//           samples_k = esc_kernel_rcpp_arma_sparse(kernel_func, samples_k, neighbour_list, edge_mat, v1,vi, li, d2, n_alpha, bw, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, new_depth, max_depth);
//         };
//       };
//     };
//   };
//   return samples_k;
// }



//' @title The worker function to calculate continuous NKDE (with ARMADILLO and integer matrix)
//' @name continuousWorker
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param samples_k a numeric vector of the actual kernel values, updates at
//' each recursion
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param v the actual node to consider for the recursion (int)
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
//' @keywords internal
arma::vec esc_kernel_rcpp_arma(fptr kernel_func, List &neighbour_list, IntegerMatrix &edge_mat, int v, double bw, NumericVector &line_weights, arma::ivec &samples_edgeid, arma::vec &samples_x, arma::vec &samples_y, arma::vec &nodes_x, arma::vec &nodes_y, int max_depth){

  // let me create the vector of kernel values
  arma::vec samples_k(samples_edgeid.n_elem);
  //samples_k.fill(0.0);

  // instead of a recursion, we will iterate over cases stored in a queue
  // a case will be an edge for which we want to recalculate the density
  // of events on its
  struct acase{
    int v1;
    int v2;
    int l;
    int depth;
    double d;
    double alpha;
  };

  //queue<acase> data_holder;
  std::vector<acase> data_holder;

  // let us prepare the first cases
  IntegerVector v_neighbours = neighbour_list[v-1];
  double alpha = 2.0/v_neighbours.length();
  for(int j = 0; j < v_neighbours.length(); j++){
    int v2 = v_neighbours[j];
    int l = edge_mat(v,v2);
    acase el = {v,v2,l,0,0.0,alpha};
    //Rcout << "Adding this cas : d="<<0.0<<", v1="<<v<<", v2="<<v2<<", l="<<l<<"\n";
    //data_holder.push(el);
    data_holder.push_back(el);
  }

  // let us start the iterations
  while(data_holder.empty()==FALSE){
    // acase cas = data_holder.front();
    // data_holder.pop();
    acase cas = data_holder.back();
    data_holder.pop_back();
    int v1 = cas.v1;
    int v2 = cas.v2;
    int depth = cas.depth;
    double d = cas.d;
    double alpha = cas.alpha;
    int l = cas.l;

    //Rcout << "Iterating on this cas : d="<<d<<", v1="<<v1<<", v2="<<v2<<", l="<<l<<" alpha="<<alpha<<"\n";

    //extracting the X and Y coordinates of the starting node
    int v_m1 = v1-1;
    double node_x = nodes_x[v_m1];
    double node_y = nodes_y[v_m1];

    // calculating the density values of the sample on that edge
    arma::uvec test = arma::find(samples_edgeid==l);
    arma::vec sampling_x = samples_x.elem(test);
    arma::vec sampling_y = samples_y.elem(test);

    // and calculate their distance to the start node
    arma::vec x_dists = arma::sqrt(arma::pow((sampling_x - node_x),2) + arma::pow((sampling_y - node_y),2)) + d;
    //Rcout << tempx <<"\n";
    // and calculating the new kernel value
    arma::vec new_k = kernel_func(x_dists,bw)*alpha;
    samples_k.elem(test) += new_k;
    //Rcout << tempk <<"\n";


    if(bw > d){

      // noice, now we can check the reachable nodes and create new cases
      IntegerVector v2_neighbours = neighbour_list[v2-1];
      int n = v2_neighbours.length();

      int new_depth;
      //updating depth
      if(n>2){
        new_depth = depth+1;
      }else{
        new_depth = depth;
      }

      if(n>1){
        // we must continue only if we are not too deep
        if(new_depth <= max_depth){
          // and iterate over the neighbours
          for(int j = 0; j < n ; j++){

            // ---- first, the simple case, we are not going back ----
            int v3 = v2_neighbours[j];
            int l2 = edge_mat(v2,v3);
            double d2 = d + line_weights[l-1];
            // first case, we must back fire
            if(v3 == v1){
              if(n>2){
                double new_alpha = alpha * ((2.0-n)/n);
                acase new_case = {v2,v3,l2,new_depth,d2,new_alpha};
                //Rcout << "  adding this cas : d="<<d2<<", v2="<<v2<<", v3="<<v3<<", l2="<<l2<<" alpha="<<new_alpha<<"\n";
                //data_holder.push(new_case);
                data_holder.push_back(new_case);
              }
            }else{
              double new_alpha = alpha * (2.0/n);
              acase new_case = {v2,v3,l2,new_depth,d2,new_alpha};
              //Rcout << "  adding this cas : d="<<d2<<", v2="<<v2<<", v3="<<v3<<", l2="<<l2<<" alpha="<<new_alpha<<"\n";
              //data_holder.push(new_case);
              data_holder.push_back(new_case);
            }
          }
        }
      }
    }
  }
  return samples_k;
}


// THIS IS THE PREVIOUS VERSION USING RECURSION
// arma::vec esc_kernel_rcpp_arma(fptr kernel_func, arma::vec samples_k, List neighbour_list, IntegerMatrix edge_mat, int v, int v1, int l1, double d,double alpha, double bw, NumericVector line_weights, arma::vec samples_edgeid, arma::vec samples_x, arma::vec samples_y, arma::vec nodes_x, arma::vec nodes_y , int depth, int max_depth){
//   //step1 find the index of the right samples
//   arma::uvec test = arma::find(samples_edgeid==l1);
//   arma::vec sampling_x = samples_x.elem(test);
//   arma::vec sampling_y = samples_y.elem(test);
//   //extracting the X and Y coordinates of the starting node
//   int v_m1 = v-1;
//   double node_x = nodes_x[v_m1];
//   double node_y = nodes_y[v_m1];
//   //step2 calculating the distances for each sample
//   arma::vec x_dists = arma::sqrt(arma::pow((sampling_x - node_x),2) + arma::pow((sampling_y - node_y),2)) + d;
//   //step3 calculating the values of the new kernel
//   arma::vec new_k = kernel_func(x_dists,bw)*alpha;
//   //note:here we seems to always have 0 values
//   samples_k.elem(test) += new_k;
//   //mettre a jour d
//   double d2 = line_weights[l1-1] + d;
//   if((bw>d2) && (depth < max_depth)){
//     //on veut trouver toutes les lignes emannant de v (Lv)
//     IntegerVector v_neighbours = neighbour_list[v1-1];
//     int n = v_neighbours.length();
//     int new_depth;
//     //updating depth
//     if(n>2){
//       new_depth = depth+1;
//     }else{
//       new_depth = depth+0;
//     }
//     if(n>1){
//       for(int j=0; j < n; ++j){
//         //int li = Lv[j];
//         int vi = v_neighbours[j];
//         //int li = edge_dict[v1][vi];
//         int li = edge_mat(v1,vi);
//         if(li==l1){
//           //we must backfire only if we have a true intersection
//           if(n>2){
//             double p2 = (2.0-n)/n;
//             double n_alpha = alpha * p2;
//             samples_k = esc_kernel_rcpp_arma(kernel_func, samples_k, neighbour_list, edge_mat, v1,vi, li, d2, n_alpha, bw, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, new_depth, max_depth);
//           }
//         }else{
//           double n_alpha = alpha * (2.0/n);
//           samples_k = esc_kernel_rcpp_arma(kernel_func, samples_k, neighbour_list, edge_mat, v1,vi, li, d2, n_alpha, bw, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, new_depth, max_depth);
//         };
//       };
//     };
//   };
//   return samples_k;
// }




//#####################################################################################
// ######################  THE CALLED FUNCTIONS  ######################################
//#####################################################################################


//' @title The main function to calculate continuous NKDE (with ARMADILO and sparse matrix)
//' @name continuousfunction2
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
//' @param div The divisor to use for the kernel. Must be "n" (the number of events within the radius around each sampling point), "bw" (the bandwidth) "none" (the simple sum).
//' @return a DataFrame with two columns : the kernel values (sum_k) and the number of events for each sample (n)
//' @export
//' @keywords internal
//'
// [[Rcpp::export]]
 DataFrame continuous_nkde_cpp_arma_sparse(List neighbour_list, NumericVector events,
                                           NumericVector weights, DataFrame samples,
                                           NumericVector bws, std::string kernel_name,
                                           DataFrame nodes, DataFrame line_list,
                                           int max_depth, bool verbose, std::string div = "bw"){

   //continuous_nkde_cpp_arma_sparse(neighbour_list,events$vertex_id, events$weight, samples@data, bws, kernel_name, nodes@data, graph_result$linelist, max_depth, tol, verbose)
   //selecting the kernel function
   fptr kernel_func = select_kernel(kernel_name);

   //step0 extract the columns of the dataframe
   NumericVector line_weights = line_list["weight"];
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
   arma::vec samples_k(samples.nrow());
   arma::vec count(samples.nrow());
   //base_k.fill(0.0);
   //base_count.fill(0.0);
   //count.fill(0.0);
   //samples_k.fill(0.0);

   //calculer le dictionnaire des lignes
   //IntegerMatrix edge_mat = make_matrix(line_list,neighbour_list);
   arma::sp_imat edge_mat = make_imatrix_sparse(line_list, neighbour_list);
   //step2 : iterer sur chaque event
   int cnt_e = events.length()-1;
   Progress p(cnt_e, verbose);
   for(int i=0; i <= cnt_e; ++i){
     p.increment(); // update progress
     //preparer les differentes valeurs de departs pour l'event y
     int y = events[i];
     double w = weights[i];
     double bw = bws[i];

     // on peut maintenant calculer la densite emanant de l'event y
     samples_k = esc_kernel_rcpp_arma_sparse(kernel_func, neighbour_list, edge_mat, y, bw, line_weights, samples_edgeid, samples_coords,nodes_coords, max_depth);
     if(div == "bw"){
       base_k += (samples_k*w) / bw;
     }else{
       base_k += samples_k*w;
     }
     // calculating the new value
     count.fill(0.0);
     arma::uvec ids = arma::find(samples_k>0);
     count.elem(ids).fill(1.0);
     base_count += count*w;

   };
   NumericVector v1 = Rcpp::NumericVector(base_k.begin(), base_k.end());
   NumericVector v2 = Rcpp::NumericVector(base_count.begin(), base_count.end());
   DataFrame df =  DataFrame::create( Named("sum_k") = v1 ,Named("n") = v2);
   return df;
 }




// DataFrame continuous_nkde_cpp_arma_sparse_old(List neighbour_list, NumericVector events,
//                                           NumericVector weights, DataFrame samples,
//                                           NumericVector bws, std::string kernel_name,
//                                           DataFrame nodes, DataFrame line_list,
//                                           int max_depth, bool verbose, std::string div = "bw"){
//
//   //continuous_nkde_cpp_arma_sparse(neighbour_list,events$vertex_id, events$weight, samples@data, bws, kernel_name, nodes@data, graph_result$linelist, max_depth, tol, verbose)
//   //selecting the kernel function
//   fptr kernel_func = select_kernel(kernel_name);
//
//   //step0 extract the columns of the dataframe
//   NumericVector line_weights = line_list["weight"];
//   arma::ivec samples_edgeid = as<arma::ivec>(samples["edge_id"]);
//   arma::vec samples_x =  as<arma::vec>(samples["X_coords"]);
//   arma::vec samples_y =  as<arma::vec>(samples["Y_coords"]);
//   arma::vec nodes_x = as<arma::vec>(nodes["X_coords"]);
//   arma::vec nodes_y = as<arma::vec>(nodes["Y_coords"]);
//
//   //step 1 : mettre toutes les valeurs a 0
//   arma::vec base_k(samples.nrow());
//   arma::vec base_count(samples.nrow());
//   arma::vec samples_k(samples.nrow());
//   arma::vec count(samples.nrow());
//   //base_k.fill(0.0);
//   //base_count.fill(0.0);
//   //count.fill(0.0);
//   //samples_k.fill(0.0);
//
//   //calculer le dictionnaire des lignes
//   //IntegerMatrix edge_mat = make_matrix(line_list,neighbour_list);
//   arma::sp_mat edge_mat = make_matrix_sparse(line_list, neighbour_list);
//   //step2 : iterer sur chaque event
//   int cnt_e = events.length()-1;
//   Progress p(cnt_e, verbose);
//   for(int i=0; i <= cnt_e; ++i){
//     p.increment(); // update progress
//     //preparer les differentes valeurs de departs pour l'event y
//     int y = events[i];
//     double w = weights[i];
//     double bw = bws[i];
//
//     // on peut maintenant calculer la densite emanant de l'event y
//     samples_k = esc_kernel_rcpp_arma_sparse(kernel_func, neighbour_list, edge_mat, y, bw, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, max_depth);
//     if(div == "bw"){
//       base_k += (samples_k*w) / bw;
//     }else{
//       base_k += samples_k*w;
//     }
//     // calculating the new value
//     count.fill(0.0);
//     arma::uvec ids = arma::find(samples_k>0);
//     count.elem(ids).fill(1.0);
//     base_count += count*w;
//
//   };
//   NumericVector v1 = Rcpp::NumericVector(base_k.begin(), base_k.end());
//   NumericVector v2 = Rcpp::NumericVector(base_count.begin(), base_count.end());
//   DataFrame df =  DataFrame::create( Named("sum_k") = v1 ,Named("n") = v2);
//   return df;
// }


// THIS IS THE PREVIOUS VERSION USING RECURSION
// DataFrame continuous_nkde_cpp_arma_sparse(List neighbour_list, NumericVector events, NumericVector weights, DataFrame samples, NumericVector bws, std::string kernel_name, DataFrame nodes, DataFrame line_list, int max_depth, bool verbose){
//
//   //continuous_nkde_cpp_arma_sparse(neighbour_list,events$vertex_id, events$weight, samples@data, bws, kernel_name, nodes@data, graph_result$linelist, max_depth, tol, verbose)
//   //selecting the kernel function
//   fptr kernel_func = select_kernel(kernel_name);
//
//   //step0 extract the columns of the dataframe
//   NumericVector line_weights = line_list["weight"];
//   arma::vec samples_edgeid = as<arma::vec>(samples["edge_id"]);
//   arma::vec samples_x =  as<arma::vec>(samples["X_coords"]);
//   arma::vec samples_y =  as<arma::vec>(samples["Y_coords"]);
//   arma::vec nodes_x = as<arma::vec>(nodes["X_coords"]);
//   arma::vec nodes_y = as<arma::vec>(nodes["Y_coords"]);
//
//   //step 1 : mettre toutes les valeurs a 0
//   arma::vec base_k(samples.nrow());
//   arma::vec base_count(samples.nrow());
//   arma::vec samples_k(samples.nrow());
//   arma::vec count(samples.nrow());
//   base_k.fill(0.0);
//   base_count.fill(0.0);
//   count.fill(0.0);
//   samples_k.fill(0.0);
//
//   //calculer le dictionnaire des lignes
//   //IntegerMatrix edge_mat = make_matrix(line_list,neighbour_list);
//   arma::sp_mat edge_mat = make_matrix_sparse(line_list, neighbour_list);
//   //step2 : iterer sur chaque event
//   int cnt_e = events.length()-1;
//   Progress p(cnt_e, verbose);
//   for(int i=0; i <= cnt_e; ++i){
//     p.increment(); // update progress
//     //preparer les differentes valeurs de departs pour l'event y
//     int y = events[i];
//     double w = weights[i];
//     double bw = bws[i];
//     //on veut trouver toutes les voisins emannant de y
//     IntegerVector y_neighbours = neighbour_list[y-1];
//     int cnt_y = y_neighbours.length();
//     double alpha = 2.0/cnt_y;
//     for(int j=0; j < cnt_y; ++j){
//       //preparing all values for this loop
//       int vi = y_neighbours[j];
//       //int li = edge_dict[y][vi];
//       int li = edge_mat(y,vi);
//       double d = 0.0 ;
//       int depth = 0 ;
//       // launching recursion
//       samples_k.fill(0.0);
//       samples_k = esc_kernel_rcpp_arma_sparse(kernel_func, samples_k, neighbour_list, edge_mat ,y,vi,li,d,alpha,bw, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, depth,max_depth);
//       base_k += samples_k*w;
//       // calculating the new value
//       count.fill(0.0);
//       arma::uvec ids = arma::find(samples_k>0);
//       count.elem(ids).fill(1.0);
//       base_count += count*w;
//     };
//
//   };
//   NumericVector v1 = Rcpp::NumericVector(base_k.begin(), base_k.end());
//   NumericVector v2 = Rcpp::NumericVector(base_count.begin(), base_count.end());
//   DataFrame df =  DataFrame::create( Named("sum_k") = v1 ,Named("n") = v2);
//   return df;
// }


//' @title The main function to calculate continuous NKDE (with ARMADILO and integer matrix)
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
//' @param div The divisor to use for the kernel. Must be "n" (the number of events within the radius around each sampling point), "bw" (the bandwidth) "none" (the simple sum).
//' @return a DataFrame with two columns : the kernel values (sum_k) and the number of events for each sample (n)
//' @export
//' @keywords internal
//'
// [[Rcpp::export]]
DataFrame continuous_nkde_cpp_arma(List neighbour_list, NumericVector events,
                                   NumericVector weights, DataFrame samples,
                                   NumericVector bws, std::string kernel_name,
                                   DataFrame nodes, DataFrame line_list,
                                   int max_depth, bool verbose, std::string div = "bw"){

  //continuous_nkde_cpp_arma_sparse(neighbour_list,events$vertex_id, events$weight, samples@data, bws, kernel_name, nodes@data, graph_result$linelist, max_depth, tol, verbose)
  //selecting the kernel function
  fptr kernel_func = select_kernel(kernel_name);

  //step0 extract the columns of the dataframe
  NumericVector line_weights = line_list["weight"];
  arma::ivec samples_edgeid = as<arma::ivec>(samples["edge_id"]);
  arma::vec samples_x =  as<arma::vec>(samples["X_coords"]);
  arma::vec samples_y =  as<arma::vec>(samples["Y_coords"]);
  arma::vec nodes_x = as<arma::vec>(nodes["X_coords"]);
  arma::vec nodes_y = as<arma::vec>(nodes["Y_coords"]);

  //step 1 : mettre toutes les valeurs a 0
  arma::vec base_k(samples.nrow());
  arma::vec base_count(samples.nrow());
  arma::vec samples_k(samples.nrow());
  arma::vec count(samples.nrow());
  // base_k.fill(0.0);
  // base_count.fill(0.0);
  // count.fill(0.0);
  // samples_k.fill(0.0);

  //calculer le dictionnaire des lignes
  //IntegerMatrix edge_mat = make_matrix(line_list,neighbour_list);
  IntegerMatrix edge_mat = make_matrix(line_list, neighbour_list);
  //step2 : iterer sur chaque event
  int cnt_e = events.length()-1;
  Progress p(cnt_e, verbose);
  for(int i=0; i <= cnt_e; ++i){
    p.increment(); // update progress
    //preparer les differentes valeurs de departs pour l'event y
    int y = events[i];
    double w = weights[i];
    double bw = bws[i];

    // on peut maintenant calculer la densite emanant de l'event y
    samples_k = esc_kernel_rcpp_arma(kernel_func, neighbour_list, edge_mat ,y,bw, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, max_depth);
    if(div == "bw"){
      base_k += (samples_k*w) / bw;
    }else{
      base_k += samples_k*w;
    }
    // calculating the new value
    count.fill(0.0);
    arma::uvec ids = arma::find(samples_k>0);
    count.elem(ids).fill(1.0);
    base_count += count*w;

  };
  NumericVector v1 = Rcpp::NumericVector(base_k.begin(), base_k.end());
  NumericVector v2 = Rcpp::NumericVector(base_count.begin(), base_count.end());
  DataFrame df =  DataFrame::create( Named("sum_k") = v1 ,Named("n") = v2);
  return df;
}


// THIS IS THE PREVIOUS VERSION USING RECURSION
// DataFrame continuous_nkde_cpp_arma(List neighbour_list, NumericVector events, NumericVector weights, DataFrame samples, NumericVector bws, std::string kernel_name, DataFrame nodes, DataFrame line_list, int max_depth, bool verbose){
//
//   //selecting the kernel function
//   fptr kernel_func = select_kernel(kernel_name);
//
//   //step0 extract the columns of the dataframe
//   NumericVector line_weights = line_list["weight"];
//   arma::vec samples_edgeid = as<arma::vec>(samples["edge_id"]);
//   arma::vec samples_x =  as<arma::vec>(samples["X_coords"]);
//   arma::vec samples_y =  as<arma::vec>(samples["Y_coords"]);
//   arma::vec nodes_x = as<arma::vec>(nodes["X_coords"]);
//   arma::vec nodes_y = as<arma::vec>(nodes["Y_coords"]);
//
//   //step 1 : mettre toutes les valeurs a 0
//   arma::vec base_k(samples.nrow());
//   arma::vec base_count(samples.nrow());
//   arma::vec samples_k(samples.nrow());
//   arma::vec count(samples.nrow());
//   base_k.fill(0.0);
//   base_count.fill(0.0);
//   count.fill(0.0);
//   samples_k.fill(0.0);
//
//   //calculer le dictionnaire des lignes
//   IntegerMatrix edge_mat = make_matrix(line_list,neighbour_list);
//   //step2 : iterer sur chaque event
//   int cnt_e = events.length()-1;
//   Progress p(cnt_e, verbose);
//   for(int i=0; i <= cnt_e; ++i){
//     p.increment(); // update progress
//     //preparer les differentes valeurs de departs pour l'event y
//     int y = events[i];
//     double w = weights[i];
//     double bw = bws[i];
//     //on veut trouver toutes les voisins emannant de y
//     IntegerVector y_neighbours = neighbour_list[y-1];
//     int cnt_y = y_neighbours.length();
//     double alpha = 2.0/cnt_y;
//     for(int j=0; j < cnt_y; ++j){
//       //preparing all values for this loop
//       int vi = y_neighbours[j];
//       //int li = edge_dict[y][vi];
//       int li = edge_mat(y,vi);
//       double d = 0.0 ;
//       int depth = 0 ;
//       // launching recursion
//       samples_k.fill(0.0);
//       samples_k = esc_kernel_rcpp_arma(kernel_func, samples_k, neighbour_list, edge_mat ,y,vi,li,d,alpha,bw, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, depth,max_depth);
//       base_k += samples_k*w;
//       // calculating the new value
//       count.fill(0.0);
//       arma::uvec ids = arma::find(samples_k>0);
//       count.elem(ids).fill(1.0);
//       base_count += count*w;
//     };
//
//   };
//   NumericVector v1 = Rcpp::NumericVector(base_k.begin(), base_k.end());
//   NumericVector v2 = Rcpp::NumericVector(base_count.begin(), base_count.end());
//   DataFrame df =  DataFrame::create( Named("sum_k") = v1 ,Named("n") = v2);
//   return df;
// }


//#####################################################################################
// #################  THE EXECUTION FUNCTIONS FOR TNKDE  ##############################
//#####################################################################################

//' @title The main function to calculate continuous TNKDE (with ARMADILO and sparse matrix)
//' @name tnkdecontinuousfunction
//' @param neighbour_list a list of the neighbours of each node
//' @param events a numeric vector of the node id of each event
//' @param events_time a numeric vector with the time for the events
//' @param weights a numeric vector of the weight of each event
//' @param samples a DataFrame of the samples (with spatial coordinates and belonging edge)
//' @param samples_time a NumericVector indicating when to do the samples
//' @param bws_net the network kernel bandwidths for each event
//' @param bws_time the time kernel bandwidths for each event
//' @param kernel_name the name of the kernel to use
//' @param nodes a DataFrame representing the nodes of the graph (with spatial coordinates)
//' @param line_list a DataFrame representing the lines of the graph
//' @param max_depth the maximum recursion depth (after which recursion is stopped)
//' @param verbose a boolean indicating if the function must print its progress
//' @param div a string indicating how to standardize the kernel values
//' @return a List with two matrices: the kernel values (sum_k) and the number of events for each sample (n)
//' @export
//' @keywords internal
//'
// [[Rcpp::export]]
List continuous_tnkde_cpp_arma_sparse(List neighbour_list,
                                     IntegerVector events, NumericVector events_time,NumericVector weights,
                                     DataFrame samples, arma::vec samples_time,
                                     NumericVector bws_net,
                                     NumericVector bws_time,
                                     std::string kernel_name, DataFrame nodes, DataFrame line_list,
                                     int max_depth, bool verbose, std::string div){

  //selecting the kernel function
  fptr kernel_func = select_kernel(kernel_name);

  //step0 extract the columns of the dataframe
  NumericVector line_weights = line_list["weight"];
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
  arma::vec samples_k(samples.nrows());

  // NOTE It is possible that we must recalculate the density for the same vertex
  // we will store them instead
  std::map<int, int> event_counter = count_values_intvec(events);
  std::map<int, arma::vec> saved_values;

  //calculer le dictionnaire des lignes
  //IntegerMatrix edge_mat = make_matrix(line_list,neighbour_list);
  arma::sp_imat edge_mat = make_imatrix_sparse(line_list, neighbour_list);

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

    // calculating the density value on the network
    if(event_counter[y] > 1){
      // here we have a dupplicated event, we must either use the savec nkde or calculate it and save it latter
      if(map_contains_key(saved_values,y)){
        // we already have the value !
        samples_k = saved_values[y];
      }else{
        // we have not yet the value
        samples_k = esc_kernel_rcpp_arma_sparse(kernel_func, neighbour_list, edge_mat, y, bw_net, line_weights, samples_edgeid, samples_coords, nodes_coords, max_depth);
        saved_values[y] = samples_k;
      }
    }else{
      // this is not a duplicate event, no need to store the value in memory
      samples_k = esc_kernel_rcpp_arma_sparse(kernel_func, neighbour_list, edge_mat, y, bw_net, line_weights, samples_edgeid, samples_coords, nodes_coords, max_depth);
    }

    // ok, now we must calculate the temporal densities
    arma::vec temporal_density = kernel_func(arma::abs(samples_time - et), bw_time);

    // and standardize by bw if required
    if(div == "bw"){
      temporal_density = temporal_density / bw_time;
      samples_k = samples_k / bw_net;
    }

    // and create an awesome spatio-temporal matrix
    arma::mat k_mat(samples.nrows(), samples_time.n_elem);

    int mj = temporal_density.n_elem;
    for(int j = 0; j < mj; j++){
      k_mat.col(j) = samples_k * temporal_density[j];
    }

    base_k += k_mat*w;
    // calculating the new value
    count.fill(0.0);
    arma::umat ids = arma::find(samples_k>0);
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



//' @title The main function to calculate continuous TNKDE (with ARMADILO and integer matrix)
//' @name tnkdecontinuousfunction
//' @param neighbour_list a list of the neighbours of each node
//' @param events a numeric vector of the node id of each event
//' @param events_time a numeric vector with the time for the events
//' @param weights a numeric vector of the weight of each event
//' @param samples a DataFrame of the samples (with spatial coordinates and belonging edge)
//' @param samples_time a NumericVector indicating when to do the samples
//' @param bws_net the network kernel bandwidths for each event
//' @param kernel_name the name of the kernel to use
//' @param nodes a DataFrame representing the nodes of the graph (with spatial coordinates)
//' @param line_list a DataFrame representing the lines of the graph
//' @param max_depth the maximum recursion depth (after which recursion is stopped)
//' @param verbose a boolean indicating if the function must print its progress
//' @param div a string indicating how to standardize the kernel values
//' @return a List with two matrices: the kernel values (sum_k) and the number of events for each sample (n)
//' @export
//' @keywords internal
//'
// [[Rcpp::export]]
List continuous_tnkde_cpp_arma(List neighbour_list,
                              IntegerVector events, NumericVector events_time,NumericVector weights,
                              DataFrame samples, arma::vec samples_time,
                              NumericVector bws_net,
                              NumericVector bws_time,
                              std::string kernel_name, DataFrame nodes, DataFrame line_list,
                              int max_depth, bool verbose, std::string div){

  //selecting the kernel function
  fptr kernel_func = select_kernel(kernel_name);

  //step0 extract the columns of the dataframe
  NumericVector line_weights = line_list["weight"];
  arma::ivec samples_edgeid = as<arma::ivec>(samples["edge_id"]);
  arma::vec samples_x =  as<arma::vec>(samples["X_coords"]);
  arma::vec samples_y =  as<arma::vec>(samples["Y_coords"]);
  arma::vec nodes_x = as<arma::vec>(nodes["X_coords"]);
  arma::vec nodes_y = as<arma::vec>(nodes["Y_coords"]);

  //step 1 : set all values to 0
  // NOTE : we will produce matrices here (tnkde)
  arma::mat base_k(samples_x.n_elem, samples_time.n_elem);
  arma::mat base_count(samples_x.n_elem, samples_time.n_elem);
  arma::mat count(samples_x.n_elem, samples_time.n_elem);
  arma::vec samples_k(samples_x.n_elem);

  // NOTE It is possible that we must recalculate the density for the same vertex
  // we will store them instead
  std::map<int, int> event_counter = count_values_intvec(events);
  std::map<int, arma::vec> saved_values;

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
    double bw_net = bws_net[i];
    double bw_time = bws_time[i];
    double et = events_time[i];

    // calculating the density value on the network
    if(event_counter[y] > 1){
      // here we have a dupplicated event, we must either use the savec nkde or calculate it and save it latter
      if(map_contains_key(saved_values,y)){
        // we already have the value !
        samples_k = saved_values[y];
      }else{
        // we have not yet the value
        samples_k = esc_kernel_rcpp_arma(kernel_func, neighbour_list, edge_mat, y, bw_net, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, max_depth);
        saved_values[y] = samples_k;
      }
    }else{
      // this is not a duplicate event, no need to store the value in memory
      samples_k = esc_kernel_rcpp_arma(kernel_func, neighbour_list, edge_mat, y, bw_net, line_weights, samples_edgeid, samples_x, samples_y, nodes_x, nodes_y, max_depth);
    }

    // ok, now we must calculate the temporal densities
    arma::vec temporal_density = kernel_func(arma::abs(samples_time - et), bw_time);

    // and standardize by bw if required
    if(div == "bw"){
      temporal_density = temporal_density / bw_time;
      samples_k = samples_k / bw_net;
    }

    // and create an awesome spatio-temporal matrix
    arma::mat k_mat(samples_x.n_elem, samples_time.n_elem);

    int mj = temporal_density.n_elem;

    for(int j = 0; j  <mj; j++){
      k_mat.col(j) = samples_k * temporal_density[j];
    }

    base_k += k_mat*w;
    // calculating the new value
    count.fill(0.0);
    arma::umat ids = arma::find(samples_k>0);
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


