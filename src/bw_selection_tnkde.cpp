#include "spNetwork.h"
#include "base_kernel_funtions.h"
#include "matrices_functions.h"

#include <cmath>

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// THE FUNCTIONS TO CALCULATE BW SELECTION CV LOO TEMPORAL
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// THE WORKER FUNCTION, CALCULATING FOR EACH BW the kernel value at each event point
// NOTE : this will return a cube : nrow = bws_net, ncol = bws_time, nframe = observations
arma::cube ess_kernel_loo_tnkde(fptros kernel_func, arma::sp_mat& edge_mat,
                             NumericVector& events,
                             NumericVector& time_events,
                             List& neighbour_list ,int v, arma::vec bws_net, arma::vec bws_time,
                             NumericVector& line_weights, int max_depth){

  //step0 : generate the queue
  int depth = 0;
  queue <List> data_holder;
  arma::cube kvalues(bws_net.n_elem, bws_time.n_elem, events.length());
  kvalues.fill(0.0);
  //step1 : generate the first case
  List cas1 = List::create(Named("d")=0.0,
                           Named("v") = v,
                           Named("prev_node") = -999,
                           Named("depth") = 0
  );
  data_holder.push(cas1);

  double max_bw_net = arma::max(bws_net);

  arma::mat k_mat(bws_net.n_elem, bws_time.n_elem);
  double bw_net, bw_time;

  //lancement des iterations
  while(data_holder.empty()==FALSE){
    k_mat.fill(0.0);
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
            for(int i = 0; i < bws_net.n_elem ; i++){
              bw_net = bws_net(i);
              double kernel_net =  kernel_func(d2,bw_net);
              for(int j = 0 ; j < bws_time.n_elem; j ++){
                bw_time = bws_time(j);
                double kernel_time = kernel_func(std::abs(time_events[v]-time_events[index]), bw_time);
                k_mat(i,j) = k_mat(i,j) + (kernel_net * kernel_time);
              }
            }
            // calculating the kernel_net value (here adding iterations on bws_net and bw_time)

            // updating the slice
            //kvalues[index] = kvalues[index] + (kernel_net * kernel_time);
            kvalues.slice(index) += k_mat;
          }
          //evaluating for the next move
          if (d2< max_bw_net and new_depth<max_depth){
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


// THE WORKER FUNCTION, CALCULATING FOR EACH BW the kernel value at each event point
// will return a cube (nrow = bws_net, ncol = bws_time, slices = events.length)
NumericVector get_loo_values_simple(List neighbour_list, NumericVector samples, NumericVector sweights,
                                    NumericVector events, NumericVector events_time,
                                    NumericVector weights,
                                    arma::vec bws_net, arma::vec bws_time, std::string kernel_name,
                                    DataFrame line_list, int max_depth){

  //selecting the kernel function
  fptros kernel_func = select_kernelos(kernel_name);

  //step0 extract the columns of the dataframe
  NumericVector line_weights = line_list["weight"];

  //step 1 : mettre toutes les valeurs a 0 (bw_net 0 et bw_time 0)
  arma::cube base_k(bws_net.n_elem, bws_time.n_elem, events.length());

  double k_net, k_time;
  arma::mat startValues(bws_net.n_elem, bws_time.n_elem);
  for(int i = 0; i < bws_net.n_elem; i ++){
    k_net = kernel_func(0,bws_net[i]);
    for(int j = 0; j < bws_time; j++){
      k_time = kernel_func(0,bws_time[j]);
      startValues(i,j) = k_net * k_time;
    }
  }
  // filling the base cube
  for (int i = 0; i < events.length(); i++){
    base_k.slice(i) = startValues * sweights(i)
  }


  //calculer la matrice des lignes
  //IntegerMatrix edge_mat = make_matrix(line_list,neighbour_list);
  arma::sp_mat edge_mat = make_matrix_sparse(line_list,neighbour_list);

  //step2 : iterer sur chaque event
  int cnt_e = events.length()-1;
  for(int i=0; i <= cnt_e; ++i){
    //preparer les differentes valeurs de departs pour l'event y
    int y = events[i];
    double w = weights[i];
    // launching recursion
    // here we got the the influences of the vertex y on each other vertices
    arma::cube k = ess_kernel_loo_tnkde(kernel_func, edge_mat,
                                           samples, events_time, neighbour_list,
                                           y, bws_net, bws_time,
                                           line_weights, max_depth);
    // getting the actual base_k values (summed at each iteration)
    base_k += (k*w);
  };
  return base_k;
}

