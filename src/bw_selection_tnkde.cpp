#include "spNetwork.h"
#include "base_kernel_funtions.h"
#include "matrices_functions.h"

#include <cmath>


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// THE FUNCTIONS TO CALCULATE BW SELECTION CV LOO TEMPORAL FOR SIMPLE NKDE
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//' @title The worker function to calculate simple TNKDE likelihood cv
//' @name ess_kernel_loo_tnkde
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param events a NumericVector indicating the nodes in graph beeing events
//' @param time_events a NumericVector indicating the timestamp of each event
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param v the actual node to consider (int)
//' @param v_time the time of v (double)
//' @param bws_net an arma::vec with the network bandwidths to consider
//' @param bws_time an arma::vec with the time bandwidths to consider
//' @param line_weights a vector with the length of the edges
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a cube with the impact of the event v on each other event for
//' each pair of bandwidths (cube(bws_net, bws_time, events))
//' @keywords internal
arma::cube ess_kernel_loo_tnkde(fptros kernel_func, arma::sp_mat& edge_mat,
                                NumericVector& events,
                                NumericVector& events_wid,
                                NumericVector& time_events,
                                List& neighbour_list,
                                int v, int wid, double v_time,
                                arma::vec bws_net, arma::vec bws_time,
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
  double bw_net, bw_time;

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
        // comment retrouver le potentiel wid du noeud ?
        //on ne veut pas revenir en arriere !
        if(v2!=prev_node){
          // find the edge between the two nodes
          int edge_id = edge_mat(v,v2);
          double d2 = line_weights[edge_id-1] + d;
          std::vector<int> index = get_all_indeces(events,v2);
          if(index.size() >0 ){
            // il semble que v2 soit un noeud pour lequel au moins un evenement est present
            for(int ii = 0; ii < bws_net.n_elem ; ii++){
              bw_net = bws_net(ii);
              double kernel_net =  kernel_func(d2,bw_net);
              for(int j = 0 ; j < bws_time.n_elem; j ++){
                bw_time = bws_time(j);
                // NOTE : we are doing the bw scaling here
                for (int zz : index){
                  double kernel_time = kernel_func(std::abs(v_time-time_events[zz]), bw_time);
                  kvalues(ii,j,zz) = kvalues(ii,j,zz) + ((kernel_net * kernel_time) * (1.0/(bw_net*bw_time)));
                }

              }
            }
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

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// THE FUNCTION TO CALCULATE BW SELECTION CV LOO TEMPORAL FOR DISCONTINUOUS NKDE
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title The worker function to calculate discontinuous TNKDE likelihood cv
//' @name esd_kernel_loo_tnkde
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param events a NumericVector indicating the nodes in graph beeing events
//' @param time_events a NumericVector indicating the timestamp of each event
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param v the actual node to consider (int)
//' @param v_time the time of v (double)
//' @param bws_net an arma::vec with the network bandwidths to consider
//' @param bws_time an arma::vec with the time bandwidths to consider
//' @param line_weights a vector with the length of the edges
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a cube with the impact of the event v on each other event for
//' each pair of bandwidths (cube(bws_net, bws_time, events))
arma::cube esd_kernel_loo_tnkde(fptros kernel_func, arma::sp_mat& edge_mat,
                                 NumericVector& events,
                                 NumericVector& events_wid,
                                 NumericVector& time_events,
                                 List& neighbour_list,
                                 int v, int wid, double v_time,
                                 arma::vec bws_net, arma::vec bws_time,
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
                           Named("depth") = 0,
                           Named("alpha") = 1
  );
  data_holder.push(cas1);

  double max_bw_net = arma::max(bws_net);

  double bw_net, bw_time;

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

    //step 3 prepare the new alpha value
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

          //est ce que v2 est un evenement pour lequel on doit garder la valeur

          std::vector<int> index = get_all_indeces(events,v2);
          if(index.size() >0 ){

            for(int ii = 0; ii < bws_net.n_elem ; ii++){
              bw_net = bws_net(ii);
              double kernel_net =  kernel_func(d2,bw_net) * new_alpha;
              for(int j = 0 ; j < bws_time.n_elem; j ++){
                bw_time = bws_time(j);
                for (int zz : index){
                  double kernel_time = kernel_func(std::abs(v_time-time_events[zz]), bw_time);
                  kvalues(ii,j,zz) = kvalues(ii,j,zz) + ((kernel_net * kernel_time) * (1.0/(bw_net*bw_time)));
                }
              }
            }
          }

          //evaluating for the next move
          if (d2< max_bw_net and new_depth<max_depth){
            List new_cas = List::create(Named("d")=d2,
                                        Named("v") = v2,
                                        Named("prev_node") = v,
                                        Named("depth") = new_depth,
                                        Named("alpha")=new_alpha
            );
            data_holder.push(new_cas);
          }
        }
      }
    }
  }

  return kvalues;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// THE FUNCTION TO CALCULATE BW SELECTION CV LOO TEMPORAL FOR CONTINUOUS NKDE
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title The worker function to calculate continuous TNKDE likelihood cv
//' @name esc_kernel_loo_tnkde
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param events a NumericVector indicating the nodes in graph beeing events
//' @param time_events a NumericVector indicating the timestamp of each event
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param v the actual node to consider (int)
//' @param v_time the time of v (double)
//' @param bws_net an arma::vec with the network bandwidths to consider
//' @param bws_time an arma::vec with the time bandwidths to consider
//' @param line_weights a vector with the length of the edges
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a cube with the impact of the event v on each other event for
//' each pair of bandwidths (cube(bws_net, bws_time, events))
arma::cube esc_kernel_loo_tnkde(fptros kernel_func, arma::sp_mat& edge_mat,
                                 NumericVector& events,
                                 NumericVector& events_wid,
                                 NumericVector& time_events,
                                 List& neighbour_list,
                                 int v, int wid, double v_time,
                                 arma::vec bws_net, arma::vec bws_time,
                                 NumericVector& line_weights, int max_depth){
  //step0 : generate the queue
  int depth = 0;
  struct acase{
    int v;
    int prev_node;
    int depth;
    double d;
    double alpha;
  };

  // prepare the data cube and matrix
  arma::cube kvalues(bws_net.n_elem, bws_time.n_elem, events.length());
  kvalues.fill(0.0);
  double max_bw_net = arma::max(bws_net);
  double bw_net, bw_time;

  //queue<acase> data_holder;
  std::vector<acase> data_holder;

  // let us prepare the first cases
  IntegerVector v_neighbours = neighbour_list[v-1];
  double alpha = 2.0/v_neighbours.length();

  acase el = {v,-999,0,0.0,alpha};
  data_holder.push_back(el);


  //lancement des iterations
  while(data_holder.empty()==FALSE){

    //unpacking (imagine some loop unrolling here with a function to deal with.)
    acase cas = data_holder.back();
    data_holder.pop_back();
    int v = cas.v;
    int prev_node = cas.prev_node;
    int depth = cas.depth;
    double d = cas.d;
    double alpha = cas.alpha;


    // we will update the densities on v
    // but only if v is a vertex on wich I can find an event

    std::vector<int> index = get_all_indeces(events,v);

    if(index.size() >0 ){
      for(int ii = 0; ii < bws_net.n_elem ; ii++){
        bw_net = bws_net(ii);
        double kernel_net =  kernel_func(d,bw_net) * alpha;
        for(int j = 0 ; j < bws_time.n_elem; j ++){
          bw_time = bws_time(j);
          // NOTE : we are not doing the bw scaling here but later
          for (int zz : index){
            double kernel_time = kernel_func(std::abs(v_time-time_events[zz]), bw_time);
            kvalues(ii,j,zz) = kvalues(ii,j,zz) + ((kernel_net * kernel_time));
          }
        }
      }
    }


    // we can now prepare the next steps
    if(max_bw_net >= d){
      IntegerVector v_neighbours = neighbour_list[v-1];
      int n = v_neighbours.length();
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
            int v2 = v_neighbours[j];
            int l2 = edge_mat(v,v2);
            double d2 = d + line_weights[l2-1];
            // first case, we must back fire
            if(v2 == prev_node){
              if(n>2){
                double p2 = (2.0-n)/n;
                double new_alpha = alpha * p2;
                acase new_case = {v2,v,new_depth,d2,new_alpha};
                //Rcout << "  adding this cas : d="<<d2<<", v2="<<v2<<", v3="<<v3<<", l2="<<l2<<" alpha="<<new_alpha<<"\n";
                //data_holder.push(new_case);
                data_holder.push_back(new_case);
              }
            }else{
              double new_alpha = alpha * (2.0/n);
              acase new_case = {v2,v,new_depth,d2,new_alpha};
              //Rcout << "  adding this cas : d="<<d2<<", v2="<<v2<<", v3="<<v3<<", l2="<<l2<<" alpha="<<new_alpha<<"\n";
              data_holder.push_back(new_case);
            }
          }
        }
      }
    }
  }

  // and now we can apply the scaling
  arma::mat scale_mat(bws_net.n_elem, bws_time.n_elem);

  for(int i = 0; i < bws_time.n_elem; i++){
    scale_mat.col(i) = 1.0/(bws_net * bws_time(i));
  }

  for(int i = 0; i < events.length(); i++){
    kvalues.slice(i) = kvalues.slice(i) % scale_mat;
  }


  return kvalues;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// THE EXPOSED FUNCTION
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//' @title The exposed function to calculate TNKDE likelihood cv
//' @name tnkde_get_loo_values
//' @param method a string, one of "simple", "continuous", "discontinuous"
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param sel_events a Numeric vector indicating the selected events (id of nodes)
//' @param sel_events_wid a Numeric Vector indicating the unique if of the selected events
//' @param sel_events_time a Numeric Vector indicating the time of the selected events
//' @param events a NumericVector indicating the nodes in graph beeing events
//' @param events_wid a NumericVector indicating the unique id of all the events
//' @param events_time a NumericVector indicating the timestamp of each event
//' @param weights a cube with the weights associated with each events for each
//' bws_net and bws_time.
//' @param bws_net an arma::vec with the network bandwidths to consider
//' @param bws_time an arma::vec with the time bandwidths to consider
//' @param kernel_name a string with the name of the kernel to use
//' @param line_list a DataFrame describing the lines
//' @param max_depth the maximum recursion depth
//' @param min_tol a double indicating by how much 0 in densities values must be replaced
//' @return a matrix with the CV score for each pair of bandiwdths
//' @export
//' @examples
//' # no example provided, this is an internal function
// [[Rcpp::export]]
arma::mat tnkde_get_loo_values(std::string method, List neighbour_list,
                               NumericVector sel_events, NumericVector sel_events_wid, NumericVector sel_events_time,
                               NumericVector events, NumericVector events_wid, NumericVector events_time,
                               arma::cube weights,
                               arma::vec bws_net, arma::vec bws_time, std::string kernel_name,
                               DataFrame line_list, int max_depth, double min_tol){

  //selecting the kernel function
  fptros kernel_func = select_kernelos(kernel_name);

  //step0 extract the columns of the dataframe
  NumericVector line_weights = line_list["weight"];

  //step 1 : mettre toutes les valeurs a 0 (bw_net 0 et bw_time 0)
  // NOTE WE calculate the values only for the events in sel_events
  arma::cube base_k(bws_net.n_elem, bws_time.n_elem, sel_events.length());

  //calculer la matrice des lignes
  //IntegerMatrix edge_mat = make_matrix(line_list,neighbour_list);
  arma::sp_mat edge_mat = make_matrix_sparse(line_list,neighbour_list);
  arma::cube k;
  //step2 : iterer sur chaque event de la zone d'etude
  int cnt_e = events.length()-1;
  for(int i=0; i <= cnt_e; ++i){
    //preparer les differentes valeurs de departs pour l'event y
    int y = events[i];
    int wid = events_wid[i];
    double v_time = events_time[i];
    // launching recursion
    // here we got the the influences of the vertex y on each other selected event in quadra
    if(method == "simple"){
      k = ess_kernel_loo_tnkde(kernel_func, edge_mat, sel_events,sel_events_wid,
                               sel_events_time, neighbour_list,
                               y, wid, v_time,
                               bws_net, bws_time,
                               line_weights, max_depth);

    }else if (method == "discontinuous"){
      k = esd_kernel_loo_tnkde(kernel_func, edge_mat, sel_events,sel_events_wid,
                                sel_events_time, neighbour_list,
                                y, wid, v_time,
                                bws_net, bws_time,
                                line_weights, max_depth);
    }else{
      k = esc_kernel_loo_tnkde(kernel_func, edge_mat, sel_events,sel_events_wid,
                                sel_events_time, neighbour_list,
                                y, wid, v_time,
                                bws_net, bws_time,
                                line_weights, max_depth);
    }
    // NOTE : the scaling by bws is applied above
    // and we must now remove its influence on itself
    arma::mat w = weights.slice(i);
    for(int ii = 0; ii < sel_events.length(); ii++){
      k.slice(ii) = k.slice(ii) % w;
    }

    // if y was a selected event, its own weight must be set to 0
    int index = get_first_index(sel_events_wid,wid);
    if(index >= 0){
      k.slice(index).fill(0);
    }

    // and summing the values at each iteration (% is the element wise product)
    base_k += k;
  };
  // and calculate the final values
  arma::mat result;

  arma::uvec neg_elems = arma::find(base_k <= 0);
  base_k.elem(neg_elems).fill(min_tol);

  result = arma::sum(arma::log(base_k),2);


  return result;
}





