#include "spNetwork.h"
#include "base_kernel_funtions.h"
#include "matrices_functions.h"

#include <cmath>


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// THE FUNCTIONS TO CALCULATE BW SELECTION CV LOO TEMPORAL FOR SIMPLE NKDE
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//' @title The worker function to calculate simple NKDE likelihood cv
//' @name ess_kernel_loo_nkde
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param events a NumericVector indicating the nodes in the graph being events
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param v the actual node to consider (int)
//' @param bws_net an arma::vec with the network bandwidths to consider
//' @param line_weights a vector with the length of the edges
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a matrix with the impact of the event v on each other events for
//' each pair of bandwidths (mat(event, bws_net))
//' @keywords internal
arma::mat ess_kernel_loo_nkde(fptros kernel_func, arma::sp_imat &edge_mat,
                                IntegerVector &events,
                                IntegerVector &events_wid,
                                List &neighbour_list,
                                int v, int wid,
                                arma::rowvec &bws_net,
                                NumericVector &line_weights, int max_depth){

  //step0 : generate the queue
  int depth = 0;
  std::queue <List> data_holder;
  arma::mat kvalues(events.length(), bws_net.n_elem);

    //step1 : generate the first case
  List cas1 = List::create(Named("d")=0.0,
                           Named("v") = v,
                           Named("prev_node") = -999,
                           Named("depth") = 0
  );
  data_holder.push(cas1);

  double max_bw_net = arma::max(bws_net);
  double bw_net;

  //lancement des iterations

  // NB : no speed gained here by using a struc or reducing the number of variable declaration

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
          std::vector<int> index = get_all_indeces_int(events,v2);
          if(index.size() >0 ){
            // il semble que v2 soit un noeud pour lequel au moins un evenement est present
            for(int ii = 0; ii < bws_net.n_elem ; ii++){
              bw_net = bws_net(ii);
              double kernel_net =  kernel_func(d2,bw_net);
              // NOTE : we are doing the bw scaling here
              for (int zz : index){
                kvalues(zz,ii) = kvalues(zz,ii) + (kernel_net * (1.0/bw_net));
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
//' @name esd_kernel_loo_nkde
//' @description The worker function to calculate discontinuous TNKDE likelihood cv (INTERNAL)
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param events a NumericVector indicating the nodes in the graph being events
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param v the actual node to consider (int)
//' @param bws_net an arma::vec with the network bandwidths to consider
//' @param line_weights a vector with the length of the edges
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a cube with the impact of the event v on each other events for
//' each pair of bandwidths (cube(bws_net, bws_time, events))
arma::mat esd_kernel_loo_nkde(fptros kernel_func, arma::sp_imat &edge_mat,
                                IntegerVector &events,
                                IntegerVector &events_wid,
                                List &neighbour_list,
                                int v, int wid,
                                arma::rowvec &bws_net,
                                NumericVector &line_weights, int max_depth){

  //step0 : generate the queue
  int depth = 0;
  std::queue <List> data_holder;

  arma::mat kvalues(events.length(), bws_net.n_elem);
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

  double bw_net;

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
          double d2 = line_weights[edge_id-1] + d;

          //est ce que v2 est un evenement pour lequel on doit garder la valeur

          std::vector<int> index = get_all_indeces_int(events,v2);
          if(index.size() >0 ){

            for(int ii = 0; ii < bws_net.n_elem ; ii++){
              bw_net = bws_net(ii);
              double kernel_net =  kernel_func(d2,bw_net) * new_alpha;
              for (int zz : index){
                kvalues(zz,ii) = kvalues(zz,ii) + (kernel_net * (1.0/bw_net));
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
//' @name esc_kernel_loo_nkde
//' @description The worker function to calculate continuous TNKDE likelihood cv (INTERNAL)
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param events a NumericVector indicating the nodes in the graph being events
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param v the actual node to consider (int)
//' @param bws_net an arma::vec with the network bandwidths to consider
//' @param line_weights a vector with the length of the edges
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a cube with the impact of the event v on each other events for
//' each pair of bandwidths (cube(bws_net, bws_time, events))
arma::mat esc_kernel_loo_nkde(fptros kernel_func, arma::sp_imat &edge_mat,
                                IntegerVector &events,
                                IntegerVector &events_wid,
                                List &neighbour_list,
                                int v, int wid,
                                arma::rowvec &bws_net,
                                NumericVector &line_weights, int max_depth){
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
  arma::mat kvalues(events.length(), bws_net.n_elem);
  double max_bw_net = arma::max(bws_net);
  double bw_net;

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

    std::vector<int> index = get_all_indeces_int(events,v);

    if(index.size() >0 ){
      for(int ii = 0; ii < bws_net.n_elem ; ii++){
        bw_net = bws_net(ii);
        double kernel_net =  kernel_func(d,bw_net) * alpha;
        // NOTE : we are not doing the bw scaling here but later
        for (int zz : index){
          kvalues(zz,ii) = kvalues(zz,ii) + kernel_net;

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
                double new_alpha = (alpha) * p2;
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
  arma::mat scale_mat(events.length(), bws_net.n_elem);

  for(int i = 0; i < bws_net.n_elem; i++){
    kvalues.col(i) = kvalues.col(i) * (1.0/(bws_net(i)));
  }


  return kvalues;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// THE EXPOSED FUNCTION
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//' @title The exposed function to calculate NKDE likelihood cv
//' @name nkde_get_loo_values
//' @description The exposed function to calculate NKDE likelihood cv (INTERNAL)
//' @param method a string, one of "simple", "continuous", "discontinuous"
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param sel_events a Numeric vector indicating the selected events (id of nodes)
//' @param sel_events_wid a Numeric Vector indicating the unique if of the selected events
//' @param events a NumericVector indicating the nodes in the graph being events
//' @param events_wid a NumericVector indicating the unique id of all the events
//' @param weights a matrix with the weights associated with each event (row) for each
//' bws_net (cols).
//' @param bws_net an arma::mat with the network bandwidths to consider for each event
//' @param kernel_name a string with the name of the kernel to use
//' @param line_list a DataFrame describing the lines
//' @param max_depth the maximum recursion depth
//' @param cvl a boolean indicating if the Cronie (TRUE) or CV likelihood (FALSE) must be used
//' @return a vector with the CV score for each bandwidth and the densities if required
//' @export
//' @examples
//' # no example provided, this is an internal function
// [[Rcpp::export]]
arma::mat nkde_get_loo_values(std::string method, List &neighbour_list,
                               IntegerVector &sel_events,
                               IntegerVector &sel_events_wid,
                               IntegerVector &events,
                               IntegerVector &events_wid,
                               arma::mat &weights,
                               arma::mat &bws_net, std::string kernel_name,
                               DataFrame &line_list, int max_depth, bool cvl){

  //selecting the kernel function
  fptros kernel_func = select_kernelos(kernel_name);

  //step0 extract the columns of the dataframe
  NumericVector line_weights = line_list["weight"];

  // NOTE WE calculate the values only for the events in sel_events
  arma::mat base_k(sel_events.length(), bws_net.n_cols);

  //calculer la matrice des lignes
  //IntegerMatrix edge_mat = make_matrix(line_list,neighbour_list);
  arma::sp_imat edge_mat = make_imatrix_sparse(line_list, neighbour_list);
  arma::mat k;
  //step2 : iterer sur chaque event de la zone d'etude
  int cnt_e = events.length()-1;
  for(int i=0; i <= cnt_e; ++i){
    //preparer les differentes valeurs de departs pour l'event y
    int y = events[i];
    int wid = events_wid[i];
    arma::rowvec ibws = bws_net.row(i);
    // launching recursion
    // here we got the the influences of the vertex y on each other selected event in quadra
    if(method == "simple"){
      k = ess_kernel_loo_nkde(kernel_func, edge_mat, sel_events,sel_events_wid,
                               neighbour_list,
                               y, wid,
                               ibws,
                               line_weights, max_depth);

    }else if (method == "discontinuous"){
      k = esd_kernel_loo_nkde(kernel_func, edge_mat, sel_events,sel_events_wid,
                               neighbour_list,
                               y, wid,
                               ibws,
                               line_weights, max_depth);
    }else{
      k = esc_kernel_loo_nkde(kernel_func, edge_mat, sel_events,sel_events_wid,
                               neighbour_list,
                               y, wid,
                               ibws,
                               line_weights, max_depth);
    }
    //Rcout << "Here is the k vector value : \n\n" << k <<"\n\n";
    //Rcout << "for the observation : \n\n" << i <<"\n\n";
    // NOTE : the scaling by bws is applied above
    // and we must now add the weight of the event
    // this weight can change according to the bw

    // Rcout << "Here are the weights in c++ : \n" <<  weights << "\n";

    for(int ii = 0; ii < k.n_cols; ii++){
      k.col(ii) = k.col(ii) * weights(i,ii);
      // Rcout << "applying weight : " << weights(i,ii) << "\n";
    }

    // if y was a selected event, its own weight must be set to 0 (if cvl == FALSE)
    // otherwise I must add to its own weight (if simple or discontinuous)
    int index = get_first_index_int(sel_events_wid,wid);
    if(index >= 0){
      if(cvl == false){
        k.row(index).fill(0.0);
      }else{
        if(method == "simple" || method == "discontinuous"){
          for(int ii = 0; ii < ibws.n_elem; ii++){
            //Rcout << "  Adding a bonus !"<<i<<"\n";
            k(index,ii) = k(index,ii) + (kernel_func(0,ibws(ii)) * (1/ibws(ii)));
          }
        }
      }
    }


    // and summing the values at each iteration (% is the element wise product)
    base_k += k;
  };

  // replacing 0 densities by mintol
  return base_k;

  // arma::uvec neg_elems = arma::find(base_k <= 0);
  // if(zero_strat == "min_double"){
  //   base_k.elem(neg_elems).fill(min_tol);
  // }else{
  //   base_k.elem(neg_elems).fill(1);
  // }

  // and calculating the final score
  // arma::colvec result;
  // if(cvl == false){
  //   result = arma::sum(arma::log(base_k),0).t();
  // }else{
  //   result = arma::sum(arma::pow(base_k,-1.0),0).t();
  // }
  // return result;
}


