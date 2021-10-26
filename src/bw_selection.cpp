#include "spNetwork.h"
#include "base_kernel_funtions.h"
#include "matrices_functions.h"


//########################################################################
//#### bandwidth selection by likelihood loo cross validation #####
//########################################################################

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
//' @keywords internal
NumericVector esc_kernel_loo(fptros kernel_func, NumericVector samples_k, List& neighbour_list, arma::sp_mat& edge_mat, int v, int v1, int l1, double d,double alpha, double bw, NumericVector& line_weights, NumericVector& events, int depth, int max_depth){

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
//' @keywords internal
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
//' @keywords internal
NumericVector esd_kernel_loo(fptros kernel_func, arma::sp_mat& edge_mat, NumericVector& events,
                             List& neighbour_list ,int v, double bw,
                             NumericVector& line_weights, int depth, int max_depth){
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
//' @keywords internal
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
//' @keywords internal
NumericVector ess_kernel_loo(fptros kernel_func, arma::sp_mat& edge_mat, NumericVector& events,
                             List& neighbour_list ,int v, double bw,
                             NumericVector& line_weights, int depth, int max_depth){
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
//' @keywords internal
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
    // here we got the the influences of the vertex y on each other vertices
    NumericVector k = ess_kernel_loo(kernel_func,edge_mat, samples, neighbour_list ,y, bw, line_weights, depth,max_depth);
    // getting the actual base_k values (summed at each iteration)
    base_k += (k*w);
  };
  return base_k;
}
