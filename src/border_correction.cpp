#include "spNetwork.h"
#include "base_kernel_funtions.h"
#include "matrices_functions.h"


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
//' @keywords internal
// [[Rcpp::export]]
List corrfactor_discontinuous_sparse(List &neighbour_list,
                                     IntegerVector &events,
                                     DataFrame &line_list,
                                     NumericVector &bws,
                                     int max_depth){

  //extraire le poids des lignes
  arma::vec line_weights =  as<arma::vec>(line_list["weight"]);
  //preparer la matrice de voisinage
  arma::sp_imat edge_mat = make_imatrix_sparse(line_list,neighbour_list);
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
    //queue<List> data_holder;
    std::vector<List> data_holder;
    List start = List::create(Named("v") = y ,
                              Named("d") = 0.0,
                              Named("depth") = 0,
                              Named("prev_node") = -1,
                              Named("alpha") = 1.0);
    data_holder.push_back(start);

    while(data_holder.empty()==FALSE){

      //extraire le 1er cas
      List cas = data_holder.back();
      data_holder.pop_back();
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
              data_holder.push_back(new_cas);
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
//' @keywords internal
// [[Rcpp::export]]
List corrfactor_discontinuous(List &neighbour_list,
                              IntegerVector &events,
                              DataFrame &line_list,
                              NumericVector &bws,
                              int max_depth){

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
    //queue<List> data_holder;
    std::vector<List> data_holder;
    List start = List::create(Named("v") = y ,
                              Named("d") = 0.0,
                              Named("depth") = 0,
                              Named("prev_node") = -1,
                              Named("alpha") = 1.0);
    data_holder.push_back(start);

    while(data_holder.empty()==FALSE){

      //extraire le 1er cas
      List cas = data_holder.back();
      data_holder.pop_back();
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
              data_holder.push_back(new_cas);
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
//' @keywords internal
// [[Rcpp::export]]
List corrfactor_continuous_sparse(List &neighbour_list,
                                  IntegerVector &events,
                                  DataFrame &line_list,
                                  NumericVector &bws,
                                  int max_depth){

  //extraire le poids des lignes
  arma::vec line_weights =  as<arma::vec>(line_list["weight"]);
  //preparer la matrice de voisinage
  arma::sp_imat edge_mat = make_imatrix_sparse(line_list,neighbour_list);
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
    //queue<List> data_holder;
    std::vector<List> data_holder;
    List start = List::create(Named("v") = y ,
                              Named("d") = 0.0,
                              Named("depth") = 0,
                              Named("prev_node") = -1,
                              Named("alpha") = 1.0);
    data_holder.push_back(start);

    while(data_holder.empty()==FALSE){

      //extraire le 1er cas
      List cas = data_holder.back();
      data_holder.pop_back();
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
              data_holder.push_back(new_cas);

            }
          }// cas special pour le back fire !
          else{
            //on ne peut faire un backfire qu'a une vraie intersection
            if(cnt_n>2){
              //IntegerVector l_neighbours = neighbour_list[v2-1];
              //int n2 = l_neighbours.length();
              //double p2 = (n2-2.0)/n2;
              double p2 = (cnt_n-2.0)/cnt_n;
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
                data_holder.push_back(new_cas);
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
//' @keywords internal
// [[Rcpp::export]]
List corrfactor_continuous(List &neighbour_list,
                           IntegerVector &events,
                           DataFrame &line_list,
                           NumericVector &bws,
                           int max_depth){

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
    //queue<List> data_holder;
    std::vector<List> data_holder;
    List start = List::create(Named("v") = y ,
                              Named("d") = 0.0,
                              Named("depth") = 0,
                              Named("prev_node") = -1,
                              Named("alpha") = 1.0);
    data_holder.push_back(start);

    while(data_holder.empty()==FALSE){

      //extraire le 1er cas
      List cas = data_holder.back();
      data_holder.pop_back();
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
              data_holder.push_back(new_cas);

            }
          }// cas special pour le back fire !
          else{
            //on ne peut faire un backfire qu'a une vraie intersection
            if(cnt_n>2){
              //IntegerVector l_neighbours = neighbour_list[v2-1];
              //int n2 = l_neighbours.length();
              //double p2 = (n2-2.0)/n2;
              double p2 = (cnt_n-2.0)/cnt_n;
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
                data_holder.push_back(new_cas);
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
