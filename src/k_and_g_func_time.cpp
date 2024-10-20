#include "spNetwork.h"
#include "matrices_functions.h"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### k and g couting space-time function ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title c++ k and g function counting worker
//' @name kgfunc_time_counting
//' @description c++ k function counting (INTERNAL)
//' @param dist_mat_net A matrix with the distances between points on the network
//' @param dist_mat_time A matrix with the distances between points in time
//' @param wc The weight of the points represented by the columns (destinations)
//' @param wr The weight of the points represented by the rows (origins)
//' @param breaks_net A numeric vector with the distance to consider on network
//' @param breaks_time A numeric vector with the distance to consider in time
//' @param width_net The width of each donut for the network dimension
//' @param width_time The width of each donut for the time dimension
//' @param cross A boolean indicating if we are calculating a cross k function or not (default is FALSE)
//' @return A list  of two numeric matrices with the values of the k and g function evaluated at the required distances
//' @export
// [[Rcpp::export]]
List kgfunc_time_counting(arma::mat dist_mat_net,
                          arma::mat dist_mat_time,
                          arma::rowvec wc,
                          NumericVector wr,
                          NumericVector breaks_net,
                          NumericVector breaks_time,
                          float width_net,
                          float width_time,
                          double cross = false){

  float end_net = max(breaks_net);
  float end_time = max(breaks_time);

  // we need an extra dimension here to do the counting. For each point, for each bandwidth on the network
  // and in time, we will count how many other point it can reach.
  arma::cube counting_k(breaks_net.size(), breaks_time.size(), dist_mat_net.n_rows);
  arma::cube counting_g(breaks_net.size(), breaks_time.size(), dist_mat_net.n_rows);

  int cross_malus = 1;
  if(cross){
    cross_malus = 0;
  }

  // we iterate over each row of dist_mat
  for(int i = 0; i < dist_mat_net.n_rows; ++i) {
    Rcout << "iterating over the row : " << i << "\n";

    float wi = wr(i);
    arma::rowvec row_net = dist_mat_net.row(i);
    arma::rowvec row_time = dist_mat_time.row(i);

    // debugging
    NumericVector r1 = wrap(row_net);
    NumericVector r2 = wrap(row_time);
    Rcout << "distances over time : \n" << r2 <<"\n" ;
    Rcout << "distances over network : \n" << r1 <<"\n" ;

    // we do a first test to remove all the distances that are way to big
    arma::uvec test = arma::find( (row_net <= end_net + width_net) && (row_time <= end_time + width_time) )  ;

    arma::vec ok_d_net = row_net.elem(test) ;
    arma::vec ok_d_time = row_time.elem(test) ;

    // debugging
    r1 = wrap(ok_d_net);
    r2 = wrap(ok_d_time);
    Rcout << "ok distances over time : \n" << r2 <<"\n" ;
    Rcout << "ok distances over network : \n" << r1 <<"\n" ;

    arma::vec ok_w = wc.elem(test);

    Rcout << "ok weights : \n" << ok_w <<"\n" ;

    arma::vec ok_d2_net = row_net.elem(test);
    arma::vec ok_d2_time = row_time.elem(test);


    arma::vec ok_w2 = wc.elem(test);
    arma::vec ok_w1 = wc.elem(test);

    // we start by iterating on the network distance
    for(int z1 = 0; z1 < breaks_net.size(); ++z1) {

      float dist_net = breaks_net[z1];

      Rcout << "    iterating over the breaks_net : " << dist_net << "\n";

      // at each iteration on the network distances, we start by cutting all the
      // points that are too far and cannot be reached.

      // arma::uvec test1_net = ok_d_net <= (dist_net+width_net);
      // arma::uvec test2_net = ok_d_net >= (dist_net-width_net);
      // arma::uvec test3_net = ok_d_net <= (dist_net);

      // with the first part of the test, we can reduce the next research
      // this is only based on the network distance
      // arma::vec ok_d_net_z1 = ok_d_net.elem(arma::find(test1_net)) ;
      // arma::vec ok_d_time_z1 = ok_d_time.elem(arma::find(test1_net)) ;
      // arma::vec ok_w_z1 = ok_w.elem(arma::find(test1_net));

      // first col is net distance, then time distance, then weight
      arma::mat ok_z1(ok_d_net.n_elem,3);
      Rcout << "titi3 \n";
      ok_z1.col(2) = ok_w;
      Rcout << ok_z1 << "\n";
      Rcout << "titi1 \n";
      ok_z1.col(0) = ok_d_net;
      Rcout << ok_z1 << "\n";
      Rcout << "titi2 \n";
      ok_z1.col(1) = ok_d_time;
      Rcout << ok_z1 << "\n";

      // we filter it to remove points that are too far
      ok_z1 = ok_z1.rows(arma::find(ok_z1.col(0) <= dist_net+width_net));
      Rcout << "After filtering \n";
      Rcout << ok_z1 << "\n";

      // we precalculate the tests for distances and will filter
      // them later
      // arma::uvec test1_net = ok_z1.col(0) <= (dist_net+width_net);
      // arma::uvec test2_net = ok_z1.col(0) >= (dist_net-width_net);
      // arma::uvec test3_net = ok_z1.col(0) <= (dist_net);
      arma::umat mat_test_net = join_rows(ok_z1.col(0) <= (dist_net+width_net),
                                          ok_z1.col(0) >= (dist_net-width_net),
                                          ok_z1.col(0) <= (dist_net));



      // we will also iterate over the time distances
      for(int z2 = 0; z2 < breaks_time.size(); ++z2) {

        Rcout << "        iterating over the breaks_net : " << z2 << "\n";

        Rcout << "Here is the ok_z1 \n\n" << ok_z1 << "\n1";

        float dist_time = breaks_time[z2];


        // for the g func, we must check that dist - width is > 0, otherwise we must apply
        // a malus to not do self counting
        int malus = 0;
        if((dist_net-width_net <= 0) | (dist_time-width_time <= 0)){
          malus = 1;
        }

        // here I do the necessary checks for the temporal dimension
        arma::uvec test1_time = ok_z1.col(1) <= (dist_time+width_time);
        arma::uvec test2_time = ok_z1.col(1) >= (dist_time-width_time);
        arma::uvec test3_time = ok_z1.col(1) <= (dist_time);


        // we can the use them to do the test conjointly
        arma::uvec trueTest = arma::find((mat_test_net.col(0)) && (test1_time) && (mat_test_net.col(1)) && (test2_time));

        // with the part test I can get the value of the k function
        ok_w1 = ok_w.elem(arma::find(test3_time && mat_test_net.col(2))) ;
        // minus one here is important to remove self weight
        counting_k(z1, z2, i) = (arma::sum(ok_w1)-cross_malus) * wi;

        // with the full test we can get the values for the g function
        ok_w2 = ok_w.elem(trueTest) ;

        // for g we do not include minus one because the donut form will prevent self count
        counting_g(z1, z2, i) = (arma::sum(ok_w2) - (malus * cross_malus)) * wi;


        // we do not want to retest all the distances so we will reduce the data in ok_z1
        arma::uvec t1_filter = arma::find(ok_z1.col(1) <= dist_time+width_time);
        ok_z1 = ok_z1.rows(t1_filter);
        mat_test_net = mat_test_net.rows(t1_filter);

      }
    }

  }

  List results = List::create(counting_k, counting_g);

  return results;
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### base k space-time function ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title c++ k space-time function
//' @name k_nt_func_cpp
//' @param dist_mat_net A square matrix with the distances between points (network)
//' @param dist_mat_time A square matrix with the distances between points (time)
//' @param start_net A float, the start value for evaluating the k-function (network)
//' @param end_net A float, the last value for evaluating the k-function (network)
//' @param step_net A float, the jump between two evaluations of the k-function (network)
//' @param start_time A float, the start value for evaluating the k-function (time)
//' @param end_time A float, the last value for evaluating the k-function (time)
//' @param step_time A float, the jump between two evaluations of the k-function (time)
//' @param Lt The total length of the network
//' @param Tt The total duration of study area
//' @param n The number of points
//' @param w The weight of the points (coincident points)
//' @param cross a boolean indicating of we are calculating a cross k or g function
//' @keywords internal
// [[Rcpp::export]]
List k_g_nt_func_cpp2(arma::mat dist_mat_net, arma::mat dist_mat_time,
                            float start_net,float end_net, float step_net,
                            float start_time,float end_time, float step_time,
                            float width_net,
                            float width_time,
                            float Lt, float Tt, int n,
                            arma::rowvec wc,
                            NumericVector wr,
                            double cross = true
                            ){

  // we prepare the breaks for the network distances
  std::vector<float> breaks0 = seq_num3(start_net,end_net,step_net);
  std::reverse(breaks0.begin(), breaks0.end());
  NumericVector breaks_net = wrap(breaks0);

  // we prepare the breaks for the time distances
  breaks0 = seq_num3(start_time,end_time,step_time);
  std::reverse(breaks0.begin(), breaks0.end());
  NumericVector breaks_time = wrap(breaks0);

  float t1;
  if(cross){
    t1 = 1.0/((n)/ (Lt*Tt));
  }else{
    t1 = 1.0/((n-1)/(Lt*Tt));
  }

  width_net = width_net/2.0;
  width_time = width_time/2.0;

  // we start here by counting for each distance how many points are
  // reachable from each point
  // the produced matrice have a number of column equal to the number
  // of breaks, and a number of row equal to the number of point

  List elements = kgfunc_time_counting(dist_mat_net, dist_mat_time, wc, wr, breaks_net, breaks_time, width_net, width_time, cross);

  arma::cube counting_k = elements[0];
  arma::cube counting_g = elements[1];

  int pts_number = wr.size();

  arma::mat k_values(breaks_net.size(), breaks_time.size());
  arma::mat g_values(breaks_net.size(), breaks_time.size());

  // here, we must calculate the mean for each tube in the above elements

  for(int i = 0; i <  breaks_net.size(); i++){

    for(int j = 0; j <  breaks_time.size(); j++){

      float total = arma::accu(counting_k.tube(i,j));
      k_values(i,j) = total;

    }

  }


  // we must apply here the weights to the rows when calculating the means
  float div = sum(wr);

  k_values = (k_values / div) * t1;
  g_values = (g_values / div) * t1;

  // we must also put everything back to the good order !
  k_values = arma::flipud(arma::fliplr(k_values));
  g_values = arma::flipud(arma::fliplr(g_values));


  List final = List::create(k_values, g_values);

  return final;

}


