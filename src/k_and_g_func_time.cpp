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
//' @return A list  of two numeric cubes with the values of the k and g function evaluated at the required distances
//' @keywords internal
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
                          bool cross = false){

  float end_net = max(breaks_net);
  float end_time = max(breaks_time);

  // we need an extra dimension here to do the counting. For each point, for each bandwidth on the network
  // and in time, we will count how many other point it can reach.

  arma::cube counting_k(breaks_net.size(), breaks_time.size(), dist_mat_net.n_rows);
  arma::cube counting_g(breaks_net.size(), breaks_time.size(), dist_mat_net.n_rows);

  // NOTE : when we are in cross mode, we must negate the malus
  int cross_malus = 1;
  if(cross){
    cross_malus = 0;
  }


  // we iterate over each row of dist_mat
  for(int i = 0; i < dist_mat_net.n_rows; ++i) {

    float wi = wr(i);

    arma::mat idata(3, dist_mat_net.n_cols);

    idata.row(0) = dist_mat_net.row(i);
    idata.row(1) = dist_mat_time.row(i);
    idata.row(2) = wc;
    idata = idata.t();



    // we do a first test to remove all the distances that are way to big
    idata = idata.rows(arma::find( (idata.col(0) <= end_net + width_net) && (idata.col(1) <= end_time + width_time) ));


    // we start by iterating on the network distance
    for(int z1 = 0; z1 < breaks_net.size(); ++z1) {

      float dist_net = breaks_net[z1];

      // at each iteration on the network distances, we start by cutting all the
      // points that are too far and cannot be reached.

      // we filter it to remove points that are too far
      arma::mat ok_z1 = idata.rows(arma::find(idata.col(0) <= dist_net+width_net));

      // we might end up with an empty ok_z1...
      // we must just skip the remaining part in that case

      if(ok_z1.n_rows > 0){

        // NOTE : we are precalculating the tests.
        // for each other point in the network distance matrix, we check if
        // it is in the donught range and below the max range (for g and k function)
        // the results are stored in a matrix. Each col correspond to a test
        arma::umat mat_test_net = join_rows(ok_z1.col(0) <= (dist_net+width_net),
                                            ok_z1.col(0) >= (dist_net-width_net),
                                            ok_z1.col(0) <= (dist_net));

        // we will also iterate over the time distances
        for(int z2 = 0; z2 < breaks_time.size(); ++z2) {

          float dist_time = breaks_time[z2];

          // here I do the necessary checks for the temporal dimension
          arma::uvec test1_time = ok_z1.col(1) <= (dist_time+width_time);
          arma::uvec test2_time = ok_z1.col(1) >= (dist_time-width_time);
          arma::uvec test3_time = ok_z1.col(1) <= (dist_time);

          // HERE I HAVE TO DEAL WITH THE MALUS FOR THE G FUNCTION TO AVOID SELF COUNTING
          // It is uncertain that we will do self-couting when we are using the g function
          // because of its donught form.
          // we must apply the malus only if the lower bound is below zero for both
          // time and network
          int malus = 0;
          if( (dist_time-width_time <= 0) & (dist_net-width_net <= 0)){
            malus = 1;
          }

          // we can the use them to do the test conjointly
          arma::uvec trueTest = arma::find(
            (mat_test_net.col(0)) && (test1_time) &&
              (mat_test_net.col(1)) && (test2_time));


          // with the part test I can get the value of the k function
          arma::vec ok_w = ok_z1.col(2);
          arma::vec ok_w1 = ok_w.elem(arma::find(test3_time && mat_test_net.col(2)));

          // idem, we must test that we do not have an empty vector here
          if(ok_w1.n_elem > 0){
            counting_k(z1, z2, i) = (arma::sum(ok_w1) - cross_malus) * wi;
          }

          // with the full test we can get the values for the g function
          arma::vec  ok_w2 = ok_w.elem(trueTest) ;

          if(ok_w2.n_elem > 0){
            // for g we do not include minus one because the donut form will prevent self count
            counting_g(z1, z2, i) = (arma::sum(ok_w2) - (malus * cross_malus)) * wi;
          }


          // we do not want to retest all the distances so we will reduce the data in ok_z1
          arma::uvec t1_filter = arma::find(ok_z1.col(1) <= dist_time+width_time);
          ok_z1 = ok_z1.rows(t1_filter);
          mat_test_net = mat_test_net.rows(t1_filter);


          // we can break the loop if there is nothing left in ok_z1
          if(ok_z1.n_rows == 0){
            break;
          }

        }
      }


    }

  }

  List results = List::create(counting_k, counting_g);

  return results;
}



//' @title c++ k function counting worker
//' @name kgfunc_time_counting
//' @description c++ k function counting (INTERNAL)
//' @param dist_mat_net A matrix with the distances between points on the network
//' @param dist_mat_time A matrix with the distances between points in time
//' @param wc The weight of the points represented by the columns (destinations)
//' @param wr The weight of the points represented by the rows (origins)
//' @param breaks_net A numeric vector with the distance to consider on network
//' @param breaks_time A numeric vector with the distance to consider in time
//' @param cross A boolean indicating if we are calculating a cross k function or not (default is FALSE)
//' @return A list  of two numeric cubes with the values of the k and g function evaluated at the required distances
//' @export
// [[Rcpp::export]]
arma::cube kfunc_time_counting(arma::mat dist_mat_net,
                          arma::mat dist_mat_time,
                          arma::rowvec wc,
                          NumericVector wr,
                          NumericVector breaks_net,
                          NumericVector breaks_time,
                          bool cross = false){

  float end_net = max(breaks_net);
  float end_time = max(breaks_time);

  // we need an extra dimension here to do the counting. For each point, for each bandwidth on the network
  // and in time, we will count how many other point it can reach.
  arma::cube counting_k(breaks_net.size(), breaks_time.size(), dist_mat_net.n_rows);

  int cross_malus = 1;
  if(cross){
    cross_malus = 0;
  }

  // we iterate over each row of dist_mat
  for(int i = 0; i < dist_mat_net.n_rows; ++i) {

    float wi = wr(i);
    arma::rowvec row_net = dist_mat_net.row(i);
    arma::rowvec row_time = dist_mat_time.row(i);

    // we do a first test to remove all the distances that are way to big
    arma::uvec test = arma::find( (row_net <= end_net) && (row_time <= end_time) )  ;

    arma::vec ok_d_net = row_net.elem(test) ;
    arma::vec ok_d_time = row_time.elem(test) ;
    arma::vec ok_w = wc.elem(test);


    arma::vec ok_d2_net = row_net.elem(test);
    arma::vec ok_d2_time = row_time.elem(test);

    arma::vec ok_w2 = wc.elem(test);
    arma::vec ok_w1 = wc.elem(test);


    // we start by iterating on the network distance
    for(int z1 = 0; z1 < breaks_net.size(); ++z1) {

      float dist_net = breaks_net[z1];

      // we can break the calculation here if ok_d_net is empty
      // this means that no distances within the two distance matrices are below the maximum given time or network distance
      if(row_net.n_elem == 0){
        break;
      }


      // first col is net distance, then time distance, then weight
      arma::mat ok_z1(ok_d_net.n_elem,3);
      ok_z1.col(0) = ok_d_net;
      ok_z1.col(1) = ok_d_time;
      ok_z1.col(2) = ok_w;

      // we filter it to remove points that are too far on the network
      ok_z1 = ok_z1.rows(arma::find(ok_z1.col(0) <= dist_net));


      // we will also iterate over the time distances
      for(int z2 = 0; z2 < breaks_time.size(); ++z2) {

        float dist_time = breaks_time[z2];

        // we filter it to remove points that are too far in time
        ok_z1 = ok_z1.rows(arma::find(ok_z1.col(1) <= (dist_time)));

        // we can break the calculation here if ok_d_net is empty
        // this means that no distances within the two distance matrices are below the maximum given time or network distance
        if(ok_z1.n_elem == 0){
          break;
        }

        // now we can calculte the k value
        // minus one here is important to remove self weight
        counting_k(z1, z2, i) = (arma::sum(ok_z1.col(2)) - cross_malus) * wi;

      }
    }

  }

  return counting_k;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### base k space-time function ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title c++ k and g space-time function
//' @name k_nt_func_cpp
//' @param dist_mat_net A square matrix with the distances between points (network)
//' @param dist_mat_time A square matrix with the distances between points (time)
//' @param start_net A float, the start value for evaluating the k-function (network)
//' @param end_net A float, the last value for evaluating the k-function (network)
//' @param step_net A float, the jump between two evaluations of the k-function (network)
//' @param width_net A float indicating the width of the donught of the g-function (network)
//' @param start_time A float, the start value for evaluating the k-function (time)
//' @param end_time A float, the last value for evaluating the k-function (time)
//' @param step_time A float, the jump between two evaluations of the k-function (time)
//' @param width_time A float indicating the width of the donught of the g-function (time)
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
                            bool cross = false
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

  List elements = kgfunc_time_counting(dist_mat_net, dist_mat_time, wc, wr, breaks_net, breaks_time, width_net, width_time);

  arma::cube counting_k = elements[0];
  arma::cube counting_g = elements[1];

  // Rcout << "HERE IS THE COUNTING G\n\n";
  // Rcout << counting_g;


  int pts_number = wr.size();

  arma::mat k_values(breaks_net.size(), breaks_time.size());
  arma::mat g_values(breaks_net.size(), breaks_time.size());

  // here, we must calculate the mean for each tube in the above elements

  for(int i = 0; i <  breaks_net.size(); i++){

    for(int j = 0; j <  breaks_time.size(); j++){

      k_values(i,j) = arma::accu(counting_k.tube(i,j));
      g_values(i,j) = arma::accu(counting_g.tube(i,j));

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
arma::mat k_nt_func_cpp2(arma::mat dist_mat_net, arma::mat dist_mat_time,
                      float start_net,float end_net, float step_net,
                      float start_time,float end_time, float step_time,
                      float Lt, float Tt, int n,
                      arma::rowvec wc,
                      NumericVector wr,
                      bool cross = false){

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

  // we start here by counting for each distance how many points are
  // reachable from each point
  // the produced matrice have a number of column equal to the number
  // of breaks, and a number of row equal to the number of point

  arma::cube counting_k = kfunc_time_counting(dist_mat_net, dist_mat_time, wc, wr, breaks_net, breaks_time, cross);


  int pts_number = wr.size();

  arma::mat k_values(breaks_net.size(), breaks_time.size());

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

  // we must also put everything back to the good order !
  k_values = arma::flipud(arma::fliplr(k_values));


  return k_values;

}
