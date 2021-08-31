#include "spNetwork.h"


// A simple function to cut lines specified as matrix of coordinates and a
// vector of distances
// [[Rcpp::export]]
List cut_lines_at_distances_cpp(List lines, NumericVector dists){

  double d, dd, totald, dt, t;
  int i;
  int j;
  double x2,y2,x1,y1,x3,y3;
  List newList;

  // on commencer par iterer sur chacune des lignes
  for(i=0; i < lines.length(); ++i){
    d = dists(i);
    NumericMatrix line = lines(i);
    // on va ensuite iterer sur les points de cette ligne
    NumericVector okX;
    NumericVector okY;
    okX.push_back(line(0,0));
    okY.push_back(line(0,1));
    totald = 0;

    for(j=1; j < line.nrow(); ++j){
        // on calcule la longueur entre ce nvx point et le precedent
        x2 = line(j,0);
        y2 = line(j,1);
        x1 = line(j-1,0);
        y1 = line(j-1,1);
        dd = sqrt(pow(x1-x2,2) + pow(y1-y2,2));
        totald = totald + dd;
        // si on ne depasse pas la longueur totale, on continue
        if (totald < d){
          okX.push_back(x2);
          okY.push_back(y2);
        }else{
          // si on depasse pas la longueur totale, on calcule la position du dernier point
          dt = d - (totald - dd);
          t = dt / dd;
          x3 = (1-t)*x1 + x2*t;
          y3 = (1-t)*y1 + t*y2;
          okX.push_back(x3);
          okY.push_back(y3);
          break;
        }
    }

    //creating a new matrix
    NumericMatrix outmat(okX.length(),2);
    outmat(_,0) = okX;
    outmat(_,1) = okY;
    newList.push_back(clone(outmat));
  }
  return newList;
};
