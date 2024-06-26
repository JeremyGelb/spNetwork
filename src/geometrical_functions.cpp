// ******************************************************************
// **** We define here some functions to perform specific geometrical
// **** operations which need C++ acceleration
// ******************************************************************
#include "spNetwork.h"
#include "matrices_functions.h"


// some boost libraries used to building an rtree
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/algorithms/distance.hpp>
#include <boost/foreach.hpp>
#include <boost/geometry/index/predicates.hpp>

// **** some namespaces related to boost ****
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// ***** some global types used ****

// basic geometries
typedef bg::model::d2::point_xy<double> point_t;
//typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
typedef bg::model::linestring<point_t> linestring_t;
typedef bg::model::box<point_t> bbox;

// a simple vector of boost lines
typedef std::vector<linestring_t> lines_vector;


// *********************************************************************
// Get coordinates of a point
// *********************************************************************
double getX(point_t pt){
  return bg::get<0>(pt);
}

double getY(point_t pt){
  return bg::get<1>(pt);
}


// *********************************************************************
// Finding the nearest line for a set of points
// *********************************************************************

/*
 *
 * In this section, I implement a set of functions aiming to find for each
 * point in a set of points (defined as a 2 column numeric matrix) its nearest
 * line in a set of lines (defined as a list of 2 column numeric matrices)
 *
 */


// creating a Linestring from a matrix of coordinates
// somthing simple
linestring_t line_from_coords(NumericMatrix coords){
  // iterating on the coordinates to create a line
  linestring_t my_line;
  int j;
  int rj = coords.rows();
  for(j = 0; j < rj; j++){
    bg::append(my_line, point_t(coords(j,0),coords(j,1)));
  }
  return my_line;
}

// Creating a std::vector containing boost Linestrings from an
// Rcpp list of numeric matrices
// Note: the case of multilinestring is not supported
// If needed, I could make it work by using WKT instead of coordinates
lines_vector lines_vector_from_coordinates(List lines){
  int i;
  lines_vector my_lines;
  // iterating on the list of coordinates
  int ll = lines.length();
  for(i=0 ; i < ll; i++){
    NumericMatrix line = lines(i);
    linestring_t my_line = line_from_coords(line);
    // append the line to the vector
    my_lines.push_back(my_line);
  }
  return my_lines;

};


// Creating a spatial index for lines
// see here: https://www.boost.org/doc/libs/1_76_0/libs/geometry/doc/html/geometry/spatial_indexes/rtree_quickstart.html
// and here for enveloppe https://www.boost.org/doc/libs/1_71_0/libs/geometry/doc/html/geometry/reference/algorithms/envelope/envelope_2.html
// the basic element here will be a pair of box;int
// another lib if necessary : https://github.com/nushoin/RTree/blob/master/RTree.h
typedef std::pair<bbox, int> rtree_element;
typedef bgi::rtree< rtree_element, bgi::quadratic<16>> lines_rtree;

lines_rtree build_rtree_for_lines(lines_vector lines){

  // **** This is the classical way ****//
  // setting the empty rtree
  lines_rtree mytree;
  //boost::geometry::model::box<point_t> abox;
  bbox abox;

  // iterating over the lines
  int i = 0;
  int ls = lines.size();
  for(i = 0 ; i < ls ; i++){
    bg::envelope(lines[i], abox);
    rtree_element el = std::make_pair(abox, i);
    mytree.insert(el);
  }

  // **** If I use a range adaptor, the packing algo can be used and is supposed to be faster ****//
  // but is causes problem with the c++ compiler on ubuntu...
  // std::vector<rtree_element> boxes;
  //
  // size_t id_gen = 0;
  // std::transform(
  //   lines.begin(), lines.end(),
  //   back_inserter(boxes),
  //   [&](linestring_t const& l) {
  //     boost::geometry::model::box<point_t> abox;
  //     bg::envelope(l, abox);
  //     return std::make_pair(abox, id_gen++);
  //   }
  // );
  //bgi::rtree<rtree_element, bgi::quadratic<16> > rtree(boxes);

  // returning the filled rtree
  return mytree;

};


// Perform a query on an index created by the above function.
// The goal is to find some close lines in a given radius
// so it is necessary to check after a query if the real objects are
// in the distance calculated for their boxes.
// To do so, a while loop is used with a maximum number of iterations
// If no geometries are found within the increasing radius, an empty
// vector is returned
typedef std::vector<rtree_element> vector_rtree_element ;

vector_rtree_element find_close_lines_in_index(lines_rtree index, lines_vector lines, point_t point, double min_dist, int max_iter){

  // defining some variables
  bool ok = false;
  double actual_dist = min_dist/2.0;
  double width;
  vector_rtree_element returned_values;
  int iter = 0;

  while(ok == false){
    actual_dist = actual_dist*2.0;
    //returned_values.empty();
    returned_values.clear();
    // calculating the box
    width = actual_dist/2.0;
    point_t mypt1(point.x()-width, point.y()-width);
    //point_t mypt1(getX(point)-width, getY(point)-width);
    point_t mypt2(point.x()+width, point.y()+width);
    //point_t mypt2(getX(point)+width, getY(point)+width);
    bbox region(mypt1,mypt2);

    //querying the index
    index.query(bgi::intersects(region), std::back_inserter(returned_values));

    if(returned_values.size() > 0){
      // if we found boxes, we have to check the distance with geometries
      BOOST_FOREACH(rtree_element const& el, returned_values){
        linestring_t my_line = lines[el.second];
        double dist = bg::distance(point, my_line);
        if(dist < actual_dist){
          // if a geometry is really that close, then it is all good
          ok = true;
        }
      }
    }
    iter++;
    if(iter == max_iter){
      break;
    }
  }

  return returned_values;

}


// The global function used in this section.
// It takes as inputs a NumericMatrix to define the points and and List of
// NumericMatrix to define the LineStrings
// [[Rcpp::export]]
IntegerVector find_nearest_object_in_line_rtree(NumericMatrix pts, List lines, double min_dist, int max_iter){

  // prepare some objects
  linestring_t my_line;

  // step1: reading the lines as a vector of lines
  lines_vector vec_lines = lines_vector_from_coordinates(lines);

  // step2: creating the lines spatial index
  lines_rtree lines_index = build_rtree_for_lines(vec_lines);

  // step3: querying the spatial index for each point
  IntegerVector final_indexes;

  int i;
  for(i = 0; i<pts.nrow(); i++){
    // creating the point

    point_t mypt(pts(i,0),pts(i,1));

    vector_rtree_element returned_values = find_close_lines_in_index(lines_index, vec_lines, mypt, min_dist, max_iter);

    // if did not find candidates
    if(returned_values.size() == 0){
      final_indexes.push_back(-1);
    }else{
      // otherwise we have to iterate over them and find the best candidate
      int best_line = returned_values[0].second;
      double best_dist = bg::distance(mypt, vec_lines[best_line]);

      BOOST_FOREACH(rtree_element const& el, returned_values){
        linestring_t my_line = vec_lines[el.second];
        double dist = bg::distance(mypt, my_line);
        if(dist < best_dist){
          best_line = el.second;
          best_dist = dist;
        }
      }
      final_indexes.push_back(best_line);
    }
  }
  return final_indexes;

};


// *********************************************************************
// Line related geometric operations
// *********************************************************************


/*
 * Some functions used to work on lines
 *
 */

// A simple function to calculate the length of a line based on its coordinates
// [[Rcpp::export]]
float calc_line_length(NumericMatrix line){

  float full_length = 0.0;
  double x2,y2,x1,y1,x3,y3,dd;
  int j;
  for(j=1; j < line.nrow(); ++j){
    x2 = line(j,0);
    y2 = line(j,1);
    x1 = line(j-1,0);
    y1 = line(j-1,1);
    dd = sqrt(pow(x1-x2,2) + pow(y1-y2,2));
    full_length = full_length + dd;
  }
  return(full_length);

};


// A simple function to cut a line at a specified distance from its start
// [[Rcpp::export]]
NumericMatrix cut_line_at_dist(NumericMatrix line, float d){

  int j;
  double x2,y2,x1,y1,x3,y3,dd,dt,t;
  NumericVector okX;
  NumericVector okY;
  okX.push_back(line(0,0));
  okY.push_back(line(0,1));
  double totald = 0;

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
  return(outmat);
};


// A simple function to cut lines specified as matrix of coordinates and a
// vector of distances
// [[Rcpp::export]]
List cut_lines_at_distances_cpp(List lines, NumericVector dists){

  double d, dd, totald, dt, t;
  int i;
  int j;
  double x2,y2,x1,y1,x3,y3;
  //List newList;
  std::vector<NumericMatrix> newList;

  // on commencer par iterer sur chacune des lignes
  for(i=0; i < lines.length(); ++i){
    d = dists(i);
    NumericMatrix line = lines(i);
    NumericMatrix outmat = cut_line_at_dist(line, d);
    // on va ensuite iterer sur les points de cette ligne
    // NumericVector okX;
    // NumericVector okY;
    // okX.push_back(line(0,0));
    // okY.push_back(line(0,1));
    // totald = 0;
    //
    // for(j=1; j < line.nrow(); ++j){
    //   // on calcule la longueur entre ce nvx point et le precedent
    //   x2 = line(j,0);
    //   y2 = line(j,1);
    //   x1 = line(j-1,0);
    //   y1 = line(j-1,1);
    //   dd = sqrt(pow(x1-x2,2) + pow(y1-y2,2));
    //   totald = totald + dd;
    //   // si on ne depasse pas la longueur totale, on continue
    //   if (totald < d){
    //     okX.push_back(x2);
    //     okY.push_back(y2);
    //   }else{
    //     // si on depasse pas la longueur totale, on calcule la position du dernier point
    //     dt = d - (totald - dd);
    //     t = dt / dd;
    //     x3 = (1-t)*x1 + x2*t;
    //     y3 = (1-t)*y1 + t*y2;
    //     okX.push_back(x3);
    //     okY.push_back(y3);
    //     break;
    //   }
    // }
    //
    // //creating a new matrix
    // NumericMatrix outmat(okX.length(),2);
    // outmat(_,0) = okX;
    // outmat(_,1) = okY;
    newList.push_back(clone(outmat));
  }
  return wrap(newList);
};


// A simple function to cut lines specified as matrix of coordinates and a
// vector of distances
// [[Rcpp::export]]
NumericMatrix trim_line_for_isos(NumericMatrix line, float start_dist, float end_dist, bool donught, float d, float dd){

  // Je vais avoir besoin de la longueur de cette ligne
  float line_length = calc_line_length(line);
  float diff;
  NumericMatrix new_mat;

  // on complete les longueur si on ne les a pas encore
  if(start_dist < 0){
    start_dist = end_dist + line_length;
  }
  if(end_dist < 0){
    end_dist = start_dist + line_length;
  }

  // 1er cas, la distance pour le point de depart
  // est trop grande, il faut racourcir au depart
  if (start_dist > d){
    diff = start_dist - d;
    new_mat = reverseByRow(line);
    new_mat = cut_line_at_dist(new_mat, (line_length - diff));
    line = reverseByRow(new_mat);
  }
  // 2e cas, la distance pour le point d'arrivee
  // est trop grande, il faut racourcir au points d'arrivee
  if (end_dist > d){
    diff = end_dist - d;
    line_length = calc_line_length(line);
    new_mat = cut_line_at_dist(line, (line_length - diff));
    line = new_mat;
  }
  // si on est dans le cas des donught, il faut aussi faire la seconde inspection
  if(donught){
    // 1er cas la distance pour le point de depart est trop courte !
    // il faut couper ce qui est en trop
    if (start_dist < dd){
      diff = dd - start_dist;
      line_length = calc_line_length(line);
      new_mat = reverseByRow(line);
      new_mat = cut_line_at_dist(new_mat, (line_length - diff));
      line = reverseByRow(new_mat);
    }
    // 2 cas la distance pour le point d'arrivee est trop courte !
    // il faut couper ce qui est en trop
    if (end_dist < dd){
      diff = dd - end_dist;
      line_length = calc_line_length(line);
      new_mat = cut_line_at_dist(line, (line_length - diff));
      line = new_mat;
    }
  }
  return (line);
};


// A simple function to cut lines specified as matrix of coordinates and a
// vector of distances
// [[Rcpp::export]]
List trim_lines_for_isos_cpp(List lines, NumericVector start_dists, NumericVector end_dists, bool donught, float d, float dd){

  float start_dist, end_dist, diff, line_length;
  //List newList;
  std::vector<NumericMatrix> newList;
  int i;
  // on va commencer par iterer sur chacune des lignes
  for(i=0; i < lines.length(); ++i){
    NumericMatrix LineCoords = lines(i);
    // et on applique simplement la fonction d'au-dessus
    start_dist = start_dists(i);
    end_dist = end_dists(i);
    NumericMatrix NewLine = trim_line_for_isos(LineCoords,
                                               start_dist, end_dist,
                                               donught, d, dd);
    newList.push_back(clone(NewLine));

  }
  return wrap(newList);
};


// *********************************************************************
// Projecting a point on a line
// *********************************************************************

double points_dot_product(point_t p1, point_t p2){
  return(p1.x()*p2.x()+p1.y()*p2.y());
  //return(getX(p1)*getX(p2)+getY(p1)*getY(p2));
};


/*
 * Project a point P onto a segment L
 */
point_t project_point_on_segment(point_t p, linestring_t line){
  point_t v1 = line[0];
  point_t v2 = line[1];
  // get dot product of e1, e2
  //point_t e1(getX(v2) - getX(v1), getY(v2) - getY(v1));
  point_t e1(v2.x() - v1.x(), v2.y() - v1.y());
  //point_t e2(getX(p) - getX(v1), getY(p) - getY(v1));
  point_t e2(p.x() - v1.x(), p.y() - v1.y());
  double valDp = points_dot_product(e1, e2);
  // get length of vectors
  //double lenLineE1 = sqrt(getX(e1) * getX(e1) + getY(e1) * getY(e1));
  double lenLineE1 = sqrt(e1.x() * e1.x() + e1.y() * e1.y());
  //double lenLineE2 = sqrt(getX(e2) * getX(e2) + getY(e2) * getY(e2));
  double lenLineE2 = sqrt(e2.x() * e2.x() + e2.y() * e2.y());
  double cos = valDp / (lenLineE1 * lenLineE2);
  // length of v1P'
  double projLenOfLine = cos * lenLineE2;
  // point_t p2((getX(v1) + (projLenOfLine * getX(e1)) / lenLineE1),
  //            (getY(v1) + (projLenOfLine * getY(e1)) / lenLineE1));
  point_t p2((v1.x() + (projLenOfLine * e1.x()) / lenLineE1),
             (v1.y() + (projLenOfLine * e1.y()) / lenLineE1));
  return p2;
};


/*
 * Project a point P onto a LineString L and get the distance from the start of
 * the linestring, the result is returned as a pair : <dist,point>
 */

std::pair<double, point_t> project_point_on_Linestring_distance(point_t p, linestring_t line){

  // step1 : finding the best segment
  //Rcout << "Projecting a point on a line\n";
  int nb_seg = line.size();
  //Rcout << "The line has "<<nb_seg<<" vertices\n";
  int i;
  int best_candidate = 0;
  double best_dist = bg::distance(p,line)*2;
  //Rcout << "This is the best expected distance "<<best_dist<<"\n\n";
  double dist;
  double cumdist = 0;
  linestring_t seg;

  for(i = 0; i < (nb_seg-1); i++){
    //Rcout << "   Iterating on candidate "<<i<<"\n";
    seg.clear();
    seg.push_back(line[i]);
    seg.push_back(line[i+1]);
    dist = bg::distance(p,seg);
    //Rcout << "   the distance is "<<dist<<"\n";
    if(dist <= best_dist){
      //Rcout << "   it was a good candidate !\n\n";
      best_dist = dist;
      best_candidate = i;
      break;
    }
    cumdist += bg::distance(line[i],line[i+1]);
  }
  // the cumulative distance is calculated untill the first vertex of the selected segment !

  // step3 : projecting the point on the best candidate
  std::pair<double, point_t> result;
  linestring_t best_segment;
  point_t p1 = line[best_candidate];
  point_t p2 = line[best_candidate+1];
  float d1 = bg::distance(p1,p);
  float d2 = bg::distance(p2,p);

  if(d1 == 0.0){
    // cas 1 : le point est sur la vertex de depart
    result = std::make_pair(cumdist,p1);
  }else if (d2 == 0.0){
    // cas 2 : le point est sur la vertex de fin
    result = std::make_pair((cumdist + bg::distance(p2,p1)),p2);
  }else{
    // cas 3 : le point est quelque part entre les deux
    best_segment.push_back(p1);
    best_segment.push_back(p2);

    point_t projpt = project_point_on_segment(p, best_segment);

    //step 4 : calculating the final dist
    cumdist+=bg::distance(projpt, line[best_candidate]);
    result = std::make_pair(cumdist,projpt);
  }

  return result;
}


/*
 * Project a point P onto a LineString L and get the point as output
 */

point_t project_point_on_Linestring_point(point_t p, linestring_t line){

  // step1 : finding the best segment
  int nb_seg = line.size()-1;
  int i;
  int best_candidate = 0;
  double best_dist = bg::distance(p,line)*2;
  double dist;
  linestring_t seg;
  for(i = 0; i < nb_seg; i++){
    //seg.empty();
    seg.clear();
    seg.push_back(line[i]);
    seg.push_back(line[i+1]);
    dist = bg::distance(p,seg);
    if(dist < best_dist){
      best_dist = dist;
      best_candidate = i;
    }
  }

  // step 2 : projecting the point on the best candidate
  linestring_t best_segment;
  best_segment.push_back(line[best_candidate]);
  best_segment.push_back(line[best_candidate+1]);

  point_t projpt = project_point_on_segment(p, best_segment);

  return projpt;
}



// *********************************************************************
// Add vertices on a line
// *********************************************************************

/*
 *
 * A function to add some points as vertices on nearest lines (identified by index)
 *
 */

// [[Rcpp::export]]
List add_vertices_lines_cpp(NumericMatrix points, List lines, arma::colvec nearest_lines_idx, float mindist){

  // creating the output line list
  //List new_lines;
  std::vector<NumericMatrix> new_lines;
  // creating an arma matrix and an arma colvec for subsetting
  mat Xmat(points.begin(), points.nrow(), points.ncol(), false);

  nearest_lines_idx = nearest_lines_idx-1;
  // iterating over the lines
  int i,j;
  for(i = 0; i < lines.length() ; i++){
    //Rcout << "Iterating on line " << i << '\n';
    NumericMatrix line = lines(i);
    // finding the points to match on that line
    mat ok_pts = Xmat.rows(find(nearest_lines_idx == i));
    //Rcout << "These points will be projected on it \n\n" << ok_pts << "\n\n";

    // if their is not points to add to the line
    if(ok_pts.n_rows == 0){
      new_lines.push_back(line);
    }else{
      // other wise, we have to add the vertices to the line !
      // creating the line
      linestring_t line_geom = line_from_coords(line);

      // step1 : calculating the cumulative distance for the points on the line
      NumericVector X1 = line(_,0);
      NumericVector Y1 = line(_,1);
      int last_elem = X1.length()-2;
      NumericVector X2 = X1[Rcpp::seq(0, last_elem)];
      X2.push_front(X2(0));
      NumericVector Y2 = Y1[Rcpp::seq(0, last_elem)];
      Y2.push_front(Y2(0));
      NumericVector dists = cumsum(Rcpp::sqrt(Rcpp::pow(X1-X2,2) + Rcpp::pow(Y1-Y2,2)));
      mat line_mat(X1.length(),3);
      line_mat.col(0) = as<arma::vec>(X1);
      line_mat.col(1) = as<arma::vec>(Y1);
      line_mat.col(2) = as<arma::vec>(dists);

      //Rcout << "Here is the matrix of the line \n\n" << line_mat << "\n\n";

      // step2 : create a matrix of point dists with the snapped points
      arma::mat distMat(ok_pts.n_rows, 3);
      int nr = ok_pts.n_rows;
      for(j = 0; j < nr; j++){
        //Rcout << "Iterating on point "<<j<<"\n";
        point_t org_pt(ok_pts(j,0),ok_pts(j,1));
        std::pair<double, point_t> point_dits = project_point_on_Linestring_distance(org_pt, line_geom);
        point_t pt = point_dits.second;
        //distMat(j,0) = getX(pt);
        //distMat(j,1) = getY(pt);
        distMat(j,0) = pt.x();
        distMat(j,1) = pt.y();
        distMat(j,2) = point_dits.first;
      }
      //Rcout << "Here is the matrix of the projected pts \n\n" << distMat << "\n\n";

      // subsetting the matrix for dists < mindist
      mat ok_distMat = distMat.rows(find(distMat.col(2) >= mindist));
      //combining the two matrices
      mat total_mat = join_cols(line_mat,ok_distMat);
      //mat total_mat = line_mat;

      //Rcout << "Here is the total matrix \n\n" << total_mat << "\n\n";

      // ordering with the distances
      uvec order = sort_index(total_mat.col(2));
      total_mat = total_mat.rows(order);

      //Rcout << "Here is the ordered total matrix \n\n" << total_mat << "\n\n";

      // adding to the list of lines
      new_lines.push_back(wrap(total_mat.cols(0,1)));
    }


  }
  return wrap(new_lines);

}

/*
 *
 * A function to add the center of lines as a vertex on these lines
 *
 */

// [[Rcpp::export]]
List add_center_lines_cpp(List lines){

  // creating the final container
  std::vector<NumericMatrix> new_lines;
  point_t p(0,0);

  // start the iterations
  int i,j;
  int ll = lines.length();
  for (i = 0; i < ll ; i++){
    NumericMatrix line = lines(i);

    // determining the location of the new point
    linestring_t line_geom = line_from_coords(line);
    double middle = bg::length(line_geom)/2.0;
    bg::line_interpolate(line_geom, middle, p);

    // inserting the new point in new coords
    double cum_dist = 0;
    double prevX = line(0,0);
    double prevY = line(0,1);
    NumericMatrix new_line(line.nrow()+1,2);
    int ad = 0;

    for(j=0; j < line.nrow(); j++){
      double X = line(j,0);
      double Y = line(j,1);
      cum_dist += sqrt(pow(X-prevX,2) + pow(Y-prevY,2));
      if((cum_dist>middle) & (ad == 0)){
        // new_line(j,0) = getX(p);
        // new_line(j,1) = getY(p);
        new_line(j,0) = p.x();
        new_line(j,1) = p.y();
        ad = 1;
      }
      new_line(j+ad,0) = X;
      new_line(j+ad,1) = Y;
    }

    // saving this new matrix
    new_lines.push_back(new_line);
  }

  List final_list = wrap(new_lines);
  return(final_list);

}


// *********************************************************************
// Split lines at points
// *********************************************************************

/*
 *
 * A function to split some lines at indicated points (nearest line must be identified
 * by their index)
 * We return here a list with the new geometries and the index of duplicated lines
 *
 */
// [[Rcpp::export]]
List split_lines_at_points_cpp(arma::mat Xmat, List lines, arma::colvec nearest_lines_idx, float mindist){

  //List new_lines;
  std::vector<NumericMatrix> new_lines;
  std::vector<int> new_index;
  nearest_lines_idx = nearest_lines_idx-1;

  // iterating over the lines
  int i,j,nr;
  for(i = 0; i < lines.length() ; i++){
    NumericMatrix line = lines(i);
    // finding the points to match on that line
    mat ok_pts = Xmat.rows(find(nearest_lines_idx == i));
    // if their is not points to add to the line
    if(ok_pts.n_rows == 0){
      new_lines.push_back(line);
      new_index.push_back(i+1);
    }else{
      // other wise, we have to add the vertices to the line !
      // creating the line
      linestring_t line_geom = line_from_coords(line);

      // step1 : calculating the cumulative distance for the points on the line
      arma::mat line_mat(line.nrow(),4);
      int ln = line.nrow();
      for (j=0 ; j < ln ; j++){
        line_mat(j,0) = line(j,0);
        line_mat(j,1) = line(j,1);
        line_mat(j,3) = 1;
        if(j == 0){
          line_mat(j,2) = 0;
        }else{
          line_mat(j,2) = sqrt(pow(line(j,0) - line(j-1,0),2) + pow(line(j,1) - line(j-1,1),2));
        }
      }
      line_mat.col(2) = arma::cumsum(line_mat.col(2));
      double line_length = line_mat(line.nrow()-1,2);

      // step2 : create a matrix of point dists with the snapped points
      arma::mat distMat(ok_pts.n_rows, 4);
      nr = ok_pts.n_rows;
      for(j = 0; j < nr; j++){
        point_t org_pt(ok_pts(j,0),ok_pts(j,1));
        std::pair<double, point_t> point_dits = project_point_on_Linestring_distance(org_pt, line_geom);
        // distMat(j,0) = getX(point_dits.second);
        // distMat(j,1) = getY(point_dits.second);
        point_t pt = point_dits.second;
        distMat(j,0) = pt.x();
        distMat(j,1) = pt.y();
        distMat(j,2) = point_dits.first;
      }
      // subsetting the matrix for dists < mindist
      mat ok_distMat = distMat.rows(find(distMat.col(2) >= mindist &&  distMat.col(2)<= line_length - mindist));

      ok_distMat.col(3) = colvec(ok_distMat.n_rows,fill::zeros);//as<arma::vec>(original);
      //combining the two matrices
      mat total_mat = join_cols(line_mat,ok_distMat);
      // ordering with the distances
      uvec order = sort_index(total_mat.col(2));
      total_mat = total_mat.rows(order);

      // and now splitting on the points with the column 3 = 0
      int prev = 0;
      nr = total_mat.n_rows;
      for(j = 0; j < nr ; j++){
        if(total_mat(j,3)==0){
          new_lines.push_back(wrap(total_mat.rows(prev,j).cols(0,1)));
          new_index.push_back(i+1);
          prev = j;
        }
      }
      // if we have one more line
      j= total_mat.n_rows-1;
      if(prev != j){
        new_lines.push_back(wrap(total_mat.rows(prev,j).cols(0,1)));
        new_index.push_back(i+1);
      }
    }
  }
  //new_index = new_index+1;
  List final = List::create(wrap(new_lines),wrap(new_index));
  return final;
}


// *********************************************************************
// Lixelize lines
// *********************************************************************

/*
 *
 * A function to split some lines to obtain lixels of equal lengths (but the extremity)
 *
 */
// [[Rcpp::export]]
List lixelize_lines_cpp(List lines, double lx_length, double mindist){

  //creating the containers (a list for the lines and a vector for the new indices)
  std::vector<NumericMatrix> new_lines;
  std::vector<int> new_index;

  // starting the iterations
  int i,j;
  for(i=0; i<lines.length(); i++){
    NumericMatrix line = lines(i);
    linestring_t line_geom = line_from_coords(line);
    double line_length = bg::length(line_geom);

    // if the length is to short
    if(line_length <= lx_length + mindist){
      new_lines.push_back(line);
      new_index.push_back(i);
    }else{
      // otherwise, we have things to do

      // fist, generating the points to add
      std::vector<double> break_lengths;
      std::vector<int> breaker;
      std::vector<double> new_x;
      std::vector<double> new_y;
      point_t p;
      bool continue_loop = true;
      double actual_dist = 0;
      while(continue_loop){
        if((line_length - actual_dist-lx_length) > mindist){
          actual_dist += lx_length;
          bg::line_interpolate(line_geom,actual_dist,p);
          // new_x.push_back(getX(p));
          // new_y.push_back(getY(p));
          new_x.push_back(p.x());
          new_y.push_back(p.y());
          break_lengths.push_back(actual_dist);
          breaker.push_back(1);
        }else{
          continue_loop = false;
        }
      }
      // now adding the real line vertices
      double prevX = line(0,0);
      double prevY = line(0,1);
      double X, Y;
      double tot_dist = 0;
      for (j = 0; j<line.nrow();j++){
        X = line(j,0);
        Y = line(j,1);
        new_x.push_back(X);
        new_y.push_back(Y);
        breaker.push_back(0);
        tot_dist += sqrt(pow(X - prevX,2) + pow(Y - prevY,2));
        break_lengths.push_back(tot_dist);
        prevX = X;
        prevY = Y;
      }

      // creating a nice arma table to sort it after
      mat coords(new_x.size(), 4);
      coords.col(0) = conv_to<colvec>::from(new_x);
      coords.col(1) = conv_to<colvec>::from(new_y);
      coords.col(2) = conv_to<colvec>::from(break_lengths);
      coords.col(3) = conv_to<colvec>::from(breaker);
      // sorting the matrix
      uvec order = sort_index(coords.col(2));
      coords = coords.rows(order);

      // finding the points on which we must do the cut
      std::vector<int> break_idx;
      break_idx.push_back(0);
      int nrs = coords.n_rows;
      for(j=0; j<nrs;j++){
        if(coords(j,3) == 1){
          break_idx.push_back(j);
        }
      }

      break_idx.push_back(coords.n_rows-1);

      // and FINALLY doing the cut
      nrs = (break_idx.size()-1);
      for(j=0; j<nrs; j++){
        new_lines.push_back(wrap(coords.submat(break_idx.at(j), 0, break_idx.at(j+1), 1 )));
        new_index.push_back(i);
      }
    }

  }
  List final = List::create(wrap(new_lines),wrap(new_index));
  return final;

}



// *********************************************************************
// Points along lines
// *********************************************************************

/*
 *
 * A function to create some points along some lines
 *
 */
// [[Rcpp::export]]
NumericMatrix points_along_lines_cpp(List lines, double dist){

  //creating the containers (a list for the points and a vector for the lines indices)
  std::vector<double> new_X;
  std::vector<double> new_Y;
  std::vector<int> new_index;
  point_t p(0,0);

  // starting the iterations
  int i,j;
  for(i=0; i<lines.length(); i++){
    NumericMatrix line = lines(i);
    linestring_t line_geom = line_from_coords(line);
    double line_length = bg::length(line_geom);

    // determining how many points I will create
    std::vector<double> dists = seq_num(0,line_length,dist);

    if(dists.size() == 0){
      dists.push_back(line_length/2.0);
    }

    // creating the points

    for(double & d : dists){
      bg::line_interpolate(line_geom, d, p);
      // new_X.push_back(getX(p));
      // new_Y.push_back(getY(p));
      new_X.push_back(p.x());
      new_Y.push_back(p.y());
      new_index.push_back(i);
    }

    // for (j = 0; j < dists.length(); j++){
    //   double d = dists.get(j);
    //   bg::line_interpolate(line_geom, d, p);
    //   new_X.push_back(p.x());
    //   new_Y.push_back(p.y());
    //   new_index.push_back(i);
    // }

  }

  // creating the final matrix
  NumericMatrix pts_coords(new_X.size(),3);

  NumericVector vec1 = wrap(new_X);
  pts_coords(_,0) = vec1;
  NumericVector vec2 = wrap(new_Y);
  pts_coords(_,1) = vec2;
  NumericVector vec3 = wrap(new_index);
  pts_coords(_,2) = vec3;

  return pts_coords;

}


// *********************************************************************
// Points at the center of lines
// *********************************************************************

/*
 *
 * A function to create some points at the center of some lines
 *
 */
// [[Rcpp::export]]
NumericMatrix points_at_lines_centers_cpp(List lines){

  //creating the containers (a list for the points and a vector for the lines indices)
  std::vector<double> new_X;
  std::vector<double> new_Y;
  point_t p(0,0);

  // starting the iterations
  int i,j;
  for(i=0; i<lines.length(); i++){
    NumericMatrix line = lines(i);
    linestring_t line_geom = line_from_coords(line);
    double line_length = bg::length(line_geom);

    float d = line_length/2.0;
    bg::line_interpolate(line_geom, d, p);
    // new_X.push_back(getX(p));
    // new_Y.push_back(getY(p));
    new_X.push_back(p.x());
    new_Y.push_back(p.y());
  }

  // creating the final matrix
  NumericMatrix pts_coords(new_X.size(),2);

  NumericVector vec1 = wrap(new_X);
  pts_coords(_,0) = vec1;
  NumericVector vec2 = wrap(new_Y);
  pts_coords(_,1) = vec2;

  return pts_coords;

}


