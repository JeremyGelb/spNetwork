#include "spatial_indexing.h"



// --- spatial_index class constructor ----
spatial_index::spatial_index (NumericMatrix x){

  this -> original_coords = x;
  // setting the empty rtree
  MyRtree this_rtree;
  // iterating over the boxes
  int i = 0;
  original_coords = x;
  for(i = 0 ; i <  original_coords.nrow(); i++){
    bbox abox(point_t(original_coords(i,0),original_coords(i,1)),
              point_t(original_coords(i,2),original_coords(i,3)));
    rtree_element el = std::make_pair(abox, i);
    this_rtree.insert(el);
  }
  this -> spindex = this_rtree;
}


// --- setting up spatial index method ----
IntegerVector spatial_index::tree_request(NumericVector reqBbox){
  // defining the vector that will contain the results
  vector_rtree_element returned_values;
  std::vector<int> idx_vector;

  //reqBbox is a vector like c(minX,minY,maxX,maxY)
  bbox region(point_t(reqBbox(0),reqBbox(1)),
              point_t(reqBbox(2),reqBbox(3)));

  //querying the index
  spindex.query(bgi::intersects(region), std::back_inserter(returned_values));

  if(returned_values.size() > 0){
    // if we found boxes, we have to check the distance with geometries
    BOOST_FOREACH(rtree_element const& el, returned_values){
      idx_vector.push_back(el.second);
    }
  }
  // Note that we need to use the R indexing starting at 1
  return Rcpp::IntegerVector(idx_vector.begin(), idx_vector.end()) + 1 ;;
}

