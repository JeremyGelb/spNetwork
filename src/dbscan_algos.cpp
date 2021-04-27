// // [[Rcpp::depends(BH)]]
//
// #include "spNetwork.h"
// #include "base_kernel_funtions.h"
// #include <boost/config.hpp>
// #include <boost/graph/adjacency_list.hpp>
// #include <boost/lexical_cast.hpp>
//
// using namespace boost;
//
// //quick example https://www.boost.org/doc/libs/1_75_0/libs/graph/example/family_tree.cpp
// // for RCPP help : https://teuder.github.io/rcpp4everyone_en/080_vector.html
//
// // create a typedef for the Graph type
// typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
// struct Vertex { int id; int isevent; };
//
// typedef adjacency_list<vecS, vecS, undirectedS , Vertex, EdgeWeightProperty> Graph;
// typedef boost::graph_traits<Graph>::edge_iterator edge_iterator;
//
//
//
// // Definition de certains types
// typedef std::pair<int, int> Edge;
// //Type utilise pour ce que retourne build_cpp_graph
// // OID ==> Graph vertex descriptor
// typedef std::map<int, boost::graph_traits<Graph>::vertex_descriptor> VertexIndex;
// typedef std::pair<Graph, VertexIndex> GraphMaper;
//
// //function to create an adjency list form a LinesDataFrame
// //must also return a map of the vertex (oid => vertex)
// GraphMaper build_cpp_graph(DataFrame vertices, DataFrame linesdf){
//
//   Graph G;
//   VertexIndex vertex_map;
//   // extraire les attributs du dataframe de vertices
//   NumericVector vertidx = vertices["oid"];
//   NumericVector isevent = vertices["isevent"];
//
//   // ajoutons chaque vertex dans le graph
//   int n = vertidx.length();
//   for(int j=0; j < n; ++j){
//     vertex_map.emplace(vertidx[j],add_vertex(Vertex{vertidx[j], isevent[j]}, G));
//   }
//
//   // extraire les attributs du dataframe de lignes
//   NumericVector start_ids = linesdf["start_oid"];
//   NumericVector end_ids = linesdf["end_oid"];
//   NumericVector line_weight = linesdf["weight"];
//
//   // // ecrivons chacune des lignes dans le graph
//
//   n = start_ids.length();
//   for(int j=0; j < n; ++j){
//     add_edge(start_ids[j], end_ids[j], line_weight[j], G);
//   }
//   return GraphMaper (G, vertex_map);
// };
//
//
// //adjacent_vertices(v, g);
// NumericVector LSPD_alog(GraphMaper GM, int V, double r, NumericVector events){
//
//   // recuperation du graph et de la map de vertex
//   Graph G = GM.first;
//   VertexIndex vertex_idx = GM.second;
//
//   // definition de certain type pour mes manipulations (WTF AM I DOING ?)
//   typename graph_traits < Graph >::out_edge_iterator ei, ei_end;
//
//   // extracting a the edgeweightmap (WTF is this shit ?)
//   boost::property_map<Graph, boost::edge_weight_t>::type EdgeWeightMap = get(boost::edge_weight_t(), G);
//
//   // Initialisation des valeurs
//
//   // CDV = inf (-1) values for all other vertices
//   std::map<int, double> CDV_values;
//   Graph::vertex_iterator v, vend;
//   int inf = std::numeric_limits<int>::max();
//   for (boost::tie(v, vend) = vertices(G); v != vend; ++v) {
//     int oid = G[*v].id;
//     if(oid != V){
//       CDV_values[G[*v].id] = inf;
//     }else{
//       CDV_values[G[*v].id] = 0;
//     }
//
//   }
//
//   NumericVector Neps;
//   queue <int> Q;
//   Q.push(V);
//
//   // lancement des iterations
//   while(Q.empty()==FALSE){
//     int p = Q.front();
//     Q.pop();
//     // on test (1) si la vertex p est bien un event
//     // on test (2) si on a pas encore inserer cette valeur dans Neps
//     int test1 = G[vertex_idx[p]].isevent;
//     int test2 = get_first_index(Neps, p);
//     if(test1 == 1  and test2 < 0){
//       Neps.push_back(p);
//     }
//     boost::graph_traits <Graph>::adjacency_iterator ai, a_end;
//     // iterating on the edges connected to p
//     for (boost::tie(ei, ei_end) = out_edges(p, G); ei != ei_end; ++ei) {
//       //ei is an edge_descriptor
//       auto source = boost::source ( *ei, G );
//       //getting the weight of the edge
//       double w = EdgeWeightMap[*ei];
//
//       if CDV_values[p]
//
//     }
//   }
//
//   return Neps;
// }
//
//
// // TESTING HERE
// // [[Rcpp::export]]
// void easytest(DataFrame linesdf, DataFrame verticesdf, int v){
//   Rcout << "Building the Graph ! \n";
//   GraphMaper GM = build_cpp_graph(verticesdf, linesdf);
//   Graph G = GM.first;
//   VertexIndex vertidxmap = GM.second;
//   Rcout << G[vertidxmap[v]].isevent << "\n";
//
//   // typename graph_traits < Graph >::out_edge_iterator ei, ei_end;
//   // boost::property_map<Graph, boost::edge_weight_t>::type EdgeWeightMap = get(boost::edge_weight_t(), G);
//   //
//   // for (boost::tie(ei, ei_end) = out_edges(vert, G); ei != ei_end; ++ei) {
//   //          //ei is an edge_descriptor
//   //          auto source = boost::source ( *ei, G );
//   //          auto target = boost::target ( *ei, G );
//   //          //getting the weight of the edge
//   //          double w = EdgeWeightMap[*ei];
//   //          Rcout << "There is an edge from " << source <<  " to " << target << " with length " << w << "\n" ;
//   // }
//
//   };
