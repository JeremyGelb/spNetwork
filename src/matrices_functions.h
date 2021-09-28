#ifndef matrices_functions
#define matrices_functions

IntegerMatrix make_matrix(DataFrame df, List neighbour_list);

arma::sp_mat make_matrix_sparse(DataFrame df, List neighbour_list);

arma::sp_mat make_edge_weight_sparse(DataFrame df, List neighbour_list);

#endif
