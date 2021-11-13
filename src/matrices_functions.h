#ifndef matrices_functions
#define matrices_functions

std::vector<double> seq_num(double start, double end, double step);

int get_first_index(NumericVector& v1, double x);

std::vector<int> get_all_indeces(NumericVector& v1, double x);

IntegerMatrix make_matrix(DataFrame df, List neighbour_list);

arma::sp_mat make_matrix_sparse(DataFrame df, List neighbour_list);

arma::sp_mat make_edge_weight_sparse(DataFrame df, List neighbour_list);

#endif
