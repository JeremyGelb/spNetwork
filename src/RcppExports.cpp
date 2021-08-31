// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// get_loo_values_continuous
DataFrame get_loo_values_continuous(List neighbour_list, NumericVector samples, NumericVector sweights, NumericVector events, NumericVector weights, NumericVector bws, std::string kernel_name, DataFrame line_list, int max_depth);
RcppExport SEXP _spNetwork_get_loo_values_continuous(SEXP neighbour_listSEXP, SEXP samplesSEXP, SEXP sweightsSEXP, SEXP eventsSEXP, SEXP weightsSEXP, SEXP bwsSEXP, SEXP kernel_nameSEXP, SEXP line_listSEXP, SEXP max_depthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type neighbour_list(neighbour_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sweights(sweightsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bws(bwsSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel_name(kernel_nameSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type line_list(line_listSEXP);
    Rcpp::traits::input_parameter< int >::type max_depth(max_depthSEXP);
    rcpp_result_gen = Rcpp::wrap(get_loo_values_continuous(neighbour_list, samples, sweights, events, weights, bws, kernel_name, line_list, max_depth));
    return rcpp_result_gen;
END_RCPP
}
// get_loo_values_discontinuous
NumericVector get_loo_values_discontinuous(List neighbour_list, NumericVector samples, NumericVector sweights, NumericVector events, NumericVector weights, NumericVector bws, std::string kernel_name, DataFrame line_list, int max_depth);
RcppExport SEXP _spNetwork_get_loo_values_discontinuous(SEXP neighbour_listSEXP, SEXP samplesSEXP, SEXP sweightsSEXP, SEXP eventsSEXP, SEXP weightsSEXP, SEXP bwsSEXP, SEXP kernel_nameSEXP, SEXP line_listSEXP, SEXP max_depthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type neighbour_list(neighbour_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sweights(sweightsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bws(bwsSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel_name(kernel_nameSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type line_list(line_listSEXP);
    Rcpp::traits::input_parameter< int >::type max_depth(max_depthSEXP);
    rcpp_result_gen = Rcpp::wrap(get_loo_values_discontinuous(neighbour_list, samples, sweights, events, weights, bws, kernel_name, line_list, max_depth));
    return rcpp_result_gen;
END_RCPP
}
// get_loo_values_simple
NumericVector get_loo_values_simple(List neighbour_list, NumericVector samples, NumericVector sweights, NumericVector events, NumericVector weights, NumericVector bws, std::string kernel_name, DataFrame line_list, int max_depth);
RcppExport SEXP _spNetwork_get_loo_values_simple(SEXP neighbour_listSEXP, SEXP samplesSEXP, SEXP sweightsSEXP, SEXP eventsSEXP, SEXP weightsSEXP, SEXP bwsSEXP, SEXP kernel_nameSEXP, SEXP line_listSEXP, SEXP max_depthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type neighbour_list(neighbour_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sweights(sweightsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bws(bwsSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel_name(kernel_nameSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type line_list(line_listSEXP);
    Rcpp::traits::input_parameter< int >::type max_depth(max_depthSEXP);
    rcpp_result_gen = Rcpp::wrap(get_loo_values_simple(neighbour_list, samples, sweights, events, weights, bws, kernel_name, line_list, max_depth));
    return rcpp_result_gen;
END_RCPP
}
// timesTwo
NumericVector timesTwo(NumericVector x);
RcppExport SEXP _spNetwork_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(timesTwo(x));
    return rcpp_result_gen;
END_RCPP
}
// continuous_nkde_cpp_arma_sparse
DataFrame continuous_nkde_cpp_arma_sparse(List neighbour_list, NumericVector events, NumericVector weights, DataFrame samples, NumericVector bws, std::string kernel_name, DataFrame nodes, DataFrame line_list, int max_depth, bool verbose);
RcppExport SEXP _spNetwork_continuous_nkde_cpp_arma_sparse(SEXP neighbour_listSEXP, SEXP eventsSEXP, SEXP weightsSEXP, SEXP samplesSEXP, SEXP bwsSEXP, SEXP kernel_nameSEXP, SEXP nodesSEXP, SEXP line_listSEXP, SEXP max_depthSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type neighbour_list(neighbour_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bws(bwsSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel_name(kernel_nameSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type line_list(line_listSEXP);
    Rcpp::traits::input_parameter< int >::type max_depth(max_depthSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(continuous_nkde_cpp_arma_sparse(neighbour_list, events, weights, samples, bws, kernel_name, nodes, line_list, max_depth, verbose));
    return rcpp_result_gen;
END_RCPP
}
// continuous_nkde_cpp_arma
DataFrame continuous_nkde_cpp_arma(List neighbour_list, NumericVector events, NumericVector weights, DataFrame samples, NumericVector bws, std::string kernel_name, DataFrame nodes, DataFrame line_list, int max_depth, bool verbose);
RcppExport SEXP _spNetwork_continuous_nkde_cpp_arma(SEXP neighbour_listSEXP, SEXP eventsSEXP, SEXP weightsSEXP, SEXP samplesSEXP, SEXP bwsSEXP, SEXP kernel_nameSEXP, SEXP nodesSEXP, SEXP line_listSEXP, SEXP max_depthSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type neighbour_list(neighbour_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bws(bwsSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel_name(kernel_nameSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type line_list(line_listSEXP);
    Rcpp::traits::input_parameter< int >::type max_depth(max_depthSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(continuous_nkde_cpp_arma(neighbour_list, events, weights, samples, bws, kernel_name, nodes, line_list, max_depth, verbose));
    return rcpp_result_gen;
END_RCPP
}
// discontinuous_nkde_cpp_arma_sparse
DataFrame discontinuous_nkde_cpp_arma_sparse(List neighbour_list, NumericVector events, NumericVector weights, DataFrame samples, NumericVector bws, std::string kernel_name, DataFrame nodes, DataFrame line_list, int max_depth, bool verbose);
RcppExport SEXP _spNetwork_discontinuous_nkde_cpp_arma_sparse(SEXP neighbour_listSEXP, SEXP eventsSEXP, SEXP weightsSEXP, SEXP samplesSEXP, SEXP bwsSEXP, SEXP kernel_nameSEXP, SEXP nodesSEXP, SEXP line_listSEXP, SEXP max_depthSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type neighbour_list(neighbour_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bws(bwsSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel_name(kernel_nameSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type line_list(line_listSEXP);
    Rcpp::traits::input_parameter< int >::type max_depth(max_depthSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(discontinuous_nkde_cpp_arma_sparse(neighbour_list, events, weights, samples, bws, kernel_name, nodes, line_list, max_depth, verbose));
    return rcpp_result_gen;
END_RCPP
}
// discontinuous_nkde_cpp_arma
DataFrame discontinuous_nkde_cpp_arma(List neighbour_list, NumericVector events, NumericVector weights, DataFrame samples, NumericVector bws, std::string kernel_name, DataFrame nodes, DataFrame line_list, int max_depth, bool verbose);
RcppExport SEXP _spNetwork_discontinuous_nkde_cpp_arma(SEXP neighbour_listSEXP, SEXP eventsSEXP, SEXP weightsSEXP, SEXP samplesSEXP, SEXP bwsSEXP, SEXP kernel_nameSEXP, SEXP nodesSEXP, SEXP line_listSEXP, SEXP max_depthSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type neighbour_list(neighbour_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bws(bwsSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel_name(kernel_nameSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type line_list(line_listSEXP);
    Rcpp::traits::input_parameter< int >::type max_depth(max_depthSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(discontinuous_nkde_cpp_arma(neighbour_list, events, weights, samples, bws, kernel_name, nodes, line_list, max_depth, verbose));
    return rcpp_result_gen;
END_RCPP
}
// corrfactor_discontinuous_sparse
List corrfactor_discontinuous_sparse(List neighbour_list, NumericVector events, DataFrame line_list, NumericVector bws, int max_depth);
RcppExport SEXP _spNetwork_corrfactor_discontinuous_sparse(SEXP neighbour_listSEXP, SEXP eventsSEXP, SEXP line_listSEXP, SEXP bwsSEXP, SEXP max_depthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type neighbour_list(neighbour_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type line_list(line_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bws(bwsSEXP);
    Rcpp::traits::input_parameter< int >::type max_depth(max_depthSEXP);
    rcpp_result_gen = Rcpp::wrap(corrfactor_discontinuous_sparse(neighbour_list, events, line_list, bws, max_depth));
    return rcpp_result_gen;
END_RCPP
}
// corrfactor_discontinuous
List corrfactor_discontinuous(List neighbour_list, NumericVector events, DataFrame line_list, NumericVector bws, int max_depth);
RcppExport SEXP _spNetwork_corrfactor_discontinuous(SEXP neighbour_listSEXP, SEXP eventsSEXP, SEXP line_listSEXP, SEXP bwsSEXP, SEXP max_depthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type neighbour_list(neighbour_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type line_list(line_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bws(bwsSEXP);
    Rcpp::traits::input_parameter< int >::type max_depth(max_depthSEXP);
    rcpp_result_gen = Rcpp::wrap(corrfactor_discontinuous(neighbour_list, events, line_list, bws, max_depth));
    return rcpp_result_gen;
END_RCPP
}
// corrfactor_continuous_sparse
List corrfactor_continuous_sparse(List neighbour_list, NumericVector events, DataFrame line_list, NumericVector bws, int max_depth);
RcppExport SEXP _spNetwork_corrfactor_continuous_sparse(SEXP neighbour_listSEXP, SEXP eventsSEXP, SEXP line_listSEXP, SEXP bwsSEXP, SEXP max_depthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type neighbour_list(neighbour_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type line_list(line_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bws(bwsSEXP);
    Rcpp::traits::input_parameter< int >::type max_depth(max_depthSEXP);
    rcpp_result_gen = Rcpp::wrap(corrfactor_continuous_sparse(neighbour_list, events, line_list, bws, max_depth));
    return rcpp_result_gen;
END_RCPP
}
// corrfactor_continuous
List corrfactor_continuous(List neighbour_list, NumericVector events, DataFrame line_list, NumericVector bws, int max_depth);
RcppExport SEXP _spNetwork_corrfactor_continuous(SEXP neighbour_listSEXP, SEXP eventsSEXP, SEXP line_listSEXP, SEXP bwsSEXP, SEXP max_depthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type neighbour_list(neighbour_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type line_list(line_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bws(bwsSEXP);
    Rcpp::traits::input_parameter< int >::type max_depth(max_depthSEXP);
    rcpp_result_gen = Rcpp::wrap(corrfactor_continuous(neighbour_list, events, line_list, bws, max_depth));
    return rcpp_result_gen;
END_RCPP
}
// cut_lines_at_distances_cpp
List cut_lines_at_distances_cpp(List lines, NumericVector dists);
RcppExport SEXP _spNetwork_cut_lines_at_distances_cpp(SEXP linesSEXP, SEXP distsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type lines(linesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dists(distsSEXP);
    rcpp_result_gen = Rcpp::wrap(cut_lines_at_distances_cpp(lines, dists));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spNetwork_get_loo_values_continuous", (DL_FUNC) &_spNetwork_get_loo_values_continuous, 9},
    {"_spNetwork_get_loo_values_discontinuous", (DL_FUNC) &_spNetwork_get_loo_values_discontinuous, 9},
    {"_spNetwork_get_loo_values_simple", (DL_FUNC) &_spNetwork_get_loo_values_simple, 9},
    {"_spNetwork_timesTwo", (DL_FUNC) &_spNetwork_timesTwo, 1},
    {"_spNetwork_continuous_nkde_cpp_arma_sparse", (DL_FUNC) &_spNetwork_continuous_nkde_cpp_arma_sparse, 10},
    {"_spNetwork_continuous_nkde_cpp_arma", (DL_FUNC) &_spNetwork_continuous_nkde_cpp_arma, 10},
    {"_spNetwork_discontinuous_nkde_cpp_arma_sparse", (DL_FUNC) &_spNetwork_discontinuous_nkde_cpp_arma_sparse, 10},
    {"_spNetwork_discontinuous_nkde_cpp_arma", (DL_FUNC) &_spNetwork_discontinuous_nkde_cpp_arma, 10},
    {"_spNetwork_corrfactor_discontinuous_sparse", (DL_FUNC) &_spNetwork_corrfactor_discontinuous_sparse, 5},
    {"_spNetwork_corrfactor_discontinuous", (DL_FUNC) &_spNetwork_corrfactor_discontinuous, 5},
    {"_spNetwork_corrfactor_continuous_sparse", (DL_FUNC) &_spNetwork_corrfactor_continuous_sparse, 5},
    {"_spNetwork_corrfactor_continuous", (DL_FUNC) &_spNetwork_corrfactor_continuous, 5},
    {"_spNetwork_cut_lines_at_distances_cpp", (DL_FUNC) &_spNetwork_cut_lines_at_distances_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_spNetwork(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
