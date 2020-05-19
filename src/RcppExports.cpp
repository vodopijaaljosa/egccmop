// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// s_ring
NumericVector s_ring(NumericVector weights, int no_elevators, double prob, double alpha, double beta, int no_cycles);
RcppExport SEXP _egccmop_s_ring(SEXP weightsSEXP, SEXP no_elevatorsSEXP, SEXP probSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP no_cyclesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type no_elevators(no_elevatorsSEXP);
    Rcpp::traits::input_parameter< double >::type prob(probSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type no_cycles(no_cyclesSEXP);
    rcpp_result_gen = Rcpp::wrap(s_ring(weights, no_elevators, prob, alpha, beta, no_cycles));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_egccmop_s_ring", (DL_FUNC) &_egccmop_s_ring, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_egccmop(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
