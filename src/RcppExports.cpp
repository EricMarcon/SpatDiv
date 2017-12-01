// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// parallelCountNbd
NumericVector parallelCountNbd(NumericVector r, IntegerVector NbSpecies, NumericVector x, NumericVector y, IntegerVector Type, NumericVector Weight);
RcppExport SEXP _SpatDiv_parallelCountNbd(SEXP rSEXP, SEXP NbSpeciesSEXP, SEXP xSEXP, SEXP ySEXP, SEXP TypeSEXP, SEXP WeightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type NbSpecies(NbSpeciesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Type(TypeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Weight(WeightSEXP);
    rcpp_result_gen = Rcpp::wrap(parallelCountNbd(r, NbSpecies, x, y, Type, Weight));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SpatDiv_parallelCountNbd", (DL_FUNC) &_SpatDiv_parallelCountNbd, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_SpatDiv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
