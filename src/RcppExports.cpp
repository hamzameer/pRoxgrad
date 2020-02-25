// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// one
int one();
RcppExport SEXP _pRoxgrad_one() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(one());
    return rcpp_result_gen;
END_RCPP
}
// signC
int signC(int x);
RcppExport SEXP _pRoxgrad_signC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(signC(x));
    return rcpp_result_gen;
END_RCPP
}
// sumC
double sumC(NumericVector x);
RcppExport SEXP _pRoxgrad_sumC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(sumC(x));
    return rcpp_result_gen;
END_RCPP
}
// meanC
double meanC(NumericVector x);
RcppExport SEXP _pRoxgrad_meanC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(meanC(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _pRoxgrad_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// soft_thresholding
NumericVector soft_thresholding(NumericVector v, double lambdaL);
RcppExport SEXP _pRoxgrad_soft_thresholding(SEXP vSEXP, SEXP lambdaLSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type lambdaL(lambdaLSEXP);
    rcpp_result_gen = Rcpp::wrap(soft_thresholding(v, lambdaL));
    return rcpp_result_gen;
END_RCPP
}
// timesTwo
NumericVector timesTwo(NumericVector x);
RcppExport SEXP _pRoxgrad_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(timesTwo(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pRoxgrad_one", (DL_FUNC) &_pRoxgrad_one, 0},
    {"_pRoxgrad_signC", (DL_FUNC) &_pRoxgrad_signC, 1},
    {"_pRoxgrad_sumC", (DL_FUNC) &_pRoxgrad_sumC, 1},
    {"_pRoxgrad_meanC", (DL_FUNC) &_pRoxgrad_meanC, 1},
    {"_pRoxgrad_rcpp_hello_world", (DL_FUNC) &_pRoxgrad_rcpp_hello_world, 0},
    {"_pRoxgrad_soft_thresholding", (DL_FUNC) &_pRoxgrad_soft_thresholding, 2},
    {"_pRoxgrad_timesTwo", (DL_FUNC) &_pRoxgrad_timesTwo, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_pRoxgrad(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
