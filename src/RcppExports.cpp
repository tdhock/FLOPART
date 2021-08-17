// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// get_label_code
Rcpp::IntegerVector get_label_code();
RcppExport SEXP _FLOPART_get_label_code() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(get_label_code());
    return rcpp_result_gen;
END_RCPP
}
// FLOPART_interface
Rcpp::List FLOPART_interface(const Rcpp::IntegerVector data_vec, const Rcpp::NumericVector weight_vec, const double penalty, const Rcpp::IntegerVector label_type_vec, const Rcpp::IntegerVector label_start_vec, const Rcpp::IntegerVector label_end_vec);
RcppExport SEXP _FLOPART_FLOPART_interface(SEXP data_vecSEXP, SEXP weight_vecSEXP, SEXP penaltySEXP, SEXP label_type_vecSEXP, SEXP label_start_vecSEXP, SEXP label_end_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type data_vec(data_vecSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type weight_vec(weight_vecSEXP);
    Rcpp::traits::input_parameter< const double >::type penalty(penaltySEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type label_type_vec(label_type_vecSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type label_start_vec(label_start_vecSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type label_end_vec(label_end_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(FLOPART_interface(data_vec, weight_vec, penalty, label_type_vec, label_start_vec, label_end_vec));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FLOPART_get_label_code", (DL_FUNC) &_FLOPART_get_label_code, 0},
    {"_FLOPART_FLOPART_interface", (DL_FUNC) &_FLOPART_FLOPART_interface, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_FLOPART(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}