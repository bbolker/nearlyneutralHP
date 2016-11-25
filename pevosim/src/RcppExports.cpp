// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/pevosim.h"
#include <Rcpp.h>

using namespace Rcpp;

// mySample
int mySample(NumericVector rates, double sum_rates);
RcppExport SEXP pevosim_mySample(SEXP ratesSEXP, SEXP sum_ratesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rates(ratesSEXP);
    Rcpp::traits::input_parameter< double >::type sum_rates(sum_ratesSEXP);
    rcpp_result_gen = Rcpp::wrap(mySample(rates, sum_rates));
    return rcpp_result_gen;
END_RCPP
}
// get_ratesC
NumericVector get_ratesC(dvec betavec, dvec gamma, ivec Ivec, int S, bool& overflow);
RcppExport SEXP pevosim_get_ratesC(SEXP betavecSEXP, SEXP gammaSEXP, SEXP IvecSEXP, SEXP SSEXP, SEXP overflowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< dvec >::type betavec(betavecSEXP);
    Rcpp::traits::input_parameter< dvec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< ivec >::type Ivec(IvecSEXP);
    Rcpp::traits::input_parameter< int >::type S(SSEXP);
    Rcpp::traits::input_parameter< bool& >::type overflow(overflowSEXP);
    rcpp_result_gen = Rcpp::wrap(get_ratesC(betavec, gamma, Ivec, S, overflow));
    return rcpp_result_gen;
END_RCPP
}
// do_extinctC
void do_extinctC(List state, const String mut_var, int extinct);
RcppExport SEXP pevosim_do_extinctC(SEXP stateSEXP, SEXP mut_varSEXP, SEXP extinctSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type state(stateSEXP);
    Rcpp::traits::input_parameter< const String >::type mut_var(mut_varSEXP);
    Rcpp::traits::input_parameter< int >::type extinct(extinctSEXP);
    do_extinctC(state, mut_var, extinct);
    return R_NilValue;
END_RCPP
}
// run_stepC
void run_stepC(List state, double t_tot, double t_end, List params, double maxrate, bool debug);
RcppExport SEXP pevosim_run_stepC(SEXP stateSEXP, SEXP t_totSEXP, SEXP t_endSEXP, SEXP paramsSEXP, SEXP maxrateSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type state(stateSEXP);
    Rcpp::traits::input_parameter< double >::type t_tot(t_totSEXP);
    Rcpp::traits::input_parameter< double >::type t_end(t_endSEXP);
    Rcpp::traits::input_parameter< List >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< double >::type maxrate(maxrateSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    run_stepC(state, t_tot, t_end, params, maxrate, debug);
    return R_NilValue;
END_RCPP
}
