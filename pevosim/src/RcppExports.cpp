// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

typedef std::vector<double> dvec;
typedef std::vector<int> ivec;

// get_rates
NumericVector get_rates(dvec betavec, dvec gamma, ivec Ivec, int S);
RcppExport SEXP pevosim_get_rates(SEXP betavecSEXP, SEXP gammaSEXP, SEXP IvecSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< dvec >::type betavec(betavecSEXP);
    Rcpp::traits::input_parameter< dvec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< ivec >::type Ivec(IvecSEXP);
    Rcpp::traits::input_parameter< int >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(get_rates(betavec, gamma, Ivec, S));
    return rcpp_result_gen;
END_RCPP
}
// do_mut
void do_mut(List state, const String mut_var, dvec orig_trait, double mut_mean, double mut_sd);
RcppExport SEXP pevosim_do_mut(SEXP stateSEXP, SEXP mut_varSEXP, SEXP orig_traitSEXP, SEXP mut_meanSEXP, SEXP mut_sdSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type state(stateSEXP);
    Rcpp::traits::input_parameter< const String >::type mut_var(mut_varSEXP);
    Rcpp::traits::input_parameter< dvec >::type orig_trait(orig_traitSEXP);
    Rcpp::traits::input_parameter< double >::type mut_mean(mut_meanSEXP);
    Rcpp::traits::input_parameter< double >::type mut_sd(mut_sdSEXP);
    do_mut(state, mut_var, orig_trait, mut_mean, mut_sd);
    return R_NilValue;
END_RCPP
}
// do_extinct
void do_extinct(List state, const String mut_var, int extinct);
RcppExport SEXP pevosim_do_extinct(SEXP stateSEXP, SEXP mut_varSEXP, SEXP extinctSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type state(stateSEXP);
    Rcpp::traits::input_parameter< const String >::type mut_var(mut_varSEXP);
    Rcpp::traits::input_parameter< int >::type extinct(extinctSEXP);
    do_extinct(state, mut_var, extinct);
    return R_NilValue;
END_RCPP
}
// run_stepC
void run_stepC(List state, double t_tot, double t_end, List params, bool debug);
RcppExport SEXP pevosim_run_stepC(SEXP stateSEXP, SEXP t_totSEXP, SEXP t_endSEXP, SEXP paramsSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type state(stateSEXP);
    Rcpp::traits::input_parameter< double >::type t_tot(t_totSEXP);
    Rcpp::traits::input_parameter< double >::type t_end(t_endSEXP);
    Rcpp::traits::input_parameter< List >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    run_stepC(state, t_tot, t_end, params, debug);
    return R_NilValue;
END_RCPP
}
