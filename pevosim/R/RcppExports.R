# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

mySample <- function(rates, sum_rates) {
    .Call('pevosim_mySample', PACKAGE = 'pevosim', rates, sum_rates)
}

get_ratesC <- function(betavec, gamma, Ivec, S, overflow) {
    .Call('pevosim_get_ratesC', PACKAGE = 'pevosim', betavec, gamma, Ivec, S, overflow)
}

do_extinctC <- function(state, mut_var, extinct) {
    invisible(.Call('pevosim_do_extinctC', PACKAGE = 'pevosim', state, mut_var, extinct))
}

run_stepC <- function(state, t_tot, t_end, params, maxrate = 1000, debug = TRUE) {
    invisible(.Call('pevosim_run_stepC', PACKAGE = 'pevosim', state, t_tot, t_end, params, maxrate, debug))
}

