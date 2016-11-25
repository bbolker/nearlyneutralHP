library(pevosim)
library(testthat)

context("overflow")
test_that("overflow",{
    pars <- list(mu = 0.0237137370566166, mut_mean = -1.15625, 
                 mut_sd = 1.65481709994318, N = 24,
                 gamma0 = 0.00316227766016838)
    argList <- c(list(nt=1e7,R0_init=1e300,seed=101,
                      mut_var="beta",
                      discrete=FALSE,
                      useCpp=TRUE,
                      progress=TRUE,
                  debug=FALSE),
                 as.list(pars))
    expect_error(do.call(run_sim,argList),"Rate overflow")
})
