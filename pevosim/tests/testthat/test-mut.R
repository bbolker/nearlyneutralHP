library(testthat)
library(pevosim)

## multimutation tests
## devtools::load_all("~/projects/nearlyneutralHP/pevosim")
context("new mutation")
test_that("basics 2", {
    set.seed(101)
    S0 <- list(traits=list(beta=list(constr=1,unconstr=0),
                           gamma=list(constr=1,unconstr=0)),
               Ivec=1,
               S=10)
    expect_true(check_state(S0))
    set.seed(101)
    S1 <- do_mut3(S0,mut_var="beta",mut_link=make.link("log"),strain=1,
                  mut_mean=0,mut_sd=1)
    expect_true(check_state(S1))
    S2 <- do_mut3(S1,mut_var="both",mut_link=make.link("log"),strain=1,
                  mut_mean=0,mut_sd=1)
    expect_true(check_state(S2))
    S3 <- do_mut3(S2,mut_var="both",mut_link=make.link("log"),strain=1,
                  mut_mean=0,mut_sd=1)
    expect_true(check_state(S3))
    ## do more testing of what actually happens

    if (FALSE) {
        set.seed(101)
        params <- c(g1=1,g2=2)
        mm <- make_bilink("powerlaw",params,gamma=2,beta=0.5)
        with(as.list(params),curve(g1*x^(1/g2),from=0,to=5,ylim=c(0,6),n=501))
        abline(v=2,lty=2)
        abline(h=0.5,lty=2)
        with(as.list(params),(0.5/g1)^g2)  ## indeed, 16
        points(1,1)  ## this is *on* the curve - will cause problems when mutating
## do we want 
## oops.  exceed power-law limit!
## how are we scaling beta?? Not dividing by N in computing rates
##  so beta should be small (R0=beta*N/gamma)
do_mut3(S0,mut_var="gamma",mut_link=mm,strain=1,
        mut_mean=0,mut_sd=1)
## gamma mut not working??
mm$gamma$linkinv(0)
environment(mm$gamma$linkinv)$minval
curve(mm$gamma$linkinv(x),from=-5,to=5)
curve(mm$beta$linkinv(x),from=-5,to=5)
}

    })
## run_sim(tradeoff_fun="none")

## R0 = 2
## I = 500
## S = 500
## beta0 = 2*0.2/1000 = 4e-4
## beta0*N/gamma = 4e-4*1000/0.2 = 0.4/0.2 = 2 = R0 (check)

##
## we should really be using cloglog link from beta to p?
