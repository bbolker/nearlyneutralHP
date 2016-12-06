## miscellanous messing-around and tests
library(testthat)
library(pevosim)
make_testdata <- FALSE

if (FALSE) {
    path <- c("/media/sf_Documents/projects/nearlyneutralHP/pevosim")
    ## source(file.path(path,"R/funs.R"))
    ## Rcpp::sourceCpp(file.path(path,"src/funs.cpp"))
    ## devtools::load_all(path)
}

if (make_testdata) {
    res1S <- run_sim(nt=1e3,rptfreq=100,seed=101)
    res2S <- run_sim(nt=1e3,rptfreq=100,seed=101,
                    mod_params=list(mut_var="gamma"))
    res3S <- run_sim(nt=1e3,rptfreq=100,seed=101,
                  mod_params=list(mut_var="gamma",
                                  mut_link=list(beta=make.link("logit"),
                                                gamma=multlogit(0.1,0.5))))
    res4S <- run_sim(nt=100,rptfreq=10,seed=101,
                    discrete=FALSE)
    save(list=ls(pattern="res[0-9]+S"),
         file="../../inst/testdata/testruns_1.rda")
}
load(system.file("testdata","testruns_1.rda",package="pevosim"))
context("discrete models")
test_that("basics", {
    ## basic run
    res1 <- run_sim(nt=1e3,rptfreq=100,seed=101)
    expect_equal(c(res1[10,]),c(res1S[10,]))
    ## mutating gamma
    res2 <- run_sim(nt=1e3,rptfreq=100,seed=101,
                    mod_params=list(mut_var="gamma"))
    expect_equal(c(res2[10,]),c(res2S[10,]))
    ## gamma constrained
    res3 <- run_sim(nt=1e3,rptfreq=100,seed=101,
                  mod_params=list(mut_var="gamma",
                                  mut_link=list(beta=make.link("logit"),
                                                gamma=multlogit(0.1,0.5))))
    expect_equal(c(res3[10,]),c(res3S[10,]))
})

test_that("multlogit", {
    Lfun0 <- make.link("logit")
    Lfun <- multlogit(0,1)
    expect_equal(Lfun$linkinv(as.numeric(-5:5)),Lfun0$linkinv(as.numeric(-5:5)))
})

test_that("continuous", {
    res4 <- run_sim(nt=100,rptfreq=10,seed=101,
                    discrete=FALSE)
    expect_equal(c(res4[10,]),c(res4S[10,]))
    })

## UNIT TESTS, continuous-time models

test_that("single steps", {
    ## single infection
    r1 <- run_sim(nt=1e-8,
                  mod_init=list(R0=10,Ivec=1),
                  rptfreq=1e-8,
                  discrete=FALSE,seed=101,useCpp=FALSE)

    expect_equal(as.data.frame(c(r1)),
                 data.frame(time = 1e-08, S = 998, I = 2,
                            mean_lbeta = -6.21460809842219,
                            sd_lbeta=0,
                            mean_lgamma = -1.6094379124341,
                            sd_lgamma = 0,
                            n_events=1))

    if (FALSE) {
    r2 <- run_sim(nt=1e-8,R0_init=10,Ivec=1,rptfreq=1e-8,discrete=FALSE,seed=101,useCpp=TRUE)
     expect_equal(c(r1),c(r2))
    }
    ## single infection plus mutation
    if (FALSE) {
        r3 <- run_sim(nt=1e-8,R0_init=10,Ivec=1,rptfreq=1e-8,mu=1,discrete=FALSE,seed=101,useCpp=TRUE)
            expect_equal(c(r3),
                 structure(list(time = 1e-08, S = 998, I = 2, mean_lbeta = -6.88334405941115, 
                        sd_lbeta = 0.668735960988958),
          .Names = c("time", "S", "I", "mean_lbeta", "sd_lbeta")))
    }
    r4 <- run_sim(nt=1e-8,
                  mod_inits=list(R0=10,Ivec=1),
                  mod_params=list(mu=1),
                  rptfreq=1e-8,discrete=FALSE,seed=101,useCpp=FALSE)


    expect_equal(as.data.frame(c(r4)),
                 data.frame(time = 1e-08, S = 998, I = 2,
                            mean_lbeta = -6.61306611533102,
                            sd_lbeta=0.39845801690882,
                            mean_lgamma = -1.6094379124341,
                            sd_lgamma = 0,
                            n_events=1))


    ## single recovery
    if (FALSE) {
        expect_message(r5 <- run_sim(nt=1e-8,rptfreq=1e-8,R0_init=0.01,Ivec=1,
                                     discrete=FALSE,seed=101,useCpp=TRUE),"system went extinct")
        expect_equal(c(r5),
                     structure(list(time = 1e-08, S = 1000, I = 0, mean_lbeta = NaN, 
                                    sd_lbeta = NaN), .Names = c("time", "S", "I", "mean_lbeta", "sd_lbeta")))
    }
    expect_message(r6 <- run_sim(nt=1e-8,rptfreq=1e-8,
                            mod_init=list(R0=0.01,Ivec=1),
                            discrete=FALSE,seed=101,useCpp=FALSE),
                  "system went extinct")
    expect_equal(as.data.frame(c(r6)),
                 data.frame(time = 1e-08, S = 1000, I = 0, mean_lbeta = NaN,
                            sd_lbeta = NaN,
                            mean_lgamma=NaN,
                            sd_lgamma = NaN,
                            n_events=1))
})
