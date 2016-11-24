## miscellanous messing-around and tests
library(testthat)
library(pevosim)

if (FALSE) {
    path <- c("/media/sf_Documents/projects/nearlyneutralHP/pevosim")
    source(file.path(path,"R/funs.R"))
    Rcpp::sourceCpp(file.path(path,"src/funs.cpp"))
}

context("discrete models")
test_that("basics", {
    ## basic run
    res1 <- run_sim(nt=1e4,rptfreq=100,seed=101)
    expect_equal(c(res1[100,]),
  structure(list(time = 10000, S = 160, I = 840, mean_lbeta = 1.50094828644975, 
  sd_lbeta = 0.11040310840802), .Names = c("time", "S", "I", "mean_lbeta", "sd_lbeta")))
    ## mutating gamma
    res2 <- run_sim(nt=1e4,rptfreq=100,seed=101,mut_var="gamma")
    expect_equal(c(res2[100,]),
 structure(list(time = 10000, S = 0, I = 1000, mean_lgamma = -8.14914364919813, 
 sd_lgamma = 0.547470943260739), .Names = c("time", "S", "I", "mean_lgamma", "sd_lgamma")))

  res3 <- run_sim(nt=1e4,rptfreq=100,seed=101,mut_var="gamma",
                mut_link=multlogit(0.1,0.5))
  expect_equal(c(res3[100,]),
 structure(list(time = 10000, S = 332, I = 668, mean_lgamma = -14.202599029542, 
 sd_lgamma = 1.03017299599263), .Names = c("time", "S", "I", 
"mean_lgamma", "sd_lgamma")))

}

test_that("multlogit", {
    Lfun0 <- make.link("logit")
    Lfun <- multlogit(0,1)
    all.equal(Lfun$linkinv(as.numeric(-5:5)),Lfun0$linkinv(as.numeric(-5:5)))
})

test_that("continuous", {
    res4 <- run_sim(nt=1e3,rptfreq=10,seed=101,discrete=FALSE)
    expect_equal(c(res4[100,]),
 structure(list(time = 1000, S = 137, I = 863, mean_lbeta = -6.36886258791232, 
  sd_lbeta = 0.159849030801711), .Names = c("time", "S", "I", 
"mean_lbeta", "sd_lbeta")))
}

if (FALSE) {
    t1A <- system.time(res1A <- run_sim(nt=1e4,rptfreq=10,seed=101,progress=TRUE))
    t1B <- system.time(res1B <- run_sim(nt=1e4,rptfreq=10,seed=101,discrete=FALSE,progress=TRUE))
    t1C <- system.time(res1C <- run_sim(nt=1e4,rptfreq=10,seed=101,discrete=FALSE,
                                 useCpp=TRUE,progress=TRUE))
    t1D <- system.time(res1D <- run_sim(nt=1e4,rptfreq=10,seed=101,discrete=FALSE,
                                        mut_link=make.link("logit"),progress=TRUE))

    save(list=ls(pattern=("(res|t)1[A-D]")),file="../simdata/simdisc.rda")
    load("../simdata/simdisc.rda")
    library(ggplot2); theme_set(theme_bw())
    library(tidyr)
    library(dplyr)
    rL <- lapply(list(res1A,res1B,res1C,res1D),gather,key=var,value=val,-time) %>%
        setNames(c("disc","cont","contCpp","cont_logit")) %>% bind_rows(.id="type")
    ggplot(rL,aes(time,val,colour=type))+geom_line()+
        facet_wrap(~var,scale="free")+
        scale_colour_brewer(palette="Set1")
}

## UNIT TESTS, continuous-time models

test_that("single steps", {
    ## single infection
    r1 <- run_sim(nt=1e-8,R0_init=10,Ivec=1,rptfreq=1e-8,discrete=FALSE,seed=101,useCpp=TRUE)
    expect_equal(c(r1),
                 structure(list(time = 1e-08, S = 998, I = 2, mean_lbeta = -6.21460809842219, 
                                sd_lbeta = 0), .Names = c("time", "S", "I", "mean_lbeta", "sd_lbeta")))
    r2 <- run_sim(nt=1e-8,R0_init=10,Ivec=1,rptfreq=1e-8,discrete=FALSE,seed=101,useCpp=FALSE)
    expect_equal(c(r1),c(r2))
    ## single infection plus mutation
    r3 <- run_sim(nt=1e-8,R0_init=10,Ivec=1,rptfreq=1e-8,mu=1,discrete=FALSE,seed=101,useCpp=TRUE)
    expect_equal(c(r3),
         structure(list(time = 1e-08, S = 998, I = 2, mean_lbeta = -6.88334405941115, 
                        sd_lbeta = 0.668735960988958),
          .Names = c("time", "S", "I", "mean_lbeta", "sd_lbeta")))
    r4 <- run_sim(nt=1e-8,R0_init=10,Ivec=1,rptfreq=1e-8,mu=1,discrete=FALSE,seed=101,useCpp=FALSE)
    expect_equal(c(r4),structure(list(time = 1e-08, S = 998, I = 2, mean_lbeta = -6.61306611533102, 
    sd_lbeta = 0.398458016908828), .Names = c("time", "S", "I", "mean_lbeta", "sd_lbeta")))
    ## single recovery
    expect_message(r5 <- run_sim(nt=1e-8,rptfreq=1e-8,R0_init=0.01,Ivec=1,
                                 discrete=FALSE,seed=101,useCpp=TRUE),"system went extinct")
    expect_equal(c(r5),
                 structure(list(time = 1e-08, S = 1000, I = 0, mean_lbeta = NaN, 
                  sd_lbeta = NaN), .Names = c("time", "S", "I", "mean_lbeta", "sd_lbeta")))
    expect_message(r6 <- run_sim(nt=1e-8,rptfreq=1e-8,R0_init=0.01,
                                Ivec=1,discrete=FALSE,seed=101,useCpp=FALSE),
                  "system went extinct")
    expect_equal(c(r5),c(r6))
})


r <- c(0.6,0.3,0.05,0.05)
my_sample(r,sum(r))
tt <- prop.table(table(replicate(10000,my_sample(r,sum(r)))))
