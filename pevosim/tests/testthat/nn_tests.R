## miscellanous messing-around and tests
library(testthat)
library(pevosim)

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

r1 <- run_sim(nt=1e2,rptfreq=1,discrete=FALSE,seed=101)
sink("cppdebug.txt")
r2 <- run_sim(nt=1e2,rptfreq=1,discrete=FALSE,seed=101,useCpp=TRUE,debug=TRUE)
sink()

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
    system.time(res1A <- run_sim(nt=1e4,rptfreq=10,seed=101))
    system.time(res1B <- run_sim(nt=1e4,rptfreq=10,seed=101,discrete=FALSE))
    save("res1A","res1B",file="../simdata/simdisc.rda")
    load("../simdata/simdisc.rda")
    library(ggplot2); theme_set(theme_bw())
    library(tidyr)
    library(dplyr)
    rL <- lapply(list(res1A,res1B),gather,key=var,value=val,-time) %>%
        setNames(c("disc","cont")) %>% bind_rows(.id="type")
    ggplot(rL,aes(time,val,colour=type))+geom_line()+
        facet_wrap(~var,scale="free")+
        scale_colour_brewer(palette="Set1")
}
