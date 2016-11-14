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

test_that("multlogit", {
    Lfun0 <- make.link("logit")
    Lfun <- multlogit(0,1)
    all.equal(Lfun$linkinv(as.numeric(-5:5)),Lfun0$linkinv(as.numeric(-5:5)))
})
