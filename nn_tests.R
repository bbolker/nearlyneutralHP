## miscellanous messing-around and tests
library(testthat)

source("nn2_funs.R")

## basic run
res1 <- run_sim(nt=1e4,rptfreq=100,seed=101)
expect_equal(res1[100,],
     structure(list(10000, 160, 840, 1.50094828644975, 
                     0.11040310840802),
               .Names = c("time", "S", "I", "mean_lbeta", "sd_lbeta"),
               row.names = 100L, class = "data.frame"))

## mutating gamma
res2 <- run_sim(nt=1e4,rptfreq=100,seed=101,mut_var="gamma")
expect_equal(res2[100,],
             structure(list(10000, 0, 1000, -8.14914364919813, 
                            0.547470943260739),
                       .Names = c("time", "S", "I", "mean_lgamma", "sd_lgamma"),
                       row.names = 100L, class = "data.frame"))

Lfun0 <- make.link("logit")
Lfun <- multlogit(0,1)
all.equal(Lfun$linkinv(as.numeric(-5:5)),Lfun0$linkinv(as.numeric(-5:5)))
res3 <- run_sim(nt=1e4,rptfreq=100,seed=101,mut_var="gamma",
                mut_link=multlogit(0.1,0.5))
## runs, but ... ?? seems to push up against minimum value
  structure(list(time = 10000, S = 332, I = 668, mean_lgamma = -14.202599029542, 
    sd_lgamma = 1.03017299599263), .Names = c("time", "S", "I", 
"mean_lgamma", "sd_lgamma"), row.names = 100L, class = "data.frame")
