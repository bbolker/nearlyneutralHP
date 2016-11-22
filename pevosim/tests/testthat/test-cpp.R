library(pevosim)
library(testthat)
expect_equal(get_rates(betavec=1:3,gamma=2:4,Ivec=0:2,S=10),
             c(0,20,60,0,3,8))

L <- list(x=4,y=3,ltraitvec=numeric(0),Ivec=integer(0))
set.seed(101)
do_mutC(L,mut_var="y",
       orig_trait=3,mut_mean=0,mut_sd=1)
expect_equal(L,
             structure(list(x = 4, y = c(3, 14.4973157193337),
                            ltraitvec = 2.67396350948461, 
                            Ivec = 1L),
                       .Names = c("x", "y", "ltraitvec", "Ivec")))

do_extinct(L,"y",1)
expect_equal(L,
             structure(list(x = 4, y = 14.4973157193337, ltraitvec = numeric(0), 
    Ivec = integer(0)), .Names = c("x", "y", "ltraitvec", "Ivec"
)))

sane_result <- function(state,N=101) {
    stopifnot(sum(state$Ivec)+state$S==N)
}
    
S0 <- list(ltraitvec=log(2/101),beta=2/101,gamma=1,Ivec=1L,S=100)
P0 <- list(mu=0,mut_sd=0,mut_mean=0,mut_var="beta",mut_link="log")
set.seed(101)
run_stepC(S0,t_tot=0,t_end=1,params=P0,debug=TRUE)
sane_result(S0,N=101)

S0 <- list(ltraitvec=log(2/1001),beta=2/1001,gamma=1,Ivec=1L,S=1000)

## non-zero mutation but mutation sd=0
##S0 <- list(ltraitvec=1,beta=2,gamma=1,Ivec=1L,S=100)
P1 <- list(mu=0.5,mut_sd=0,mut_mean=0,mut_var="beta",mut_link="log")
run_stepC(S0,t_tot=0,t_end=1,params=P0,debug=FALSE)
S0$S ## 13
