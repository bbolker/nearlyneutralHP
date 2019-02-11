#######
## Run run_sim many times across either parameter values or, first,
## a single set of parameter values to capture magnitude of the effect of stochasticity 
## as a debug like tool
#######

library(arm); library(ggplot2); library(gridExtra); library(dplyr) 
source("ggplot_theme.R")

nt         <- 1e6
num_points <- 1000
rptfreq    <- max(nt / num_points, 1) 
nrpt       <- nt %/% rptfreq
nruns      <- 7
randseed   <- sample(1:1e5, nruns, replace = FALSE)

## Set up parameters in a manner to correspond to output data frame
params   <- data.frame(
   nt                  = rep(nt, nruns)
 , rptfreq             = rep(rptfreq, nruns)
 , nrpt                = rep(nrpt, nruns)
 , mut_var             = rep("beta", nruns)
 , run                 = seq(1, nruns, 1)
 , seed                = randseed
 , d                   = rep(0.01, nruns)
 , mu                  = rep(0.001, nruns)
 , mut_mean            = rep(-0.05, nruns)    
 , mut_sd              = rep(0.05, nruns)   
 , alpha0              = rep(0.05, nruns)
 , tol0                = rep(1, nruns)
 , res0                = rep(1, nruns)
 , mut_host_mean_shift = rep(2, nruns)        
 , mut_host_sd_shift   = rep(1, nruns)
 , mut_host_mu_shift   = rep(4, nruns) ## To make sure the host never evolves set this to a *really* big number
 , mut_host_res_bias   = rep(0.5, nruns)
 , host_dyn_only       = rep(FALSE, nruns)
 , mut_mean            = rep(-1, nruns)
 , mut_sd              = rep(0.5, nruns)
 , power_c             = rep(0.75, nruns)
 , power_exp           = rep(2, nruns)
 , b_decay             = rep(2.3, nruns)
 , b                   = rep(0.4, nruns)
 , N                   = rep(1000, nruns)
 , balance_birth       = rep(TRUE, nruns)
 , stochastic_birth    = rep(FALSE, nruns)
  )

for (i in 1:nruns) {
  
  print(i / nruns)
  
  time_check <- print(system.time(res_1000 <- with(params
        , run_sim(
   nt                  = nt[i]
 , rptfreq             = rptfreq[i]
 , mut_var             = mut_var[i]
 , seed                = seed[i]
 , d                   = d[i]
 , mu                  = mu[i]
 , alpha0              = 0.25
 , tol0                = tol0[i]
 , res0                = res0[i]
 , mut_host_mean_shift = mut_host_mean_shift[i]
 , mut_host_sd_shift   = mut_host_sd_shift[i]
 , mut_host_mu_shift   = mut_host_mu_shift[i]
 , mut_host_res_bias   = mut_host_res_bias[i]
 , host_dyn_only       = host_dyn_only[i]
 , mut_mean            = mut_mean[i]
 , mut_sd              = mut_sd[i]
 , power_c             = power_c[i]
 , power_exp           = power_exp[i]
 , b_decay             = b_decay[i]
 , b                   = b[i]
 , N                   = N[i]
 , balance_birth       = balance_birth[i]
 , stochastic_birth    = stochastic_birth[i]
 , debug2              = FALSE
 , debug3              = TRUE
 , progress            = TRUE
 , R0_init             = 400 #652.3262
      ))))
 
  ## clean up run i (add parameters and remove NA if the system went extinct)   
res_1000 <- res_1000 %>% mutate(run = i, elapsed_time = time_check[3])  
res_1000 <- left_join(res_1000, params, by = "run")
res_1000 <- res_1000[complete.cases(res_1000), ]

  if (i == 1) {
res_1000_all <- res_1000
  } else {
res_1000_all <- rbind(res_1000_all, res_1000)
  }

}

#saveRDS(res_1000_all, "res_1000.Rds")

#### ggplots of run

  gg1.1 <- ggplot(res_1000_all, aes(time / nrpt, plogis(mean_plbeta)))  + geom_line() + xlab("Time (steps)") + ylab("Parasite Intrinsic Beta")
  gg1.2 <- ggplot(res_1000_all, aes(time / nrpt, power_tradeoff(plogis(mean_plalpha), params$power_c[1], params$power_exp[1]) * 
      plogis(mean_plbeta))) + geom_line() + xlab("Time (steps)") + ylab("Parasite Realized Beta")
  gg1.3 <- ggplot(res_1000_all, aes(time / nrpt, plogis(mean_plalpha))) + geom_line() + xlab("Time (steps)") + ylab("Parasite alpha")
  gg1.4 <- ggplot(res_1000_all, aes(time / nrpt
    , (( power_tradeoff(plogis(mean_plalpha), params$power_c[1], params$power_exp[1]) * plogis(mean_plbeta) ) / 
    ( plogis(mean_plalpha) + 0.2 + params$d[1] )))) + geom_line() + xlab("Time (steps)") + ylab("Parasite R0")
  gg1.5 <- ggplot(res_1000_all, aes(time / nrpt, exp(mean_hlres)))      + geom_line() + xlab("Time (steps)") + ylab("Host Resistance")
  gg1.6 <- ggplot(res_1000_all, aes(time / nrpt, exp(mean_hltol)))      + geom_line() + xlab("Time (steps)") + ylab("Host Tolerance")

  gg2.1 <- ggplot(res_1000_all, aes(time / nrpt, num_I_strains))        + geom_line() + xlab("Time (steps)") + ylab("Number of I Strains")
  gg2.2 <- ggplot(res_1000_all, aes(time / nrpt, num_S_strains))        + geom_line() + xlab("Time (steps)") + ylab("Number of S Strains")
 
  gg3.1 <- ggplot(res_1000_all, aes(time / nrpt, num_S))                + geom_line() + xlab("Time (steps)") + ylab("Number of S Individuals")
  gg3.2 <- ggplot(res_1000_all, aes(time / nrpt, num_I))                + geom_line() + xlab("Time (steps)") + ylab("Number of I Individuals")
  gg3.3 <- ggplot(res_1000_all, aes(time / nrpt, pop_size))             + geom_line() + xlab("Time (steps)") + ylab("Total population Size")


#### ggplots tracking progress on tradeoff space
power_trade_dat <- data.frame(
  alpha = seq(0.01, 1.0, by = 0.01)
, beta  = power_tradeoff(alpha = seq(0.01, 1.0, by = 0.01), c = params$power_c[1], curv = params$power_exp[1])
, R0    = power_tradeoff(alpha = seq(0.01, 1.0, by = 0.01), c = params$power_c[1], curv = params$power_exp[1]) / 
    ( seq(0.01, 1.0, by = 0.01) + 0.2 + params$d[1] )
)

gg4.1 <- ggplot(power_trade_dat, aes(alpha, beta)) + geom_line() + 
  geom_point(data = res_1000_all, aes(
    plogis(mean_plalpha)
  , power_tradeoff(plogis(mean_plalpha), params$power_c[1], params$power_exp[1]) * plogis(mean_plbeta)
    )) + 
  geom_line(data = res_1000_all, aes(
    plogis(mean_plalpha)
  , power_tradeoff(plogis(mean_plalpha), params$power_c[1], params$power_exp[1]) * plogis(mean_plbeta)
  ))

gg4.2 <- ggplot(power_trade_dat, aes(alpha, R0)) + geom_line() + 
  geom_point(data = res_1000_all, aes(
    plogis(mean_plalpha)
  , ( power_tradeoff(plogis(mean_plalpha), params$power_c[1], params$power_exp[1]) * plogis(mean_plbeta) ) / 
    ( plogis(mean_plalpha) + 0.2 + params$d[1] )  
    )) + 
  geom_line(data = res_1000_all, aes(
    plogis(mean_plalpha)
  , ( power_tradeoff(plogis(mean_plalpha), params$power_c[1], params$power_exp[1]) * plogis(mean_plbeta) ) /
    ( plogis(mean_plalpha) + 0.2 + params$d[1] )
  ))

grid.arrange(gg1.2, gg4.1, gg3.1, gg1.1, gg4.2, gg3.2, gg1.3, gg1.4, gg2.1, ncol = 3, nrow = 3)


gg1.2b <- ggplot(res_1000_all, aes(time / nrpt, mean_plbeta)) + geom_line(aes(colour = as.factor(mut_mean))) +
  xlab("Time (steps)") + ylab("Parasite Intrinsic Beta")
gg1.3b <- ggplot(res_1000_all, aes(time / nrpt, mean_plalpha)) + geom_line(aes(colour = as.factor(mut_mean))) +
  xlab("Time (steps)") + ylab("Parasite Alpha") + geom_hline(aes(yintercept = 0.855))
grid.arrange(gg1.2b, gg1.3b, ncol = 1)

## Numbers at the end of the name correspond to slide numbers in Res_to_explore keynote
res_1000_all_48 <- res_1000_all
