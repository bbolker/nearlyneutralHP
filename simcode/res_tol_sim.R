#######
## Run run_sim many times across either parameter values or, first,
## a single set of parameter values to capture magnitude of the effect of stochasticity 
## as a debug like tool
#######

library(arm); library(ggplot2); library(gridExtra); library(dplyr) 
source("ggplot_theme.R")

nt         <- 1e5
num_points <- 500
rptfreq    <- max(nt / num_points, 1) 
nrpt       <- nt %/% rptfreq
nruns      <- 50
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
 , mu                  = rep(0.01, nruns)
 , alpha0              = rep(0.01, nruns)
 , tol0                = rep(1, nruns)
 , res0                = rep(1, nruns)
 , mut_host_mean_shift = rep(2, nruns)
 , mut_host_sd_shift   = rep(1, nruns)
 , mut_host_mu_shift   = rep(5, nruns)
 , mut_host_res_bias   = rep(0.5, nruns)
 , host_dyn_only       = rep(FALSE, nruns)
 , mut_mean            = rep(-1, nruns)
 , mut_sd              = rep(0.5, nruns)
 , power_c             = rep(0.75, nruns)
 , power_exp           = rep(2, nruns)
 , b_decay             = rep(2.3, nruns)
 , b                   = rep(0.2, nruns)
 , N                   = rep(1000, nruns)
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
 , alpha0              = alpha0[i]
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
 , debug2              = FALSE
 , debug3              = FALSE
 , progress            = TRUE
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

saveRDS(res_1000_all, "res_1000.Rds")

  (gg1.1 <- ggplot(res_1000_all, aes(time / nrpt, plogis(mean_plbeta)))  + geom_line() + xlab("Time (steps)") + ylab("Parasite Intrinsic Beta"))
  (gg1.2 <- ggplot(res_1000_all, aes(time / nrpt, plogis(mean_plalpha))) + geom_line() + xlab("Time (steps)") + ylab("Parasite Intrinsic Alpha"))
  (gg1.3 <- ggplot(res_1000_all, aes(time / nrpt, exp(mean_hlres)))      + geom_line() + xlab("Time (steps)") + ylab("Host Resistance"))
  (gg1.4 <- ggplot(res_1000_all, aes(time / nrpt, exp(mean_hltol)))      + geom_line() + xlab("Time (steps)") + ylab("Host Tolerance"))

  (gg2.1 <- ggplot(res_1000_all, aes(time / nrpt, num_I_strains))        + geom_line() + xlab("Time (steps)") + ylab("Number of I Strains"))
  (gg2.2 <- ggplot(res_1000_all, aes(time / nrpt, num_S_strains))        + geom_line() + xlab("Time (steps)") + ylab("Number of S Strains"))
 
  (gg3.1 <- ggplot(res_1000_all, aes(time / nrpt, num_S))                + geom_line() + xlab("Time (steps)") + ylab("Number of S Individuals"))
  (gg3.2 <- ggplot(res_1000_all, aes(time / nrpt, num_I))                + geom_line() + xlab("Time (steps)") + ylab("Number of I Individuals"))
  (gg3.3 <- ggplot(res_1000_all, aes(time / nrpt, pop_size))             + geom_line() + xlab("Time (steps)") + ylab("Total population Size"))


  
  ggplot(res_1000c, aes(time / 200, num_S)) + geom_line() + geom_line(aes(time / 200, num_I), colour = "red") +
    xlab("Time (x200 steps)") + ylab("Number of S and I Individuals")

