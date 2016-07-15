## nearly neutral sims ...

source("nn2_funs.R")

fn <- "nn_runs_1.rda"
if (!file.exists(fn)) {
    set.seed(101)
    ## 3-4 minutes for 10^6 steps: profile?
    print(system.time(res_1000 <- run_sim(nt=1e6,rptfreq=1000,
                                     seed=101,progress=TRUE)))
    ## 140 seconds ...
    print(system.time(res_100 <- run_sim(nt=1e6,rptfreq=1000,N=100,
                                      seed=101,progress=TRUE)))
    ## switched to R0_init=4; goes extinct otherwise
    print(system.time(res_33 <- run_sim(nt=1e6,rptfreq=1000,N=33,
                                      R0_init=4,seed=101,
                                      progress=TRUE)))
    save("res_1000","res_100","res_33",file=fn)
}

