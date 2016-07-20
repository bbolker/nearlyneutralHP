source("nn2_funs.R")

fn <- "nn_runs_4.rda"
if (!file.exists(fn)) {
    set.seed(101)
    Nvals <- c(1000,100,33)
    resG <- lapply(Nvals,
                   run_sim,
                   nt=1e6,rptfreq=1000,
                   mut_var="gamma",
                   seed=101,
                   progress=TRUE,
                   R0_init=4)
    save("resG",file=fn)
}

