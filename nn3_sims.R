## nearly neutral sims ... Latin hypercube

source("nn2_funs.R")
nsim <- 1000
lhs_df <- cbind(
    mu=10^seq(-3,-1,length.out=nsim),
    mut_mean=seq(-4,-0.5,length.out=nsim),
    mut_sd=10^seq(-1,0.5,length.out=nsim),
    N=round(10^seq(0.5,3,length.out=nsim)),
    gamma=1/round(seq(2,10,length.out=nsim)))
set.seed(101)
for (i in 2:ncol(lhs_df)) {
    lhs_df[,i] <- sample(lhs_df[,i])
}
lhs_df <- data.frame(run=1:nrow(lhs_df),lhs_df)
saveRDS(lhs_df,file="lhs_df.rds")
## pairs(lhs_df,pch=".",gap=0)
fn <- "nn_runs_3.rda"
 resList <- list()
 for (i in 1:nrow(lhs_df)) {
     pars <- lhs_df[i,]
     cat(c(i,pars,"\n")
         argList <- c(list(nt=1e6,rptfreq=1000,R0_init=4,seed=101,
                           progress=TRUE),
                      as.list(pars))
         resList[[i]] <- try(do.call(run_sim,argList))
         save("lhs_df","resList",file=fn)
     }
 }
