library(pevosim)
library(qrng)
fn <- "nn_runs_8.rda"
nsim <- 100
ranges <- list(mu=c(-3,-1),
               mut_mean=c(-4,-0.5),
               mut_sd=c(-1,0.5),
               N=c(1,3),
               gamma0=c(2,10))
npar <- length(ranges)
pow10 <- function(x) 10^x
rpow10 <- function(x) round(10^x)
inv10 <- function(x) 10^(-x)
trans <- list(pow10,identity,pow10,rpow10,inv10)

ss <- sobol(nsim,npar)
for (i in 1:npar) {
    ss[,i] <- ss[,i]*diff(ranges[[i]])+ranges[[i]][1]
    ss[,i] <- trans[[i]](ss[,i])
}
colnames(ss) <- names(ranges)
ss_df <- data.frame(run=1:nrow(ss),ss)
resList <- list()
options(digits=3)
for (i in 1:nrow(ss_df)) {
    pars <- ss_df[i,-1]
    cat(i,paste(names(pars),format(unlist(pars)),sep="="),"\n")
    argList <- c(list(nt=1e7,R0_init=20,seed=101,
                      mut_var="beta",
                      discrete=FALSE,
                      useCpp=TRUE,
                      progress=TRUE,
                      debug=FALSE),
                 as.list(pars))
    resList[[i]] <- try(do.call(run_sim,argList))
    cat("done\n")
    save("ss_df","resList",file=fn)
}



