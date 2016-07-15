## nearly neutral sims ...

## going to model logit-beta, for now (lbeta)
## can easily switch to another scale

## params
## N: population size
## mu: per-infection mutation probability
## gamma: recovery rate
## lbeta0: beta of initial

run_sim <- function(R0_init=2,
                    gamma=1/5,
                    N=1000,
                    mu=0.01,
                    mut_type="shift",
                    mut_mean=-1,
                    mut_sd=0.5,
                    nt=100000,
                    rptfreq=100) {

    ## R0 = beta*N/gamma
    lbeta0 <- qlogis(R0_init*gamma/N)

    ## start with a vector of size 1 ...
    lbetavec <- c(lbeta0)
    Ivec <- N/2
    S <- N-sum(Ivec)

    nrpt <- nt %/% rptfreq
    nq <- 9
    qvec <- seq(0,1,length=nq+2)[-c(1,nq+2)]
    res <- as.data.frame(matrix(NA,nrow=nrpt,ncol=5,
            dimnames=list(NULL,c("time","S","I","mean_lbeta","sd_lbeta"))))
    for (i in 1:nt) {
        ## cat("time ",i,"\n")
        betavec <- plogis(lbetavec)
        ## cat("betavec:",betavec,"\n")
        ## cat("Ivec:",Ivec,"\n")
        uninf <- rbinom(1,size=S,prob=prod((1-betavec)^Ivec))
        ## 'prob' is internally normalized
        newinf <- drop(rmultinom(1,size=S-uninf,prob=Ivec*betavec))
        stopifnot(uninf+sum(newinf)==S)
        recover <- rbinom(length(Ivec),size=Ivec,prob=gamma)
        mutated <- rbinom(length(newinf),size=newinf,prob=mu)
        stopifnot(length(recover)==length(mutated))
        stopifnot(length(newinf)==length(Ivec))
        Ivec <- Ivec-recover+(newinf-mutated)
        ## cat("I,uninf, unmutated, mutated, recovered",
        ## c(sum(Ivec),uninf,sum(newinf-mutated),sum(mutated),sum(recover)),
        ## "\n")
        ## now do mutation ...
        if (mut_type=="shift") {
            tot_mut <- sum(mutated)
            newlbeta <- rep(lbetavec,mutated)+
                rnorm(tot_mut,mut_mean,mut_sd)
            Ivec <- c(Ivec,rep(1,tot_mut))
            lbetavec <- c(lbetavec,newlbeta)
        }
        if (all(Ivec==0)) break
        if (length(extinct <- which(Ivec==0))>0) {
            lbetavec <- lbetavec[-extinct]
            Ivec <- Ivec[-extinct]
        }
        S <- S + sum(recover) - sum(newinf)
        stopifnot(length(S)==1)
        stopifnot(sum(Ivec)+S == N)
        lbeta_mean <- sum(Ivec*lbetavec)/sum(Ivec)
        lbeta_sd <- sqrt(sum(Ivec*(lbetavec-lbeta_mean)^2)/sum(Ivec))
        if (i %% rptfreq == 0)
            res[i %/% rptfreq,] <- c(i,S,sum(Ivec),lbeta_mean,lbeta_sd)
    }
    return(res)
}

set.seed(101)
system.time(res <- run_sim(nt=1e6,rptfreq=1000))
## 3 minutes for 10^6 steps: profile?
library(ggplot2); theme_set(theme_bw())
library(tidyr)
library(dplyr)
vtype <- data.frame(var=c("S","I","mean_lbeta","sd_lbeta"),
                    type=rep(c("pop","lbeta"),c(2,2)))
resm <- gather(res,var,val,-time) %>% full_join(vtype)
g0 <- ggplot(resm,aes(time,val,colour=var))+
    geom_line()+facet_wrap(~var,scale="free")

g0 %+% subset(resm,time>1e5)+
    geom_smooth()

## hard to tell whether stable ...


system.time(res2 <- run_sim(nt=1e6,rptfreq=1000,N=100))
res2m <- gather(res2,var,val,-time) %>% full_join(vtype)
g0 %+% res2m
