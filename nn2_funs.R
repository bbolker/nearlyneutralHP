## going to model logit-beta, for now (lbeta)
## can easily switch to another scale

## params
## N: population size
## mu: per-infection mutation probability
## gamma: recovery rate
## lbeta0: beta of initial

run_sim <- function(R0_init=2,  ## >1
                    gamma =1/5,  ## >0
                    N=1000,     ## integer >0
                    mu=0.01,    ## >0
                    mut_type="shift",
                    mut_mean=-1,## <0 (for sensibility)
                    mut_sd=0.5, ## >0
                    mut_var="beta",
                    mut_link=NULL,
                    nt=100000,
                    rptfreq=100, ## divides nt?
                    seed=NULL,
                    progress=FALSE) {

    ## TO DO:
    ##  - profile
    ##  - add progress bar?
    ##  - allow immigration?
    ##  - SIRS model?
    ##  - mutation of gamma?

    if (!is.null(seed)) set.seed(seed)

    Ivec <- round(N*(1-1/R0_init))
    S <- N-sum(Ivec)

    ## R0 = beta*N/gamma
    beta0 <- R0_init*gamma/N

    if (mut_var=="beta") {
        mut_link <- make.link("logit")
        ltraitvec <- mut_link$linkfun(beta0)
    } else {
        mut_link <- make.link("log")
        ltraitvec <- mut_link$linkfun(gamma)
    }
    
    nrpt <- nt %/% rptfreq
    nq <- 9
    qvec <- seq(0,1,length=nq+2)[-c(1,nq+2)]
    res <- as.data.frame(matrix(NA,nrow=nrpt,ncol=5,
            dimnames=list(NULL,c("time","S","I",
                                 paste0(c("mean_l","sd_l"),mut_var)))))

    gammavec <- gamma
    betavec <- beta0
    
    for (i in 1:nt) {
        ## cat("time ",i,"\n")
        traitvec <- mut_link$linkinv(ltraitvec)
        if (mut_var=="beta") {
            betavec <- traitvec
        } else if (mut_var=="gamma") {
            gamma <- traitvec
        }
        
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
            if (tot_mut>0) {
                newltrait <- rep(ltraitvec,mutated)+
                    rnorm(tot_mut,mut_mean,mut_sd)
                Ivec <- c(Ivec,rep(1,tot_mut))
                ltraitvec <- c(ltraitvec,newltrait)
            }
        }
        if (all(Ivec==0)) {
            message(sprintf("system went extint prematurely (t=%d)",i))
            break
        }
        if (length(extinct <- which(Ivec==0))>0) {
            ltraitvec <- ltraitvec[-extinct]
            Ivec <- Ivec[-extinct]
        }
        S <- S + sum(recover) - sum(newinf)
        stopifnot(length(S)==1)
        stopifnot(sum(Ivec)+S == N)
        ltrait_mean <- sum(Ivec*ltraitvec)/sum(Ivec)
        ltrait_sd <- sqrt(sum(Ivec*(ltraitvec-ltrait_mean)^2)/sum(Ivec))
        if (i %% rptfreq == 0) {
            if (progress) cat(".")
            res[i %/% rptfreq,] <- c(i,S,sum(Ivec),ltrait_mean,ltrait_sd)
        }
    }
    return(res)
}

