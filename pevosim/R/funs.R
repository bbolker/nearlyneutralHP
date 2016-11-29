## going to model logit-beta, for now (lbeta)
## can easily switch to another scale

#' @useDynLib pevosim
#' @importFrom Rcpp evalCpp
## params
## N: population size
## mu: per-infection mutation probability
## gamma: recovery rate
## lbeta0: beta of initial
get_mut <- function(orig_trait,
                    mut_mean,mut_sd,mut_type="shift") {
    if (mut_type=="shift") {
        new_trait <- orig_trait+
            rnorm(length(orig_trait),mut_mean,mut_sd)
        if (any(is.na(new_trait))) stop("??")
        return(new_trait)
    } else stop("unknown mut_type")
}

do_mut <- function(state,mut_var,mut_link,orig_trait,...) {
    new_trait <- get_mut(orig_trait,...)
    state[[mut_var]] <- c(state[[mut_var]],mut_link$linkinv(new_trait))
    state$ltraitvec <- c(state$ltraitvec,new_trait)
    state$Ivec <- c(state$Ivec,rep(1,length(orig_trait)))
    return(state)
}

do_extinct <- function(state,mut_var,extinct) {
    state[[mut_var]] <- state[[mut_var]][-extinct]
    state$ltraitvec <- state$ltraitvec[-extinct]
    state$Ivec <- state$Ivec[-extinct]
    return(state)
}

##' Logit transform on (min,max) rather than (0,1)
##' @param minval minimum value
##' @param maxval maximum value
##' @export
multlogit <- function(minval=0,maxval=1,scale=1) {
    delta <- maxval-minval
    ## R CMD check will always complain about this.
    ## C_logit_link is not otherwise externally accessible.
    ##  could write my own ...
    linkfun <- function(mu) .Call(stats:::C_logit_link, (mu-minval)/delta)
    linkinv <- function(eta) .Call(stats:::C_logit_linkinv, eta)*delta + minval
    mu.eta <- function(eta) .Call(stats:::C_logit_mu_eta, eta)*delta
    valideta <- function(eta) TRUE
    ## need to be able to find C_logit_mu_eta ...
    ## environment(linkfun) <- environment(linkinv) <- environment(mu.eta) <- environment(valideta) <- asNamespace("stats")
    structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
                   valideta = valideta, name = "multlogit"), class = "link-glm")
}

#' Evolutionary simulations
#' 
#' @param R0_init starting value of R0
#' @param gamma0 recovery rate
#' @param gamma_max  ?? not used ??
#' @param N population size
#' @param mut_type mutation model (only "shift" implemented)
#' @param mut_mean mean value of mutations
#' @param mut_sd mutational standard deviations
#' @param mut_var which variable(s) mutate?
#' @param mut_link link function for mutation
#' @param Ivec initial vector of infected numbers
#' @param mu mutation probability per replication
#' @param discrete (logical) run discrete-time sim?
#' @param dt time step
#' @param seed random-number seed
#' @param nt number of time steps
#' @param rptfreq reporting frequency (should divide nt)
#' @param progress draw progress bar?
#' @param debug (logical) debugging output?
#' @param useCpp (logical) use C++ code for continuous-time models?
#' @importFrom stats make.link rnorm rbinom rmultinom rexp runif
#' @export
run_sim <- function(R0_init=2,  ## >1
                    gamma0 =1/5,  ## >0
                    gamma_max = Inf,
                    N=1000,     ## integer >0
                    mu=0.01,    ## >0
                    mut_type="shift",
                    mut_mean=-1,## <0 (for sensibility)
                    mut_sd=0.5, ## >0
                    mut_var="beta",
                    mut_link=NULL,
                    Ivec=NULL,
                    dt=1,
                    nt=100000,
                    rptfreq=max(nt/500,1), ## divides nt?
                    discrete=TRUE, ## discrete-time?
                    seed=NULL,
                    progress=FALSE,
                    debug=FALSE,
                    useCpp=FALSE) {

    if (round(N)!=N) {
        warning("rounding N")
        N <- round(N)
    }
    if (mut_mean>0) {
        warning("positive mutation bias")
    }
    stopifnot(R0_init>0,
              gamma0>0,
              N>0,
              mu>0,
              mut_sd>0,
              (nt/rptfreq) %% 1 ==0
              )
    dfun <- function(lab="") {
        if (debug) {
            cat(lab,"\n")
            print(state)
        }
        return(NULL)
    }

    ## TO DO:
    ##  - profile
    ##  - add progress bar?
    ##  - allow immigration?
    ##  - SIRS model?
    ##  - mutation of gamma?

    if (!is.null(seed)) set.seed(seed)

    if (is.null(Ivec)) {
        ## start at equilibrium I ...
        Ivec <- max(1,round(N*(1-1/R0_init)))
    }
    S <- N-sum(Ivec)

    ## R0 = beta*N/gamma
    beta0 <- R0_init*gamma0/N

    ## set up initial trait vector
    if (mut_var=="beta") {
        if (is.null(mut_link)) {
            mut_link <- if (discrete) make.link("logit") else make.link("log")
        }
        ltraitvec <- mut_link$linkfun(beta0)
    } else {
        if (is.null(mut_link)) mut_link <- make.link("log")
        ltraitvec <- mut_link$linkfun(gamma0)
    }
    ## parameters structure (parallel vectors), so these can
    ## be modified via function and passed back ...
    state <- list(beta=beta0,gamma=gamma0,
                  ltraitvec=ltraitvec,Ivec=Ivec,
                  S=S)

    dfun("init")
    
    nrpt <- nt %/% rptfreq
    ## not used: intended for accumulating information about a binned
    ##  distribution over time
    ## nq <- 9
    ## qvec <- seq(0,1,length=nq+2)[-c(1,nq+2)]
    res <- as.data.frame(matrix(NA,nrow=nrpt,ncol=5,
                                dimnames=list(NULL,c("time","S","I",
                                                     paste0(c("mean_l","sd_l"),mut_var)))))

    t_tot <- 0  ## for continuous model


    for (i in 1:nrpt) {
        ## cat("time ",i,"\n")
        if (discrete) {
            for (j in 1:rptfreq) {
                for (t0 in seq(round(1/dt))) {
                ## cat("betavec:",betavec,"\n")
                ## cat("Ivec:",Ivec,"\n")
                ## prob of escaping infection completely
                uninf <- rbinom(1,size=state$S,
                                prob=prod((1-state$beta)^state$Ivec))
                ## division of new infectives among strains
                ## 'prob' is internally normalized
                newinf <- drop(rmultinom(1,size=state$S-uninf,
                                         prob=state$Ivec*state$beta))
                stopifnot(uninf+sum(newinf)==state$S)
                recover <- rbinom(length(state$Ivec),
                                  size=state$Ivec,prob=state$gamma)
                ## fraction of new infections -> mutation
                mutated <- rbinom(length(newinf),size=newinf,prob=mu)
                stopifnot(length(recover)==length(mutated))
                stopifnot(length(newinf)==length(state$Ivec))
                ## update infectives
                state$Ivec <- state$Ivec-recover+(newinf-mutated)
                ## cat("I,uninf, unmutated, mutated, recovered",
                ## c(sum(Ivec),uninf,sum(newinf-mutated),sum(mutated),sum(recover)),
                ## "\n")
                dfun("before mutation")
                if (debug) print(mutated)
                
                ## now do mutation ...
                tot_mut <- sum(mutated)
                if (tot_mut>0) {
                    state <- do_mut(state,mut_var,
                                    mut_link,
                                    orig_trait=rep(state$ltraitvec,mutated),
                                    mut_mean=mut_mean,
                                    mut_sd=mut_sd,
                                    mut_type=mut_type)
                }
                if (all(state$Ivec==0)) {
                    message(sprintf("system went extinct prematurely (t=%d)",i))
                    break
                }
                if (length(extinct <- which(state$Ivec==0))>0) {
                    state <- do_extinct(state,mut_var,extinct)
                }
                dfun("after mutation")
                if (length(state$Ivec)>0) {
                    ## avoid X+numeric(0)==numeric(0) problem ...
                    state$S <- state$S + sum(recover) - sum(newinf)
                }
                stopifnot(length(state$S)==1)
                stopifnot(sum(state$Ivec)+state$S == N)
            }  ## rptfreq time steps
        } else  ## end discrete case
        {
            if (useCpp) {
                ## cat("useCpp: state: \n")
                ## print(state)
                run_stepC(state,
                          (i-1)*rptfreq,
                          i*rptfreq,
                params=list(mu=mu,
                            mut_sd=mut_sd,
                            mut_mean=mut_mean,
                            mut_var=mut_var,
                            mut_link=""), ## FIXME
                debug=debug)
            } else { ## use R
                while (t_tot<(i+1)*rptfreq) {
                    nstrain <- length(state$ltraitvec)
                    rates <- get_rates(state)
                    if (debug) cat(t_tot,rates,"\n")
                    t_tot <- t_tot + rexp(1,sum(rates))
                    w <- sample(length(rates),size=1,prob=rates)
                    event <-  ((w-1) %/% nstrain) + 1
                    strain <- ((w-1) %% nstrain) + 1
                    if (event==1) { ## infection
                        if (runif(1)<mu) {
                            state <- do_mut(state,mut_var,
                                            mut_link,
                                            state$ltraitvec[strain],
                                            mut_mean,mut_sd)
                        } else {
                            state$Ivec[strain] <- state$Ivec[strain]+1
                        }
                        state$S <- state$S - 1
                    } ## infection

                    if (event==2) { ## recovery
                        state$Ivec[strain] <- state$Ivec[strain]-1
                        state$S <- state$S+1
                        if (state$Ivec[strain]==0) {
                            state <- do_extinct(state,mut_var,strain)
                        }
                    }
                    stopifnot(length(state$S)==1)
                    stopifnot(sum(state$Ivec)+state$S == N)
                } ## loop until end of time period
            } ## use R, not Rcpp
        } ## continuous time
        ## summary statistics
        I_tot <- sum(state$Ivec)
        ltrait_mean <- sum(state$Ivec*state$ltraitvec)/I_tot
        ltrait_sd <- sqrt(sum(state$Ivec*(state$ltraitvec-ltrait_mean)^2)/I_tot)
        if (progress) cat(".")
        res[i,] <- c(i*rptfreq,state$S,I_tot,ltrait_mean,ltrait_sd)
        ## DRY ...
        if (all(state$Ivec==0)) {
            message(sprintf("system went extinct prematurely (t=%d)",i))
            break
        }
    } ## loop over reporting frequencies
    if (progress) cat("\n")                      
    attr(res,"mut_link") <- mut_link
    return(res)
}


#' Get rates
#' @param state list of state vectors
#' @export
get_rates <- function(state,dt=1) {
    inf_rates <- state$beta*state$Ivec*state$S*dt
    recover_rates <- state$gamma*state$Ivec*dt
    return(c(inf_rates,recover_rates))
}
