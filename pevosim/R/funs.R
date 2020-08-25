## going to model logit-beta, for now (lbeta)
## can easily switch to another scale

#######
## MK: search ** for new additions:
## 1) Gamma change sign of mutation mean and link scale
## 2) Host evolution in an extremely naive way (no tracking of host type, so always assume selective sweep,
## which is really weird for disadvantageous mutations: hence option "advantageous" below). Nonetheless, a place to begin
##      Note: I decided to start all host parameters with mut_host_ (this could be changed to reduce length of names)
#######

#' @useDynLib pevosim
#' @importFrom Rcpp evalCpp
#' @param N population size
#' @param mu per-infection mutation probability
#' @param gamma recovery rate
#' @param lbeta0 initial beta
get_mut <- function(orig_trait,mut_var,mut_mean,mut_sd,mut_host,mut_host_type,mut_host_sd_shift,mut_type="shift") {
    if (mut_type=="shift") {
      ## **direction of change different for host and pathogen and different for the two parameters (beta and gamam)
        new_trait <- orig_trait +
            rnorm(length(orig_trait),
              ifelse(mut_var == "beta",                           ## if beta is evolving:
              ifelse(mut_host == FALSE, mut_mean, mut_mean*-1),   ## in the pathogen: mut_mean is negative; in host change positive to keep disadvantageous mutations common
              ifelse(mut_host == FALSE, mut_mean*-1, mut_mean)),  ## if gamma is evolving in the pathogen mut_mean is negative; in host positive
              ifelse(mut_host == FALSE, mut_sd, mut_sd/mut_host_sd_shift)) ## allow the sd for the host to be a multiple of the sd for the pathogen (1 for equal)
        ## **if only advantageous mutations are propagated (not neutral...)
        if (mut_host == TRUE & mut_host_type == "advantageous") {
        ## ** retain the minimum of the new and original trait
        ## ** Doing it this way assumes that a host has the opportunity to counter-evolve to each strain (maybe realistic maybe not depending)
        ## ** on the assumed mechanics. Alternatively could assume a single adaptation that is equivalent to all strains
        new_trait <- pmin(new_trait, orig_trait)
        }
        if (any(is.na(new_trait))) stop("??")
        return(new_trait)
    } else stop("unknown mut_type")
}

do_mut <- function(state,mut_var,mut_link,mut_host,orig_trait,...) {
    ## **If the mutation is in the parasite
    if (mut_host == FALSE) {
    new_trait        <- get_mut(orig_trait,mut_var,mut_host,...)
    state[[mut_var]] <- c(state[[mut_var]],mut_link$linkinv(new_trait))
    state$ltraitvec  <- c(state$ltraitvec,new_trait)
    state$Ivec       <- c(state$Ivec,rep(1,length(orig_trait)))
    } else {
    ##** Else if the mutation is in the host (some adjustments to how state gets updated when the host is evolving)
    new_trait        <- get_mut(orig_trait,mut_var,mut_host,...)
    state[[mut_var]] <- mut_link$linkinv(new_trait)
    state$ltraitvec  <- new_trait
    ##** No change to I
    }
    return(state)
}

do_extinct <- function(state,mut_var,extinct) {
    state[[mut_var]] <- state[[mut_var]][-extinct]
    state$ltraitvec  <- state$ltraitvec[-extinct]
    state$Ivec       <- state$Ivec[-extinct]
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
#' @param mut_host **mutate the host as well (only on if mutation in gamma)
#' @param mut_host_type **neutral (all hosts replaced by the mutant that arises) or advantageous (only advantageous host mutations spread): either way the
#' dyanmics are one of selective sweeps
#' @param mut_host_prob **treat host mutation as a per time step probability?
#' @param mut_host_freq **If not, frequency that hosts mutate (set number of time steps (e.g. generation time or something like that)
#' @param mut_host_sd_shift **Host multiple of parasite sd
#' @param mut_host_mean_shift **Host multiple of parasite mean (not set up yet)
#' @param mut_host_mu_shift **Host multiple of parasite mean (not set up yet). On the wrong scale (muliple on prob scale) but ok for now because of tiny prob....
#' @param Ivec initial vector of infected numbers
#' @param mu mutation probability per replication
#' @param discrete (logical) run discrete-time sim?
#' @param dt time step
#' @param seed random-number seed
#' @param nt number of time steps
#' @param rptfreq reporting frequency (should divide nt)
#' @param progress draw progress bar?
#' @param debug (logical) debugging output?
#' @param debug2 **(logical) print gamma as the sim runs?
#' @param useCpp (logical) use C++ code for continuous-time models?
#' @importFrom stats make.link rnorm rbinom rmultinom rexp runif
#' @export
run_sim <- function(R0_init=2,  ## >1
                    gamma0=1/5,  ## >0
                    gamma_max=Inf,
                    N=1000,     ## integer >0
                    mu=0.01,    ## >0
                    mut_type="shift",
                    mut_mean=-1,## <0 (for sensibility)
                    mut_sd=0.5, ## >0
                    mut_var="beta",
                    mut_link=NULL,
                    mut_host=FALSE,
                    mut_host_type="neutral", ## **all host adaptations take over the pop (option: neutral) or only advantageous (option: advantageous)
                                             ## **could more cleanly use TRUE/FALSE here, but wanted to leave options open for other alternatives
                    mut_host_prob=TRUE,
                    mut_host_freq="prob",
                    mut_host_sd_shift=1,     ## **1 for identical sd to the parasite (mean assumed to be the same)
                 #  mut_host_mean_shift=1,   ## **1 for identical sd to the parasite (mean assumed to be the same); not included yet
                    mut_host_mu_shift=1,
                    Ivec=NULL,
                    dt=1,
                    nt=100000,
                    rptfreq=max(nt/500,1), ## divides nt?
                    discrete=TRUE, ## discrete-time?
                    seed=NULL,
                    progress=FALSE,
                    debug=FALSE,
                    debug2=FALSE,
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
      ## **Positive with log link can push gamma overboard when mutation in gamma is on average disadvantageous.
      ## **Made to logit link for now, not sure...
        if (is.null(mut_link)) mut_link <- make.link("logit") #mut_link <- make.link("log")
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
    ## **decided I wanted to track more stuff
    res <- as.data.frame(matrix(NA,nrow=nrpt,ncol=8,
                                dimnames=list(
                                    NULL
                                  , c(
                                      "time"
                                    , "S"
                                    , "I"
                                    , paste0(c("mean_l", "sd_l"), mut_var)
                                    , "num_strains"
                                    , ifelse(mut_var == "beta", "gamma", "beta") ##** list the non-evolving param
                                    , "pop_size"))))

    t_tot <- 0  ## for continuous model

    for (i in 1:nrpt) {
        ## cat("time ",i,"\n")
        if (discrete) {
            for (j in 1:rptfreq) {
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

                ## **new host mutation, arising every time step with a given prob ("prob") or after X generations (a numeric value)
                if (mut_host_freq == "prob") {
                ## **assume probability of host mutation to be some factor of parasite mutation
                mutated_host <- rbinom(N,1,prob=mu/mut_host_mu_shift)
                ## **because host types aren't being tracked, for now just track if any mutations arise or not
                ## **was thinking about taking X number of mutations and choosing most advantageous from them but that seems to be a
                ## **very half-assed improvement. So just mechanics for now, worry about this stuff when the model changes
                mutated_host <- ifelse(sum(mutated_host) > 0, 1, 0)
                } else if (class(mut_host_freq) == "numeric") {
                ## **check if host reproduction has occurred. Prior to setting up different classes of hosts, this assumes instantaneous selective sweep
                  ## **which is pretty rediculuous, but a starting place
                if ((i / mut_host_freq) %% 1 == 0) {
                mutated_host <- 1
                } else {
                mutated_host <- 0
                }
                }

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
                ## **mut parasite first. Another choice of order that may matter (especically if the model gets more complicated)
                if (tot_mut>0) {
                    state <- do_mut(state,
                                    mut_var=mut_var,
                                    mut_link,
                                    orig_trait=rep(state$ltraitvec,mutated),
                                    mut_mean=mut_mean,
                                    mut_sd=mut_sd,
                                    mut_host=FALSE,
                                    mut_host_type=mut_host_type,
                                    mut_host_sd_shift=mut_host_sd_shift,
                                    mut_type=mut_type)
                }
                    ## **host mutation if there was one
                    if (mut_host == TRUE & mutated_host == 1) {
                    state <- do_mut(state,
                                    mut_var=mut_var,
                                    mut_link,
                                    orig_trait=rep(state$ltraitvec,1),
                                    mut_mean=mut_mean,
                                    mut_sd=mut_sd,
                                    mut_host=TRUE,
                                    mut_host_type=mut_host_type,
                                    mut_host_sd_shift=mut_host_sd_shift,
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

                ##** check progress of evolution
                if (debug2 == TRUE) {
                print(paste(j, state$gamma, sep = "  "))
                }

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
        ## **Decided I wanted to track some more stuff
        num_strains <- length(state$Ivec)
        if (progress) cat(".")
        res[i,] <- c(i*rptfreq,state$S,I_tot,ltrait_mean,ltrait_sd,num_strains,
          ifelse(mut_var == "beta", state$gamma, state$beta),sum(state$S, state$I))
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
