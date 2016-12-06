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

## For efficiency (??), we keep track of parameters on both
## the constrained and unconstrained scales ...
## state:
## state$beta$constr, state$beta$unconstr,
##    state$gamma$constr, state$gamma$unconstr

## FIXME: only works on single mutated strains at present
## allow mut_var to differ (randomly) among strains mutating
##   at each step ... ?

do_mut2 <- function(state,mut_var,mut_link,strain,mut_mean,mut_sd) {

    ## retrieve original unconstr and constr trait values
    ## vectorized over which trait(s) and strain(s) are mutating
    sMap <- function(...) unlist(Map(...))
    nmut <- length(strain)
    for (tt in names(state$traits)) {
        orig_trait <- state$traits[[tt]]$unconstr[strain]
        orig_trait_constr <- state$traits[[tt]]$constr[strain]
        new_trait <- get_mut(orig_trait,mut_mean,mut_sd)
        new_trait_constr <- sMap(function(f,o) f$linkinv(o), mut_link, new_trait)
        ## set up all traits for new strain
        ttu <- ifelse(tt==mut_var,new_trait,orig_trait)
        ttc <- ifelse(tt==mut_var,new_trait_constr,orig_trait_constr)
        ## append
        state$traits[[tt]]$unconstr <- c(state$traits[[tt]]$unconstr,ttu)
        state$traits[[tt]]$constr <- c(state$traits[[tt]]$constr,ttc)
    }
    state$Ivec <- c(state$Ivec,rep(1,nmut))
    return(state)
}

##' do mutation, allowing for mutation in either trait
##' @param state state list
##' @param mut_var which variable mutates
##' @param mut_link link/inv link function for traits
##' @param strain which strain(s) are mutating?
##' @param mut_mean mutation bias
##' @param mut_sd mutation SD
##' @importFrom methods is
##' @export
do_mut3 <- function(state,mut_var,mut_link,strain,mut_mean,mut_sd) {
    if (mut_var=="both") {
        ## for now, equal probability of beta vs gamma mutation
        mut_var <- ifelse(runif(length(strain))<0.5,
                          "beta","gamma")
    }
    mm <- if (length(mut_mean)==1) mut_mean else mut_mean[mut_var]
    ms <- if (length(mut_sd)==1) mut_sd else mut_sd[mut_var]
    ml <- if (is(mut_link,"link-glm")) list(mut_link) else mut_link[mut_var]
    if ((nstrain <- length(strain))>1) {
        mut_var <- rep(mut_var,length.out=nstrain)
        mm <- rep(mm,length.out=nstrain)
        ms <- rep(ms,length.out=nstrain)
        ml <- replicate(nstrain,ml)
    }
    ## switch sign of mutation bias for gamma
    mm <- ifelse(mut_var=="gamma",-mm,mm)
    state <- do_mut2(state,
                     mut_var,
                     ml,
                     strain,
                     mut_mean=mm,
                     mut_sd=ms)
    return(state)
}

##' process extinction (remove zero-density strains)
##' @param state state list
##' @param strain indices of strain(s) that went extinct
##' @export
do_extinct <- function(state,strain) {
    for (v in names(state$traits)) {
        for (cc in c("constr","unconstr")) {
            state$traits[[v]][[cc]] <- state$traits[[v]][[cc]][-strain]
        }
    }
    state$Ivec <- state$Ivec[-strain]
    return(state)
}

##' Logit transform on (min,max) rather than (0,1);
##' reduces to log transform if max value is unbounded
##' @param minval minimum value
##' @param maxval maximum value
##' @importFrom stats plogis qlogis
##' @export
multlogit <- function(minval=0,maxval=1) {
    if (maxval<Inf) {
        delta <- maxval-minval
        linkfun <- function(mu) qlogis((mu-minval)/delta)
        linkinv <- function(eta) plogis(eta)*delta+minval
    } else {
        linkfun <- function(mu) log(mu-minval)
        linkinv <- function(eta) exp(eta)+minval
    }
    structure(list(linkfun = linkfun, linkinv = linkinv, name = "multlogit"),
              class = "link-glm")
}


##' make link functions for both parameters, based on tradeoff
##' @param tradeoff_fun string: "powerlaw" or "none"
##' @param tradeoff_pars listof parameters for tradeoff curve (g1, g2 = scale and curvature of power law)
##' @param hazard parameterize events by hazard?
##' @param gamma value of gamma
##' @param beta value of beta
##' @param discrete (logical)
##' @export
make_bilink <- function(tradeoff_fun,tradeoff_pars,
                        gamma,beta,discrete=FALSE,
                        hazard=FALSE) {
    if (tradeoff_fun=="powerlaw") {
        maxbeta <- with(as.list(tradeoff_pars),g1*gamma^(1/g2))
        ## mingamma <- with(as.list(tradeoff_pars),(beta/g1)^g2)
        return(list(beta=multlogit(0,maxbeta),
                    gamma=multlogit(tradeoff_pars["mingamma"],Inf)))
    } else if (tradeoff_fun=="none") {
        if (discrete && !hazard) {
            ## version 0: beta constrained 0-1
            return(list(beta=make.link("logit"),
                        gamma=make.link("logit")))
        } else {
            ## unconstrained; beta, gamma can be any positive value
            return(list(beta=make.link("log"),
                        gamma=make.link("log")))
        }
    } else stop("unknown tradeoff function")
}

make_state <- function(gamma,beta,S,
                       Ivec=rep(1,n),
                       tradeoff_fun="powerlaw",
                       tradeoff_pars=c(g1=1,g2=2),
                       mut_links=NULL) {
    n <- max(length(gamma),length(beta))
    gamma <- rep(gamma,length.out=n)
    beta <- rep(beta,length.out=n)
    return(list(traits=list(beta=list(constr=beta,
                                      unconstr=mut_links$beta$linkfun(beta)),
                            gamma=list(constr=gamma,
                                      unconstr=mut_links$gamma$linkfun(gamma))),
                Ivec=Ivec,
                S=S))
}

check_state <- function (state,N=NA) {
    ok <- TRUE
    n <- length(state$traits$beta$constr)
    ok <- ok && length(state$traits$beta$constr)==n &&
        length(state$traits$beta$unconstr)==n &&
        length(state$traits$gamma$constr)==n && 
        length(state$traits$gamma$unconstr)==n && 
        length(state$Ivec)==n
    if (!is.na(N)) {
        ok <- ok && sum(state$Ivec)+state$S==N
    }
    return(ok)
}

#' Evolutionary simulations
#' 
#' @param inits list
#' \itemize{
#' \item R0 (starting value of R0)
#' \item gamma (recovery rate)
#' \item Ivec (initial infectives)
#' }
#' @param mod_inits list of initial values to modify
#' @param params list: default values
#' \itemize{
#' \item N=1000 (population size)
#' \item mu (mutation probability per replication)
#' \item mut_type mutation model (only "shift" implemented)
#' \item mut_mean mean value of mutations
#' \item mut_sd mutational standard deviations
#' \item mut_var which variable(s) mutate?
#' \item mut_link link function for mutation
#' \item tradeoff_pars
#' }
#' @param mod_params list of parameters to modify
## FIXME: do this by storing default value list in environment?? 
#' @param discrete (logical) run discrete-time sim?
#' @param tradeoff_fun tradeoff function, "powerlaw" or "none"
#' @param seed random-number seed
#' @param nt number of time steps
#' @param dt within-step granularity
#' @param hazard parameterize events by hazard? (allows time-step cutting in discrete-time models)
#' @param rptfreq reporting frequency (should divide nt)
#' @param progress draw progress bar?
#' @param debug (logical) debugging output?
#' @param useCpp (logical) use C++ code for continuous-time models?
#' @importFrom stats make.link rnorm rbinom rmultinom rexp runif
#' @importFrom utils modifyList
#' @export
run_sim <- function(inits=list(R0=2, ## >1
                              gamma=1/5,
                              Ivec=NULL),
                    mod_inits=list(),
                    params=list(N=1000,     ## integer >0
                                mu=0.01,    ## >0
                                mut_type="shift",
                                mut_mean=-1,## <0 (for sensibility)
                                mut_sd=0.5, ## >0
                                mut_var="beta",
                                mut_link=NULL,
                                tradeoff_pars=NULL),
                    mod_params=list(),
                    tradeoff_fun="none",
                    nt=100000,
                    dt=1,
                    rptfreq=max(nt/500,1), ## divides nt?
                    discrete=TRUE, ## discrete-time?
                    hazard=FALSE,
                    seed=NULL,
                    progress=FALSE,
                    debug=FALSE,
                    useCpp=FALSE) {

    if (discrete && !hazard && dt!=1) {
        warning("probs not adjusted for dt!=1 ...")
    }
    if (length(mod_params)>0) {
        params <- modifyList(params,mod_params)
    }
    if (length(mod_inits)>0) {
        inits <- modifyList(inits,mod_inits)
    }

    if (round(params$N)!=params$N) {
        warning("rounding N")
        params$N <- round(params$N)
    }
    if (params$mut_mean>0) {
        warning("positive mutation bias")
    }
    stopifnot(inits$R0>0,
              inits$gamma>0,
              params$N>0,
              params$mu>0,
              params$mut_sd>0,
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

    if (is.null(inits$Ivec)) {
        ## start at equilibrium I ...
        inits$Ivec <- max(1,round(params$N*(1-1/inits$R0)))
    }
    
    ## R0 = beta*N/gamma
    ## want

    ## raw scale, but should be small, rate ~ prob
    beta0 <- inits$R0*inits$gamma/params$N
    ## prod((1-beta0)^inits$Ivec)
    ##
    if (is.null(mut_link <- params$mut_link)) {
        mut_link <- make_bilink(tradeoff_fun,
                                params$tradeoff_pars,inits$gamma,beta0,
                                discrete=discrete,
                                hazard=hazard)
    }

    ## set up initial trait vector
    state <- make_state(inits$gamma,beta0,
                        params$N-sum(inits$Ivec),
                        Ivec=inits$Ivec,
                        tradeoff_fun=tradeoff_fun,
                        tradeoff_pars=params$tradeoff_pars,
                        mut_links=mut_link)
     

    dfun("init")
    
    nrpt <- nt %/% rptfreq
    ## not used: intended for accumulating information about a binned
    ##  distribution over time
    ## nq <- 9
    ## qvec <- seq(0,1,length=nq+2)[-c(1,nq+2)]
    tnames <- c(outer(c("mean_l","sd_l"),names(state$traits),paste0))
    res <- as.data.frame(matrix(NA,nrow=nrpt,ncol=8,
                                dimnames=list(NULL,c("time","S","I",
                                                    tnames,
                                                    "n_events"))))

    t_tot <- 0  ## for continuous model

    t_steps <- round(1/dt)

    for (i in 1:nrpt) {
        n_events <- 0
        ## cat("time ",i,"\n")
        if (discrete) {
            for (j in 1:rptfreq) {
                for (tstp in 1:t_steps) {
                ## if (i==43 && j==20) browser()
                nstrain <- length(state$Ivec)
                beta <- state$traits$beta$constr
                gamma <- state$traits$gamma$constr
                ## cat("betavec:",betavec,"\n")
                ## cat("Ivec:",Ivec,"\n")
                ## prob of escaping infection completely
                infprob <- if (hazard) 1-exp(-beta*dt) else beta
                uninf <- rbinom(1,size=state$S,
                                prob=prod((1-infprob)^state$Ivec))
                n_events <- n_events + (state$S-uninf)
                ## division of new infectives among strains
                ## 'prob' is internally normalized
                if (state$S-uninf>0) {
                  newinf <- drop(rmultinom(1,size=state$S-uninf,
                                         prob=state$Ivec*infprob))
                } else newinf <- rep(0,length(state$Ivec))
                stopifnot(uninf+sum(newinf)==state$S)
                recprob <- if (hazard) 1-exp(-gamma*dt) else gamma
                recover <- rbinom(length(state$Ivec),
                                  size=state$Ivec,prob=recprob)
                n_events <- n_events+sum(recover)
                ## fraction of new infections -> mutation
                mutated <- rbinom(length(newinf),size=newinf,prob=params$mu)
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
                    state <- do_mut3(state,
                                     params$mut_var,
                                     mut_link,
                                     strain=rep(1:nstrain,mutated),
                                     mut_mean=params$mut_mean,
                                     mut_sd=params$mut_sd)
                }
                if (discrete && !hazard &&
                    any(state$traits$beta$constr>1)) stop("beta >1")
                if (all(state$Ivec==0)) {
                    message(sprintf("system went extinct prematurely (t=%d)",i))
                    break
                }
                if (length(extinct <- which(state$Ivec==0))>0) {
                    state <- do_extinct(state,extinct)
                }
                dfun("after mutation")
                if (length(state$Ivec)>0) {
                    ## avoid X+numeric(0)==numeric(0) problem ...
                    state$S <- state$S + sum(recover) - sum(newinf)
                }
                stopifnot(length(state$S)==1)
                stopifnot(sum(state$Ivec)+state$S == params$N)
                } ## intermed time steps
            }  ## rptfreq time steps
        } else  ## end discrete case
        {
            if (useCpp) {
                ## cat("useCpp: state: \n")
                ## print(state)
                run_stepC(state,
                          (i-1)*rptfreq,
                          i*rptfreq,
                params=list(mu=params$mu,
                            mut_sd=params$mut_sd,
                            mut_mean=params$mut_mean,
                            mut_var=params$mut_var,
                            mut_link=""), ## FIXME
                debug=debug)
            } else { ## use R
                while (t_tot<(i+1)*rptfreq) {
                    nstrain <- length(state$Ivec)
                    rates <- get_rates(state)
                    if (debug) cat(t_tot,rates,"\n")
                    t_tot <- t_tot + rexp(1,sum(rates))
                    w <- sample(length(rates),size=1,prob=rates)
                    event <-  ((w-1) %/% nstrain) + 1
                    strain <- ((w-1) %% nstrain) + 1
                    if (event==1) { ## infection
                        if (runif(1)<params$mu) {
                            state <- do_mut3(state,
                                             params$mut_var,
                                             mut_link,
                                             strain,
                                             mut_mean=params$mut_mean,
                                             mut_sd=params$mut_sd)
                        } else {
                            state$Ivec[strain] <- state$Ivec[strain]+1
                        }
                        state$S <- state$S - 1
                    } ## infection

                    if (event==2) { ## recovery
                        state$Ivec[strain] <- state$Ivec[strain]-1
                        state$S <- state$S+1
                        if (state$Ivec[strain]==0) {
                            state <- do_extinct(state,strain)
                        }
                    }
                    n_events <- n_events+1
                    if (sum(state$Ivec)==0) break
                    stopifnot(length(state$S)==1)
                    stopifnot(sum(state$Ivec)+state$S == params$N)
                } ## loop until end of time period
            } ## use R, not Rcpp
        } ## continuous time
        ## summary statistics
        res[i,"time"] <- i*rptfreq
        res[i,"S"] <- state$S
        res[i,"I"] <- I_tot <- sum(state$Ivec)
        for (tt in names(state$traits)) {
            res[i,paste0("mean_l",tt)] <- tmean <- 
                sum(state$Ivec*state$traits[[tt]]$unconstr)/I_tot
            res[i,paste0("sd_l",tt)] <-
                sqrt(sum(state$Ivec*(state$traits[[tt]]$unconstr-tmean)^2)/I_tot)
        }
        res[i,"n_events"] <- n_events
        if (progress) cat(".")
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
get_rates <- function(state) {
    inf_rates <- state$traits$beta$constr*state$Ivec*state$S
    recover_rates <- state$traits$gamma$constr*state$Ivec
    return(c(inf_rates,recover_rates))
}
