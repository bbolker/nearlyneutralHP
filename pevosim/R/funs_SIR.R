#######
## Expand the work in funs.R, which considers only parasite evolution (with host evolution as a true afterthought)
## to include host evolution. Requires a switch to an SIR model or at the very least an SIS model with host death
## as well as host reproduction
#######

#######
## Search ^^ for planned changes
#######

## ^^ 1 Establish _how_ the model is going to work. What does resistance and tolerance do in hosts for the variety of parasite strains?
## ^^         A tolerance and or resistant host infected with a two different pathogen strains results in?
## ^^   1.1 For a first pass assume hosts reduce beta or increase recovery as a multiplicative change on the logit scale
## ^^         This is only resistance (but the structure of the model won't change from here so a fine place to begin). Tolerance can come later
## ^^ 2 Set up the structure of the model with a first pass such that there are no costs to host or parasite
## ^^   2.1 Convert to having a vector of S
## ^^   2.2 Which is going to require a matrix of I
## ^^   2.3 Convert to an SIS model
## ^^ 3 Add parasite evolution along a tradeoff
## ^^ 4 Add host evolution with simple costs


######
## General notes
######

## ^^I am not very happy with all of the matrix(..., ncol = X, nrow = Y) floating around. Would like to clean this up
## ^^Some checks have been commented out. They need to be reframed and added back in. (#!)
## ^^Assuming mutated hosts come from only S class will lead to few mutated hosts as the sim progresses as beta climbs to 1

## ^^state$ltraitvec can remain a vector, which defines the intrinsic strategy of each parasite strain
## ^^state$hrtraitvec can be a new vector that gives host investment in resistance
## ^^These two vectors give rise to a matrix of the intersection
## ^^Newly infected S (for now we won't assume S resist getting infected???) are determined using a
##    ^^multinomial draw over the matrix of I


#' @useDynLib pevosim
#' @importFrom Rcpp evalCpp
#' @param N population size
#' @param mu per-infection mutation probability
#' @param gamma recovery rate
#' @param lbeta0 initial beta
get_mut <- function(orig_trait, mut_var, mut_mean, mut_sd, mut_host, mut_host_type, mut_host_sd_shift, mut_type = "shift") {
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

do_mut <- function(state, mut_var, mut_link, mut_host, orig_trait, ...) {
    ## **If the mutation is in the parasite
    if (mut_host == FALSE) {
    new_trait        <- get_mut(orig_trait,mut_var,mut_host,...)
    state$ltraitvec  <- c(state$ltraitvec,new_trait)
    } else {
    ## **else if the mutation is in the host (some adjustments to how state gets updated when the host is evolving)
    ## ^^new vector for host traits
    new_trait        <- get_mut(orig_trait,mut_var,mut_host,...)
    state$hrtraitvec <- c(state$hrtraitvec,new_trait)
    }
    return(state)
}

update_mut <- function(state, mut_link, mutated, mutated_host, mut_var) {
   ## ^^Assumes mutated hosts can't get infected in the same time step. Birth or something...

   ## ^^Update all trait combinations
   state[[mut_var]] <- t(mut_link$linkinv(outer(state$ltraitvec, state$hrtraitvec,  "/")))

   ## ^^Updating alpha will follow the tradeoff curve eventually, but for now just propagate 1s
   state$alpha <- matrix(data = state$alpha[1,1], nrow = nrow(state$beta), ncol = ncol(state$beta))

   ## ^^Update Infected matrix with the new strain, maintaining which S class received that mutation.
    ## ^^For each mutated strain, Imat gets a new column with a single 1, in the row in which the mutation occurred
   if (sum(mutated) > 0) {
     ## ^^probably a rediculous way to do this, but I cant figure out a better way right now
     ## ^^first X columns of 0s
     new_mutes   <- matrix(data = 0, nrow = nrow(state$Imat), ncol = sum(mutated))
     num_mutes   <- rowSums(mutated)
     which_hosts <- rep(which(num_mutes > 0), num_mutes[which(num_mutes > 0)])
     for (z in seq_along(which_hosts)) {
      new_mutes[which_hosts[z], z] <- 1
     }
    state$Imat <- cbind(state$Imat, new_mutes)
   }

   if (sum(mutated_host) > 0) {
   ## ^^Update Svec (mutated_host is a vector that tracks which host mutated)
   state$Svec       <- cbind(state$Svec, matrix(rep(1, sum(mutated_host)), nrow = 1))
   ## ^^Add new rows with 0s for the new S genotypes
   state$Imat       <- rbind(state$Imat, matrix(data = rep(0, sum(mutated_host)*ncol(state$Imat)), ncol = ncol(state$Imat), nrow = sum(mutated_host)))
   }

   return(state)
}

do_extinct <- function(state,mut_var,extinct,parasite) {
  ## ^^If a parasite strain has gone extinct remove a column
  if (parasite == TRUE) {
    state[[mut_var]] <- state[[mut_var]][,-extinct, drop = FALSE]
  ## ^^Again, will eventually need to clean this up
    state$alpha      <- state$alpha[,-extinct, drop = FALSE]
    state$ltraitvec  <- state$ltraitvec[-extinct, drop = FALSE]
    state$Imat       <- state$Imat[,-extinct, drop = FALSE]
 ## ^^If a host strain has gone extinct remove a row
  } else {
    state[[mut_var]] <- state[[mut_var]][-extinct, , drop = FALSE]
    state$alpha      <- state$alpha[-extinct, , drop = FALSE]
    state$hrtraitvec <- state$hrtraitvec[-extinct, drop = FALSE]
    state$Imat       <- state$Imat[-extinct, , drop = FALSE]
    state$Svec       <- state$Svec[-extinct, drop = FALSE]
  }
    return(state)
}

## ^^This function is borderline hideous but I am almost positive it does what I want it to do. All ears for a better way to write this
get_inf <- function (Svec, uninf, Imat, beta) {

matrix(
data =
    apply(
  Svec - uninf
, 2
, function (x) colSums(
matrix(
data = rmultinom(
   1
 , size = x
 , prob = beta
   )
, ncol = ncol(beta)
, nrow = nrow(beta)
      )
    )
  )
, ncol = ncol(beta)
, nrow = nrow(beta)
, byrow = TRUE
)

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
#' @param host_dyn_only With SIR model check that host population doesn't go crazy in the absence of pathogen
#' @param R0_init starting value of R0
#' @param gamma0 recovery rate
#' @param alpha0 intitial parasite virulence
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
#' @param start_res ## ^^starting host resistance value
#' @param start_tol ## ^^starting host tolerance value
#' @param Imat initial vector of infected numbers
#' @param mu mutation probability per replication
#'
#' @param b host birth rate (for now assume only S births, but this may have to change)
#' @param d host death rate for S
#'
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

#' ## ^^New parameters will include
#' @param many_tradeoff_curves host and parasite cost and benefit functions
#'

run_sim <- function(R0_init=2,  ## >1
                    gamma0=1/5,  ## >0
                    gamma_max=Inf,
                    alpha0=1,
                    b=0.1,
                    d=0.1,
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
                    mut_host_mu_shift=2,
                    start_res=1,   ## ^^starting host resistance value
                    start_tol=1,   ## ^^starting host tolerance value
                    Imat=NULL,
                    host_dyn_only=FALSE,
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

    if (is.null(Imat)) {
        ## start at equilibrium I ...
        Imat <- max(1,round(N*(1-1/R0_init)))
        Imat <- as.matrix(Imat)
    }
    Svec <- N-Imat

    ## R0 = beta*N/gamma
    beta0 <- R0_init*gamma0/N

    ## ^^Assume some start for alpha, for now assume 1 (set init_alpha to 1), build in tradeoff curve later
    alpha0 <- alpha0

    ## set up initial trait vector
    if (mut_var=="beta") {
        if (is.null(mut_link)) {
            mut_link <- if (discrete) make.link("logit") else make.link("log")
        }

      ## ^^What to do with combined parasite and host traits?
        ltraitvec <- mut_link$linkfun(beta0)
    } else {
      ## **Positive with log link can push gamma overboard when mutation in gamma is on average disadvantageous.
      ## **Made to logit link for now, not sure...
        if (is.null(mut_link)) mut_link <- make.link("logit") #mut_link <- make.link("log")
        ltraitvec <- mut_link$linkfun(gamma0)
    }

    ## ^^Initial trait vec for the host
    hrtraitvec <- rep(start_res, length(Svec))
    httraitvec <- rep(start_tol, length(Svec))

    ## parameters structure (parallel vectors), so these can
    ## be modified via function and passed back ...
    state <- list(
      beta       = as.matrix(beta0)  ## ^^beta will now be dependent on both host and parasite phenotype
    , gamma      = as.matrix(gamma0) ## ^^as will gamma (make matrix after, because ltraitvec needs to remain a vector)
    , alpha      = as.matrix(alpha0) ## ^^new parameter. Virulence of parasite
    , ltraitvec  = ltraitvec
    , hrtraitvec = hrtraitvec
    , httraitvec = httraitvec
    , Imat       = Imat
    , Svec       = Svec)

    ## ^^for no infection (to check host dynamics in the absence of infection)
    ## ^^Maybe not the _most_ iffecient place/way to do this, but I think makes the least clutter
    if (host_dyn_only == TRUE) {
#      state$beta     <- 0
      state$Imat[1,1] <- 0
      state$Svec[1,1] <- N
#      mu             <- 0
    }

    dfun("init")

    nrpt <- nt %/% rptfreq
    ## not used: intended for accumulating information about a binned
    ##  distribution over time
    ## nq <- 9
    ## qvec <- seq(0,1,length=nq+2)[-c(1,nq+2)]
    ## **decided I wanted to track more stuff
    ## ^^add tracking host responses to this data frame
    res <- as.data.frame(matrix(NA, nrow = nrpt, ncol = 11,
                                dimnames = list(
                                    NULL
                                  , c(
                                      "time"
                                    , "num_S"
                                    , "num_S_strains"
                                    , "num_I"
                                    , "num_I_strains"
                                    , "pop_size"
                                    , paste0(c("mean_hl", "sd_hl"), mut_var)
                                    , paste0(c("mean_pl", "sd_pl"), mut_var)
                                    , ifelse(mut_var == "beta", "gamma", "beta") ##** list the non-evolving param
                                    ))))

    t_tot <- 0  ## for continuous model

    for (i in 1:nrpt) {
        ## cat("time ",i,"\n")
        if (discrete) {
            for (j in 1:rptfreq) {
                ## cat("betavec:",betavec,"\n")
                ## cat("Imat:",Imat,"\n")

#                print(j)

                ## ^^[Step 1]: Birth ^^ ## Not accessible to death or infection until the next time step.
                ## ^^ugly form of density dependence here for now keeping per-capita birth equal to the total death rate
                  ## of both S and I
                birth_rate <- ((sum(state$Imat) - (sum(state$Imat) * (1 - d) * (1 - d * state$alpha))) + (sum(state$Svec) * d)) / sum(state$Svec)
                birth <- matrix(rbinom(length(state$Svec),size=state$Svec,prob=birth_rate), nrow = 1)

                ## ^^[Step 2]: Death of S ^^ ##
                deathS     <- matrix(rbinom(length(state$Svec),size=state$Svec,prob=d), nrow = 1)
                state$Svec <- state$Svec - deathS

                ## ^^[Step 3]: Infection ^^ ##
                ## prob of escaping infection completely
                uninf <- matrix(rbinom(length(state$Svec),size=state$Svec,prob=prod((1-state$beta)^state$Imat)), nrow = 1)
                ## division of new infectives among strains
                ## 'prob' is internally normalized
                ## ^^See function. A bit ugly because of apply and matrices needing a transpose...
                ## ^^would like to use ifelse here, but annoyingly ifelse doesn't seem to want to pass a matrix

                newinf <- get_inf(state$Svec, uninf, state$Imat, beta = state$beta)
                stopifnot(sum(rowSums(newinf) + uninf) == sum(state$Svec))

                ## ^^[Step 4]: Death of I ^^ ##
                 ## ^^Natural death + parasite induced death
                deathI     <- matrix(data = rbinom(length(c(state$Imat)),size=c(state$Imat),prob=d)
                  , ncol = ncol(state$Imat), nrow = nrow(state$Imat))
                state$Imat <- state$Imat - deathI
                ## ^^here alpha is treated as some multiple of d
                deathI     <- matrix(data = rbinom(length(c(state$Imat)),size=c(state$Imat),prob=c(state$alpha)*d)
                  , ncol = ncol(state$Imat), nrow = nrow(state$Imat))
                state$Imat <- state$Imat - deathI
                stopifnot(((all(state$Imat) >= 0) | all(state$Svec) >= 0))

                ## ^^[Step 5]: Recovery of I ^^ ##
                recover <- matrix(data = rbinom(length(c(state$Imat)),size=c(state$Imat),prob=c(state$gamma))
                  , ncol = ncol(state$Imat), nrow = nrow(state$Imat))

                ## ^^[Step 6]: Mutation of new infections ^^ ##
                ## fraction of new infections -> mutation
                #mutated <- rbinom(ncol(newinf), size=colSums(newinf), prob=mu)
                ## ^^Would like to use ifelse, but I hate its behavior
               if (host_dyn_only == FALSE) {
                  mutated <- matrix(rbinom(newinf, size=newinf, prob=mu), nrow = nrow(newinf), ncol = ncol(newinf))
               } else {
                  mutated <- matrix(rbinom(nrow(newinf), size=newinf, prob=mu), nrow = nrow(newinf), ncol = ncol(newinf))
               }

                ## **new host mutation, arising every time step with a given prob ("prob") or after X generations (a numeric value)

                ## ^^[Step 7]: Mutation of new births ^^ ##
                if (mut_host_freq == "prob") {
                ## **assume probability of mutation is just in S (as if S hosts are the only hosts reproducing -- a common assumption)
                mutated_host <- rbinom(length(birth),size=birth, prob = mu/mut_host_mu_shift)
                birth        <- birth - mutated_host
                ## **because host types aren't being tracked, for now just track if any mutations arise or not
                ## **was thinking about taking X number of mutations and choosing most advantageous from them but that seems to be a
                ## **very half-assed improvement. So just mechanics for now, worry about this stuff when the model changes
                } else if (class(mut_host_freq) == "numeric") {
                ## **check if host reproduction has occurred. Prior to setting up different classes of hosts, this assumes instantaneous selective sweep
                  ## **which is pretty rediculuous, but a starting place
                if ((i / mut_host_freq) %% 1 == 0) {
                ## ^^This no longer makes sense for a vector of S
                mutated_host <- rep(1, length(state$Svec))
                } else {
                mutated_host <- 0
                }
                }

              ## ^^Some debugging stuff
#                state_check <- state
#                state_check$j <- j
#                assign("state_check",state_check,.GlobalEnv)
#                assign("uninf_check",uninf,.GlobalEnv)

#!              stopifnot(length(recover)==length(mutated))
#!              stopifnot(length(newinf)==length(state$Imat))
                ## update infectives
                ## ^^I am not enirely sure where I went wrong with my code, but I kinda hate having to t()...
#                if (host_dyn_only == FALSE) {
#                print(recover); print(newinf)
                state$Imat <- state$Imat - recover + newinf - mutated
#                } else {
#                state$Imat <- state$Imat - recover + t(newinf) - t(mutated)
#                }
                ## cat("I,uninf, unmutated, mutated, recovered",
                ## c(sum(Imat),uninf,sum(newinf-mutated),sum(mutated),sum(recover)),
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
                                    orig_trait=rep(state$ltraitvec,colSums(mutated)), ## ^^Just care about intrinsic nature of a strain
                                    mut_mean=mut_mean,
                                    mut_sd=mut_sd,
                                    mut_host=FALSE,
                                    mut_host_type=mut_host_type,
                                    mut_host_sd_shift=mut_host_sd_shift,
                                    mut_type=mut_type)
                }
                    ## **host mutation if there was one
                    ## ^^Mechanics of host mutation will only need to change a bit (similar mechanics to parasite vector: rbinom(1,S_i,mu))

                    if (mut_host == TRUE & sum(mutated_host) > 0) {
                    state <- do_mut(state,
                                    mut_var=mut_var,
                                    mut_link,
                                    orig_trait=rep(state$hrtraitvec,mutated_host), ## ^^For now assume mutated hosts come from S class
                                    mut_mean=mut_mean,
                                    mut_sd=mut_sd,
                                    mut_host=TRUE,
                                    mut_host_type=mut_host_type,
                                    mut_host_sd_shift=mut_host_sd_shift,
                                    mut_type=mut_type)
                    }

                ## ^^update Svec with infections and recoveries prior to the mutations (and host birth)
                if (length(state$Imat)>0) {
                    ## avoid X+numeric(0)==numeric(0) problem ...
#                if (host_dyn_only == FALSE) {
                state$Svec <- state$Svec + rowSums(recover) - rowSums(newinf) + birth
#                } else {
#                state$Svec <- state$Svec + t(rowSums(recover)) - t(rowSums(t(newinf))) + birth
#                }
                }

                ## Update state after host and parasite mutations
                if (tot_mut > 0 | sum(mutated_host) > 0) {
                state <- update_mut(
                    state        = state
                  , mut_link     = mut_link
                  , mutated      = mutated
                  , mutated_host = mutated_host
                  , mut_var      = mut_var)
                }

                if (all(state$Imat==0) & host_dyn_only == FALSE) {
                    message(sprintf("system went extinct prematurely (t=%d)",i))
                    break
                }
                if (all(state$Svec==0) & host_dyn_only == TRUE) {
                    message(sprintf("system went extinct prematurely (t=%d)",i,j))
                    break
                }

                ## ^Look over a column for an extinct strain
                extinct_p <- which(colSums(state$Imat) == 0)
                extinct_h <- which(rowSums(state$Imat) == 0 & state$Svec == 0)
                if (length(extinct_p)>0 & host_dyn_only == FALSE) {
                    state <- do_extinct(state,mut_var,extinct = extinct_p, parasite = TRUE)
                }
                if (length(extinct_h)>0) {
                    state <- do_extinct(state,mut_var,extinct = extinct_h, parasite = FALSE)
                }
                dfun("after mutation")

     #!         stopifnot(length(state$Svec)==1)
     #!         stopifnot(sum(state$Imat)+sum(state$Svec) == N)

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
                            state$Imat[strain] <- state$Imat[strain]+1
                        }
                        state$S <- state$S - 1
                    } ## infection

                    if (event==2) { ## recovery
                        state$Imat[strain] <- state$Imat[strain]-1
                        state$S <- state$S+1
                        if (state$Imat[strain]==0) {
                            state <- do_extinct(state,mut_var,strain)
                        }
                    }
                    stopifnot(length(state$S)==1)
                    stopifnot(sum(state$Imat)+state$S == N)
                } ## loop until end of time period
            } ## use R, not Rcpp
        } ## continuous time
        ## summary statistics
        I_tot        <- ncol(state$Imat)
        num_I        <- sum(state$Imat)
        ltrait_mean  <- sum(colSums(state$Imat)*state$ltraitvec)/num_I
        ltrait_sd    <- sqrt(sum(colSums(state$Imat)*(state$ltraitvec-ltrait_mean)^2)/num_I)
        num_S        <- sum(state$Svec)
        S_tot        <- length(state$Svec)
        lhtrait_mean <- sum(state$Svec*state$hrtraitvec)/num_S
        lhtrait_sd   <- sqrt(sum(state$Svec*(state$hrtraitvec-lhtrait_mean)^2)/num_S)
        pop_size     <- num_I + num_S

        ## **Decided I wanted to track some more stuff
        num_strains <- length(state$Imat)
        if (progress) cat(".")
        res[i,] <- c(
          i*rptfreq
        , num_S
        , S_tot
        , num_I
        , I_tot
        , pop_size
        , lhtrait_mean
        , lhtrait_sd
        , ltrait_mean
        , ltrait_sd
        , ifelse(mut_var == "beta", state$gamma, state$beta),sum(state$S, state$I)
          )
        ## DRY ...
        if (all(state$Imat==0) & host_dyn_only == FALSE) {
            message(sprintf("system went extinct prematurely (t=%d)",i))
            break
        }
        if (all(state$Svec==0) & host_dyn_only == TRUE) {
            message(sprintf("system went extinct prematurely (t=%d)",i,j))
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
    inf_rates <- state$beta*state$Imat*state$S*dt
    recover_rates <- state$gamma*state$Imat*dt
    return(c(inf_rates,recover_rates))
}
