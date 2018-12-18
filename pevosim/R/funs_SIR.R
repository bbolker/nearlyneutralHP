#######
## Expand the work in funs.R, which considers only parasite evolution (with host evolution as a true afterthought)
## to include host evolution. Requires a switch to an SIR model or at the very least an SIS model with host death
## as well as host reproduction
#######

library(arm); library(ggplot2); library(gridExtra); source("ggplot_theme.R")

#######
## Cautions / General notes
#######

## Don't choose gamma as the mutational parameter right now, it won't work.
## Don't give more than one host genotype as starting values, this hasn't been implemented yet.

## Decided to just remove continuous time for now to shorten script. It can be recovered from github if continuous time is desired
## Some debug notes:
  ## — See stats:::simulate.lm. get(".Random.seed", envir = .GlobalEnv)
  ## if (condition) browser()

#######
## Changes in progress and/or planned changes
#######

## 1 Host density dependent growth. For now assuming declining birth and constant death rate
##    1.1 Need to set up in the preamble of the function what b0 and decay rate are depending on pop size so that birth and death balances at N0
## 2 Parasites evolve in alpha and beta in a slightly funny way currently. Possibly update to:
  ## theta1 = log(alpha)
  ## theta2 = logit(beta / f(alpha)), where f is currently given by power_tradeoff()
##    2.1 Can also consider reparameterizing so that f(alpha) is not a function of alpha
##    2.2 Can also consider whether or not we want alpha and beta to evolve independently
## 3 Work on quasi-equilibrium within host model: exploitation rate vs resistance - what load does this end up being?
## 4 Convert to parasite exploitation and pathogenesis
## 5 Add host evolution with simple costs

######
## Accessory Functions
######

## Get host and parasite mutations
get_mut_h        <- function (state, orig_trait, mut_var, mut_mean, mut_sd, res_mut, tol_mut, mut_host_sd_shift, mut_host_mean_shift, mut_type = "shift") {

    if (mut_type=="shift") {

      ## Direction of change different for host and pathogen and different for the two parameters (beta and gamam)

      ## First get new resistant mutatnts
        new_trait_r <- orig_trait$res_mut +
            rnorm(length(orig_trait$res_mut)
            , ifelse(mut_var == "beta"
              , mut_mean*-1/mut_host_mean_shift
              , mut_mean/mut_host_mean_shift)
            , mut_sd/mut_host_sd_shift)  ## allow the sd for the host to be a multiple of the sd for the pathogen (1 for equal)

      ## Second, retrieve resistance values for the new tolerant mutatnts
        repeated_trait_r <- rep(state$hrtraitvec, tol_mut)
        new_trait_r      <- c(new_trait_r, repeated_trait_r)

         new_trait_t <- orig_trait$tol_mut +
            rnorm(length(orig_trait$tol_mut)
            , ifelse(mut_var == "beta"
              , mut_mean*-1/mut_host_mean_shift
              , mut_mean/mut_host_mean_shift)
            , mut_sd/mut_host_sd_shift)  ## allow the sd for the host to be a multiple of the sd for the pathogen (1 for equal)

        repeated_trait_t <- rep(state$httraitvec, res_mut)
        new_trait_t      <- c(repeated_trait_t, new_trait_t)

        if (any(is.na(c(new_trait_r, new_trait_t)))) stop("??")
        return(list(res_trait = new_trait_r, tol_trait = new_trait_t))
    } else stop("unknown mut_type")
}
get_mut_p        <- function (orig_trait, mut_var, mut_mean, mut_sd, mut_type = "shift") {

   if (mut_type=="shift") {
        new_trait_pos <- orig_trait$pos_trait +
            rnorm(length(orig_trait$pos_trait),
              ifelse(mut_var == "beta",  ## **beta and gamma in opposite directions
               mut_mean
            ,  mut_mean*-1)
            , mut_sd)

        new_trait_neg <- orig_trait$neg_trait +
            rnorm(length(orig_trait$neg_trait), mut_mean, mut_sd)

        if (any(is.na(c(new_trait_pos, new_trait_neg)))) stop("??")
        return(list(new_trait_pos = new_trait_pos, new_trait_neg = new_trait_neg))
    } else stop("unknown mut_type")
}
do_mut           <- function (state, mut_var, mut_host, orig_trait, ...) {
    ## If the mutation is in the parasite
    if (mut_host == FALSE) {
    new_trait        <- get_mut_p(orig_trait, mut_var, ...)
    state$ltraitvec  <- c(state$ltraitvec, new_trait$new_trait_pos)
    state$palphavec  <- c(state$palphavec, new_trait$new_trait_neg)
    } else {
    ## Else if the mutation is in the host (some adjustments to how state gets updated when the host is evolving)
    new_trait        <- get_mut_h(state, orig_trait, mut_var,...)
    state$hrtraitvec <- c(state$hrtraitvec, new_trait$res_trait)
    state$httraitvec <- c(state$httraitvec, new_trait$tol_trait)
    }
    return(state)
}
## Update mutant trait values using power-law tradeoff. Give further thought to scales of trait evolution and interaction
 ## Don't mutate gamma right now, it won't work
update_mut_pt    <- function (state, power_c, power_exp, mut_link_p, mut_link_h, mutated, mutated_host, mut_var, ...) {

   ## Scale beta according to tradeoff curve
   new_par_beta <- scale_beta_alpha(state, power_c, power_exp, mut_link_p, ...)

   ## Resistance will act to decrease parasite transmission and virulence following the shape of the tradeoff curve
     ## Need to think criticall about what scale this should be conducted on. Both logit and probability scale feel like
      ## they each have problems
   ## First calculate a tradeoff curve with the same curvature that passes through the parasite's (alpha, beta)
   cvec         <- pt_calc_c(beta = new_par_beta, alpha = mut_link_p$linkinv(state$palphavec), curv = power_exp)

   ## For each of these tradeoff curves, calculate a new alpha and beta for each host that is infected
   new_alphas   <- t(outer(mut_link_p$linkinv(state$palphavec), mut_link_h$linkinv(state$hrtraitvec), "*"))
   new_betas    <- matrix(power_tradeoff(cvec, alpha = c(new_alphas), curv = power_exp), nrow = nrow(new_alphas), ncol = ncol(new_alphas))

   state$alpha  <- new_alphas
   state$beta   <- new_betas

   ## Further adjust alpha via tolerance, which will act as a multiple to parasite intrinsic mortality rate
   state$alpha  <- sweep(state$alpha, 2, matrix(mut_link_h$linkinv(state$httraitvec), ncol = 1), "*")

   ## Update Infected matrix with the new strain, maintaining which S class received that mutation.
    ## For each mutated strain, Imat gets a new column with a single 1, in the row in which the mutation occurred
   if (sum(mutated) > 0) {
     ## First make a matrix of 0s, then add a single one in each column corresponding to the row of the host that the parasite mutated in
     new_mutes            <- matrix(data = 0, nrow = nrow(state$Imat), ncol = sum(mutated))
     num_mutes            <- matrix(c(rep(which(rowSums(mutated) > 0), rowSums(mutated)[rowSums(mutated) != 0])
       , 1:sum(mutated)), ncol = 2, nrow = sum(mutated))
     new_mutes[num_mutes] <- 1
     state$Imat           <- cbind(state$Imat, new_mutes)
   }

   if (sum(mutated_host) > 0) {
   ## Update Svec (mutated_host is a vector that tracks which host mutated)
   state$Svec       <- cbind(state$Svec, matrix(rep(1, sum(mutated_host)), nrow = 1))
   ## Add new rows with 0s for the new S genotypes
   state$Imat       <- rbind(state$Imat, matrix(data = rep(0, sum(mutated_host)*ncol(state$Imat)), ncol = ncol(state$Imat), nrow = sum(mutated_host)))
   }

   return(state)
}
## Mutation updating old version (no tradeoff -- no relationship between beta and alpha)
update_mut       <- function (state, mut_link, mutated, mutated_host, mut_var) {

   state[[mut_var]] <- t(mut_link$linkinv(outer(state$ltraitvec, state$hrtraitvec,  "/")))

   ## Parasite intrinsic mortality rate following some relationship yet to be determined
   state$palphavec  <- rep(state$palphavec[1], length(state$ltraitvec))

   ## Tolerance will act as a multiple to parasite intrinsic mortality rate
   state$alpha      <- t(outer(state$palphavec, state$hrtraitvec,  "*"))

   ## Update Infected matrix with the new strain, maintaining which S class received that mutation.
    ## For each mutated strain, Imat gets a new column with a single 1, in the row in which the mutation occurred
   if (sum(mutated) > 0) {
     ## First make a matrix of 0s, then add a single one in each column corresponding to the row of the host that the parasite mutated in
     new_mutes            <- matrix(data = 0, nrow = nrow(state$Imat), ncol = sum(mutated))
     num_mutes            <- matrix(data = c(rep(c(mutated), c(mutated)), 1:sum(mutated)), ncol = 2, nrow = sum(mutated))
     new_mutes[num_mutes] <- 1
     state$Imat           <- cbind(state$Imat, new_mutes)
   }

   if (sum(mutated_host) > 0) {
   ## Update Svec (mutated_host is a vector that tracks which host mutated)
   state$Svec       <- cbind(state$Svec, matrix(rep(1, sum(mutated_host)), nrow = 1))
   ## Add new rows with 0s for the new S genotypes
   state$Imat       <- rbind(state$Imat, matrix(data = rep(0, sum(mutated_host)*ncol(state$Imat)), ncol = ncol(state$Imat), nrow = sum(mutated_host)))
   }

   return(state)
}
## Remove extinct strains
do_extinct       <- function (state, mut_var, extinct, parasite) {
  ## If a parasite strain has gone extinct remove a column
  if (parasite == TRUE) {
    state[[mut_var]] <- state[[mut_var]][,-extinct, drop = FALSE]
  ## Again, will eventually need to clean this up
    state$alpha      <- state$alpha[,-extinct, drop = FALSE]
    state$ltraitvec  <- state$ltraitvec[-extinct, drop = FALSE]
    state$Imat       <- state$Imat[,-extinct, drop = FALSE]
 ## If a host strain has gone extinct remove a row
  } else {
    state[[mut_var]] <- state[[mut_var]][-extinct, , drop = FALSE]
    state$alpha      <- state$alpha[-extinct, , drop = FALSE]
    state$hrtraitvec <- state$hrtraitvec[-extinct, drop = FALSE]
    state$Imat       <- state$Imat[-extinct, , drop = FALSE]
    state$Svec       <- state$Svec[-extinct, drop = FALSE]
  }
    return(state)
}
## This function is borderline hideous but it does what I want it to do. All ears for a better way to write this
get_inf          <- function (Svec, uninf, Imat, beta) {

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
## Calculate death, either using a classic density dependence function, or a constant rate
get_death_con    <- function (state, d, S) {
  if (S == TRUE) {
  deathS     <- rbinom_mat(n = length(state$Svec), size = state$Svec, prob = d, nrow = nrow(state$Svec), ncol = ncol(state$Svec))
  state$Svec <- state$Svec - deathS
  } else {
  deathI     <- rbinom_mat(n = length(c(state$Imat)), size = c(state$Imat), prob = d, nrow = nrow(state$Imat), ncol = ncol(state$Imat))
  state$Imat <- state$Imat - deathI
  deathI     <- rbinom_mat(n = length(c(state$Imat)), size = c(state$Imat), prob = c(state$alpha), nrow = nrow(state$Imat), ncol = ncol(state$Imat))
  state$Imat <- state$Imat - deathI
  }
  state
}
get_death_dd     <- function (S, I, N0, d0, decay) {
  d0 * exp(-(S + I) / (N0 / decay))
}
## Calculate birth, either using a classic density dependence function, or a balancing function
get_birth_dd     <- function (state, N0, b0, decay) {
  birth_rate <- b0 * exp(-(sum(state$Svec) + sum(state$Imat)) / (N0 / decay))
  birth      <- rbinom_mat(n = length(state$Svec), size = state$Svec + rowSums(state$Imat) ## Assume S and I reproduce
    , prob = birth_rate, nrow = nrow(state$Svec), ncol = ncol(state$Svec))
  birth
}
get_birth_bal    <- function (state, d) {

 birth_rate <- (sum(state$Imat - (state$Imat * (1 - d) * (1 - d * state$alpha))) + (sum(state$Svec) * d)) / sum(state$Svec)
  if (birth_rate != 0 & sum(state$Svec > 0)) {
    birth <- rbinom_mat(n = length(state$Svec), size = state$Svec, prob = birth_rate, nrow = nrow(state$Svec), ncol = ncol(state$Svec))
  } else {
    birth <- matrix(rep(0, length(state$Svec)), nrow = 1)
  }
 birth

}
## Power law relationship between alpha and beta
power_tradeoff   <- function (alpha, c, curv) {
  c * alpha ^ (1 / curv)
}
pt_calc_c        <- function (alpha, beta, curv) {
 beta / ( alpha ^ (1 / curv) )
}
## Scale parasite beta and alpha evolution by the tradeoff
scale_beta_alpha <- function (state, power_c, power_exp, mut_link_p, ...) {

## Determine the beta of a given pathogen strain | that pathogen's current alpha value
 ## Maximum possible beta
max_beta      <- power_tradeoff(c = power_c, alpha = mut_link_p$linkinv(state$palphavec), curv = power_exp)

## realized beta
realized_beta <- mut_link_p$linkinv(state$ltraitvec) * max_beta

## Proportional change in beta relative to where it came from. (option to adjust alpha as correlated evolution instead of an evolving trait)
# beta_shift <-  realized beta / (plogis(orig_trait$pos_trait) *  max_beta)

## Could also assume either correlated evolution in one way or another.
 ## e.g an Alpha shift following beta evolution. Forces a positive-positive change or negative-negative change.
  ## For each beta shift choose an alpha shift that is in the same direction but somewhat random.
# alpha_shift <- apply(beta_shift, 2, function(x) ifelse(x > 1, runif(1, 1, x), runif(1, x, 1)))

return(realized_beta)

}
## Calculate starting trait values for parasite and host | desired starting R0
calc_startvals   <- function (alpha0, res0, tol0, gamma, d, R0_init, N, power_c, power_exp, mut_link_p) {

## alpha only considering host resistance
alpha_r    <- alpha0 * res0

## alpha given host resistance and tolerance
alpha_rt   <- alpha_r * tol0

## From these alphas back calculate joint beta,
joint_beta <- R0_init * (gamma + d + alpha_rt) / N

## then the c on which the suboptimal strain resides,
needed_c   <- pt_calc_c(alpha = alpha_r, beta = joint_beta, curv = power_exp)

## then beta without the effects of host genotype
beta_p     <- power_tradeoff(alpha = alpha0, c = needed_c, curv = power_exp)

## then parasite intrinsic beta (proportion of optimal beta | alpha)
intr_beta  <- mut_link_p$linkfun(beta_p / power_tradeoff(alpha = alpha0, c = power_c, curv = power_exp))

return(list(
  intrinsic_beta = intr_beta
, joint_beta     = joint_beta
, joint_alpha    = alpha_rt
  ))

}
## Wrappers for selecting random numbers from matrices and returning a matrix
rbinom_mat       <- function (n, size, prob, nrow, ncol) {
  matrix(rbinom(n, size, prob), nrow = nrow, ncol = ncol)
}

##' Logit transform on (min,max) rather than (0,1)
##' @param minval minimum value
##' @param maxval maximum value
##' @export
multlogit <- function(minval = 0, maxval = 1, scale = 1) {
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
#' @useDynLib pevosim
#' @importFrom Rcpp evalCpp
#' @param N population size
#' @param mu per-infection mutation probability
#' @param gamma recovery rate
#' @param lbeta0 initial beta
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
#' @param mut_link_p link function for mutation for parasite
#' @param mut_link_h link function for mutation for host
#' @param mut_host_sd_shift **Host multiple of parasite sd
#' @param mut_host_mean_shift **Host multiple of parasite mean (not set up yet)
#' @param mut_host_mu_shift **Host multiple of parasite mean (not set up yet). On the wrong scale (muliple on prob scale) but ok for now because of tiny prob....
#' @param res0 ## ^^starting host resistance value
#' @param tol0 ## ^^starting host tolerance value
#' @param Imat initial vector of infected numbers
#' @param mu mutation probability per replication
#' @param b host birth rate (for now assume only S births, but this may have to change) -- Currently birth rate isn't used, d is used to balance pop size
#' @param d host death rate for S
#' @param dt time step
#' @param seed random-number seed
#' @param nt number of time steps
#' @param rptfreq reporting frequency (should divide nt)
#' @param progress draw progress bar?
#' @param debug (logical) debugging output?
#' @param useCpp (logical) use C++ code for continuous-time models?
#' @importFrom stats make.link rnorm rbinom rmultinom rexp runif
#' @export

######
## Sim function (wrapper)
######

run_sim <- function(
   R0_init             = 2       ## >1
 , gamma0              = 1/5     ## >0
 , gamma_max           = Inf
 , alpha0              = 1       ## !! Critical piece is that this is the intrinsic parasite mortality probability (without influence of hosts)
 , d                   = 0.01
 , b                   = 0.1
 , b_decay             = 2.3
 , N                   = 1000    ## integer >0
 , mu                  = 0.01    ## >0
 , mut_type            = "shift"
 , mut_mean            = -1      ## <0 (for sensibility)
 , mut_sd              = 0.5     ## >0
 , mut_var             = "beta"
 , mut_link_p          = NULL
 , mut_link_h          = NULL
 , mut_host_sd_shift   = 1       ## **1 for identical sd to the parasite
 , mut_host_mean_shift = 1       ## **1 for identical mean to the parasite
 , mut_host_mu_shift   = 2
 , mut_host_res_bias   = 0.5     ## 0.5 means equal probability of evolving resistance or tolerance
 , res0                = 1       ## Starting host mean resistance value
 , res0_sd             = 0       ## Variation in resistance among starting host strains (if > 1 host strain). Not yet used, but plan to
 , tol0                = 1       ## Starting host mean tolerance value
 , tol0_sd             = 0       ## Variation in tolerance among starting host strains (if > 1 host strain). Not yet used, but plan to
 , power_c             = 0.75    ## Power law tradeoff scaling
 , power_exp           = 2       ## Power law tradeoff exponent
 , Imat                = NULL
 , host_dyn_only       = FALSE
 , dt                  = 1
 , nt                  = 100000
 , rptfreq             = max(nt / 500, 1)
 , seed                = NULL
 , progress            = FALSE
 , debug               = FALSE
 , useCpp              = FALSE) {

    if (round(N)!=N) {
        warning("rounding N")
        N <- round(N)
    }
    if (mut_mean>0) {
        warning("positive mutation bias")
    }
    stopifnot(
      R0_init > 0
    , gamma0  > 0
    , N       > 0
    , mu      > 0
    , mut_sd  > 0
    , (nt/rptfreq) %% 1 == 0)

    dfun <- function(lab="") {
        if (debug) {
            cat(lab,"\n")
            print(state)
        }
        return(NULL)
    }

    if (!is.null(seed)) set.seed(seed)

    if (is.null(Imat)) {
        ## start at equilibrium I ...
        Imat <- max(1, round(N*(1 - 1/R0_init)))
        Imat <- as.matrix(Imat)
    }

    Svec  <- N-Imat

   ## Set up mutation link
    if (mut_var=="beta") {
        if (is.null(mut_link_p)) mut_link_p <- make.link("logit")
        if (is.null(mut_link_h)) mut_link_h <- make.link("log")
      } else {
      ## Positive with log link can push gamma overboard when mutation in gamma is on average disadvantageous.
       ## Made to logit link for now, not sure...
        if (is.null(mut_link_p)) mut_link_p <- make.link("logit")
        if (is.null(mut_link_h)) mut_link_h <- make.link("log")
        }

    ## Determine intial alpha from parasite and host initial traits
     ## Note: Naming conventions are a bit odd here. alpha0 as a parameter in the function refers to the desired
      ## starting parasite mortality rate, irrespective of the host trait. That parameter is back transformed
       ## to the intrinsic parasite scale and used to calculate a true starting alpha (also called alpha0 here
        ## that is based on both parasite and host triats)

    startvals <- calc_startvals(alpha0, res0, tol0, gamma, d, R0_init, N, power_c, power_exp, mut_link_p)
    beta0     <- startvals$joint_beta
    ltraitvec <- startvals$intrinsic_beta

    ## Alpha (intrinsic parasite mortality pressure)
    palphavec  <- mut_link_p$linkfun(alpha0)

    ## Initial trait vectors for the host genotypes. Assumes all hosts start with identical traits (for now)
    hrtraitvec <- rep(mut_link_h$linkfun(res0), length(res0))
    httraitvec <- rep(mut_link_h$linkfun(tol0), length(tol0))

    ## Alpha0 now becomes starting host mortality rate as a function of both host and parasite gentoype
    alpha0     <-  startvals$joint_alpha

    ## Parameters structure (parallel vectors), so these can
     ## Be modified via function and passed back ...
    state <- list(
      beta       = as.matrix(beta0)
    , gamma      = as.matrix(gamma0) ## ^^as will gamma (make matrix after, because ltraitvec needs to remain a vector)
    , alpha      = as.matrix(alpha0)
    , ltraitvec  = ltraitvec
    , palphavec  = palphavec
    , hrtraitvec = hrtraitvec
    , httraitvec = httraitvec
    , Imat       = Imat
    , Svec       = Svec)

    ## For no infection (to check host dynamics in the absence of infection)
     ## Maybe not the _most_ iffecient place/way to do this, but I think makes the least clutter
    if (host_dyn_only == TRUE) {
      state$Imat[1, 1] <- 0
      state$Svec[1, 1] <- N
    }

    dfun("init")

    nrpt <- nt %/% rptfreq

     ## Added tracking host responses
    res <- as.data.frame(matrix(
      NA, nrow = nrpt, ncol = 11
    , dimnames = list(
      NULL
    , c("time"
      , "num_S"
      , "num_S_strains"
      , "num_I"
      , "num_I_strains"
      , "pop_size"
      , paste0(c("mean_hl", "sd_hl"), mut_var)
      , paste0(c("mean_pl", "sd_pl"), mut_var)
      , ifelse(mut_var == "beta", "gamma", "beta")  ##** list the non-evolving param
      ))))

    t_tot <- 0  ## for continuous model

    for (i in 1:nrpt) {
        ## cat("time ",i,"\n")

            for (j in 1:rptfreq) {
                ## cat("betavec:",betavec,"\n")
                ## cat("Imat:",Imat,"\n")

                ## [Step 1]: Birth. Not accessible to death or infection until the next time step.

                ## Ugly temp form of density dependence here for now keeping per-capita birth equal to the total death rate of both S and I
                 ## Switch to logistic growth possibly
                #birth      <- get_birth_bal(state, d = d0)
                ## dd birth currently set so natural birth and death balances at N
                 ## Later can be set up so that investment in resistance or tolerance changes b rate
                birth      <- get_birth_dd(state, N0 = N, b0 = b, decay = b_decay)

                ## [Step 2]: Death of S.
                state      <- get_death_con(state, d, S = TRUE)

                ## [Step 3]: Infection.
                ## Prob of escaping infection completely
                uninf      <- rbinom_mat(n = length(state$Svec), size = state$Svec, prob = prod((1-state$beta)^state$Imat), nrow = nrow(state$Svec), ncol = ncol(state$Svec))
                ## Division of new infectives among strains
                 ## 'prob' is internally normalized
                newinf     <- get_inf(state$Svec, uninf, state$Imat, beta = state$beta)
                stopifnot(sum(rowSums(newinf) + uninf) == sum(state$Svec))

                ## [Step 4]: Death of I.
                ## Natural death + parasite induced death
                state      <- get_death_con(state, d, S = FALSE)

                ## [Step 5]: Recovery of I. ##
                recover    <- rbinom_mat(n = length(c(state$Imat)), size = c(state$Imat), prob = c(state$gamma), nrow = nrow(state$Imat), ncol = ncol(state$Imat))

                ## [Step 6]: Mutation of new infections.
                ## Fraction of new infections -> mutation
               if (host_dyn_only == FALSE & sum(newinf) != 0) {
                  mutated  <- rbinom_mat(n = newinf, size = newinf, prob = mu, nrow = nrow(newinf), ncol = ncol(newinf))
               } else {
                  mutated  <- rbinom_mat(n = nrow(newinf), size = newinf, prob = mu, nrow = nrow(newinf), ncol = ncol(newinf))
               }

                ## [Step 7]: Mutation of new hosts (during birth).
                ## For now assume that S hosts are the only hosts reproducing -- a common assumption, unclear if it should stay.
                 ## Probably not when costs of resistance and tolerance are included
                mutated_host <- rbinom(length(birth), size = birth, prob = mu/mut_host_mu_shift)
                birth        <- birth - mutated_host

                ## Of the mutated hosts sort into resistance and tolerance mutants
                mutated_host_r <- rbinom(length(mutated_host), size = mutated_host, prob = mut_host_res_bias)
                mutated_host_t <- mutated_host - mutated_host_r

                stopifnot(length(recover)==length(mutated))
                stopifnot(length(newinf)==length(state$Imat))

                ## [Step 8]: Update Infecteds.
                state$Imat <- state$Imat - recover + newinf - mutated

                ## cat("I,uninf, unmutated, mutated, recovered",
                 ## c(sum(Imat),uninf,sum(newinf-mutated),sum(mutated),sum(recover)),
                  ## "\n")

                dfun("before mutation")
                if (debug) print(mutated)

                ## [Step 9]: Find new phenotypes for mutated parasites and hosts.
                tot_mut <- sum(mutated)

                ## Mut parasite first. Another choice of order that may/may not matter
                if (tot_mut>0) {
                    state <- do_mut(
                      state
                    , mut_var     = mut_var
                    , orig_trait  = list(
                        pos_trait = rep(state$ltraitvec, colSums(mutated))
                      , neg_trait = rep(state$palphavec, colSums(mutated)))  ## ^^Just care about intrinsic nature of a strain
                    , mut_mean    = mut_mean
                    , mut_sd      = mut_sd
                    , mut_host    = FALSE
                    , mut_type    = mut_type)
                }

                ## Host mutations
                if (sum(mutated_host) > 0) {
                    state <- do_mut(
                      state
                    , mut_var             = mut_var
                    , orig_trait          = list(
                        res_mut           = rep(state$hrtraitvec, mutated_host_r)
                      , tol_mut           = rep(state$httraitvec, mutated_host_t))
                    , res_mut             = mutated_host_r
                    , tol_mut             = mutated_host_t
                    , mut_mean            = mut_mean
                    , mut_sd              = mut_sd
                    , mut_host            = TRUE
                    , mut_host_sd_shift   = mut_host_sd_shift
                    , mut_host_mean_shift = mut_host_mean_shift
                    , mut_type            = mut_type)
                }

                ## [Step 10]: Update Svec with infections and recoveries prior to the mutations (and host birth).
                state$Svec <- state$Svec + rowSums(recover) - rowSums(newinf) + birth

                ## [Step 11]: Update state (a big component beign parasite alpha and beta) with new parasite and host evolution.
                if (tot_mut > 0 | sum(mutated_host) > 0) {
                state <- update_mut_pt(
                    state        = state
                  , power_c      = power_c
                  , power_exp    = power_exp
                  , mut_link_p   = mut_link_p
                  , mut_link_h   = mut_link_h
                  , mutated      = mutated
                  , mutated_host = mutated_host
                  , mut_var      = mut_var)
                }

                if (sum(state$Imat)==0 & host_dyn_only == FALSE) {
                    message(sprintf("system went extinct prematurely (t=%d)",i))
                    break
                }
                if (sum(state$Svec)==0 & host_dyn_only == TRUE) {
                    message(sprintf("system went extinct prematurely (t=%d)",i,j))
                    break
                }

                ## Look over columns of I for extinct parasites and over rows of I and entries of S for extinct hosts
                extinct_p <- which(colSums(state$Imat) == 0)
                extinct_h <- which(rowSums(state$Imat) == 0 & state$Svec == 0)
                if (length(extinct_p) > 0 & host_dyn_only == FALSE) {
                    state <- do_extinct(state, mut_var ,extinct = extinct_p, parasite = TRUE)
                }
                if (length(extinct_h) > 0) {
                    state <- do_extinct(state, mut_var, extinct = extinct_h, parasite = FALSE)
                }
                dfun("after mutation")

            }  ## rptfreq time steps

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
        , ifelse(mut_var == "beta", state$gamma[1,1], state$gamma[1,1]))

        ## DRY ...
        if (sum(state$Imat) == 0 & host_dyn_only == FALSE) {
            message(sprintf("system went extinct prematurely (t=%d)",i))
            break
        }
        if (sum(state$Svec) == 0 & host_dyn_only == TRUE) {
            message(sprintf("system went extinct prematurely (t=%d)",i,j))
            break
        }

    } ## loop over reporting frequencies
    if (progress) cat("\n")
    return(res)
}

#' Get rates
#' @param state list of state vectors
#' @export
get_rates <- function(state, dt = 1) {
    inf_rates <- state$beta*state$Imat*state$S*dt
    recover_rates <- state$gamma*state$Imat*dt
    return(c(inf_rates,recover_rates))
}
