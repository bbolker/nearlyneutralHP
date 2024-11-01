#######
## Expand the work in funs.R, which considers only parasite evolution (with host evolution as a true afterthought)
## to include host evolution. Requires a switch to an SIR model or at the very least an SIS model with host death
## as well as host reproduction
#######

#######
## New Cautions and notes upon returning on Feb 1
#######

### First, search *** For new changes based on updated logic

## (1) Unclear that the two step evolutionary process of first evo in alpha (while holding beta proportional to its
 ## max beta according to the tradeoff curve), followed by evo in beta is an appropriate method for modeling the evolution
## (2) Unclear that treating birth as a probability is appropraite. This mostly arises at low pop sizes, where births are
 ## set to balance deaths. As virulence increases deaths get high, which forces the probability to 1, removing a piece of the
  ## stochasticity. Maybe it is ok, but I definitely need to think about it. It may be more appropriate to actually model the
   ## number of newborns per adult...
## (3)


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

    if (mut_type == "shift") {

      ## Direction of change different for host and pathogen and different for the two parameters (beta and gamam)

      ## First get new resistant mutatnts
        new_trait_r <- orig_trait$res_mut +
            rnorm(length(orig_trait$res_mut)
            , ifelse(mut_var == "beta"
              , mut_mean/mut_host_mean_shift
              , mut_mean*-1/mut_host_mean_shift)
            , mut_sd/mut_host_sd_shift)  ## allow the sd for the host to be a multiple of the sd for the pathogen (1 for equal)

      ## Second, retrieve resistance values for the new tolerant mutatnts
        repeated_trait_r <- rep(state$hrtraitvec, tol_mut)
        new_trait_r      <- c(new_trait_r, repeated_trait_r)

         new_trait_t <- orig_trait$tol_mut +
            rnorm(length(orig_trait$tol_mut)
            , ifelse(mut_var == "beta"
#              , mut_mean/mut_host_mean_shift
              , 0
              , mut_mean*-1/mut_host_mean_shift)
            , mut_sd/mut_host_sd_shift)  ## allow the sd for the host to be a multiple of the sd for the pathogen (1 for equal)

        repeated_trait_t <- rep(state$httraitvec, res_mut)
        new_trait_t      <- c(repeated_trait_t, new_trait_t)

        if (any(is.na(c(new_trait_r, new_trait_t)))) stop("??")
        return(list(res_trait = new_trait_r, tol_trait = new_trait_t))
    } else stop("unknown mut_type")
}
get_mut_p        <- function (orig_trait, mut_var, power_c, power_exp, mut_link_p, mut_mean, mut_sd, mut_type = "shift", agg_eff_adjust) {

   if (mut_type == "shift") {

     ## Parasite first evolves in alpha
    #   new_trait_neg <- orig_trait$neg_trait + rnorm(length(orig_trait$neg_trait), mut_mean*-1, mut_sd) ## Parasites evolve to increase mortality rate

      ## Assume a neutrally evolving aggressiveness
       neg_trait_adj <- rnorm(length(orig_trait$neg_trait), 0, mut_sd)
       new_trait_neg <- orig_trait$neg_trait + neg_trait_adj

     ## but evolving in alpha necessarily changes what parasite intrinsic beta means, because parasite intrinsic beta
      ## is how close a parasite is to maximizing beta at the current alpha.

     ## To make alpha evolve as an independent trait, the beta of the new strain needs to be updated so that the parasite
      ## is the same proportional distance to its optimum on the tradeoff curve

     ## Find optimum possible beta with the new alpha, and determine how close original beta would be to the new optima
#        orig_beta <- power_tradeoff(alpha = mut_link_p$linkinv(orig_trait$neg_trait)
#          , c = power_c, curv = power_exp) * mut_link_p$linkinv(orig_trait$pos_trait)

      ## Propagate a new beta to accompany the new alpha from which beta can independently evolve
#        orig_trait$pos_trait <- logit(power_tradeoff(alpha = mut_link_p$linkinv(new_trait_neg), c = power_c, curv = power_exp) *
#            mut_link_p$linkinv(orig_trait$pos_trait))

### *** (Return to code note, Feb 1) In this step whole mutation step, because of how it is set up, the intrinsic parasite beta
 ### *** that is evolving is just how close the parasite is to its maximum beta constrained by the tradeoff curve. Thus it is able
  ### *** to evolve independently. In the later step the intrinsic beta is translated to a new realized beta given the new location
   ### *** on the tradeoff curve that it finds itslf on given by alpha

   ## if an increase in parasite aggressiveness leads to a decrease in parasite efficiency, first adjust
    ## efficiency according to how much aggressiveness changed
       if (agg_eff_adjust == TRUE) {

  #### ******
    ## Does the order of these two things matter? Is something about getting this first adjustment impacting how
      ## the biased mutation is realized in the second step?

      new_trait_pos <- orig_trait$pos_trait - neg_trait_adj

          ## Assume a negatively evolving efficiency
        new_trait_pos <- new_trait_pos +
            rnorm(length(new_trait_pos),
              ifelse(mut_var == "beta",  ## Beta and gamma in opposite directions
               mut_mean
            ,  mut_mean*-1)
            , mut_sd
     #       , 0
              )

       } else {

         new_trait_pos <- orig_trait$pos_trait +
            rnorm(length(orig_trait$pos_trait),
              ifelse(mut_var == "beta",  ## Beta and gamma in opposite directions
               mut_mean
            ,  mut_mean*-1)
            , mut_sd)

       }

### *** Can also set up to not allow any update in
#new_trait_pos <- orig_trait$pos_trait

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
update_mut_pt    <- function (state, orig_trait, power_c, power_exp, mut_link_p, mut_link_h, mutated, mutated_host, mut_var, ...) {

   ## Scale beta according to tradeoff curve
   new_par_beta <- scale_beta_alpha(state, power_c, power_exp, mut_link_p, ...)

   ## Resistance will act to decrease parasite transmission and virulence following the shape of the tradeoff curve
     ## Need to think criticall about what scale this should be conducted on. Both logit and probability scale feel like
      ## they each have problems
   ## First calculate a tradeoff curve with the same curvature that passes through the parasite's (alpha, beta)
   cvec         <- pt_calc_c(beta = new_par_beta, alpha = mut_link_p$linkinv(state$palphavec), curv = power_exp)

   ## For each of these tradeoff curves, calculate a new alpha and beta for each host that is infected
  # new_alphas  <- t(outer(mut_link_p$linkinv(state$palphavec), mut_link_h$linkinv(state$hrtraitvec), "*"))
   new_alphas   <- t(mut_link_p$linkinv(outer(state$palphavec, state$hrtraitvec, "-")))
   new_betas    <- matrix(power_tradeoff(rep(cvec, each = nrow(new_alphas)), alpha = c(new_alphas)
     , curv = power_exp), nrow = nrow(new_alphas), ncol = ncol(new_alphas))

   state$alpha  <- new_alphas
   state$beta   <- new_betas

   #######
   ## ** ^^ Important order problem??? Because of the nonlinear scaling, applying resistance first and then tolerance second,
    ## tolerance will have a smaller effect. But tolerance doesn't affect beta, so unclear how to fix this.
   #######

   ## Further adjust alpha via tolerance, which will act as a multiple to parasite intrinsic mortality rate
   state$alpha  <- get_alpha_tol(state, numcol = ncol(state$alpha), numrow = nrow(state$alpha), mut_link_h, mut_link_p)

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
   state$Imat       <- rbind(state$Imat, matrix(data = rep(0, sum(mutated_host)*ncol(state$Imat))
     , ncol = ncol(state$Imat), nrow = sum(mutated_host)))
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
    state$palphavec  <- state$palphavec[-extinct, drop = FALSE]
    state$Imat       <- state$Imat[,-extinct, drop = FALSE]
 ## If a host strain has gone extinct remove a row
  } else {
    state[[mut_var]] <- state[[mut_var]][-extinct, , drop = FALSE]
    state$alpha      <- state$alpha[-extinct, , drop = FALSE]
    state$hrtraitvec <- state$hrtraitvec[-extinct, drop = FALSE]
    state$httraitvec <- state$httraitvec[-extinct, drop = FALSE]
    state$Imat       <- state$Imat[-extinct, , drop = FALSE]
    state$Svec       <- state$Svec[, -extinct, drop = FALSE]
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
 , prob = beta * Imat
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
  return(list(state, deathS))
  } else {
  deathI1     <- rbinom_mat(n = length(c(state$Imat)), size = c(state$Imat), prob = d, nrow = nrow(state$Imat), ncol = ncol(state$Imat))
  state$Imat  <- state$Imat - deathI1
  deathI2     <- rbinom_mat(n = length(c(state$Imat)), size = c(state$Imat), prob = c(state$alpha), nrow = nrow(state$Imat), ncol = ncol(state$Imat))
  state$Imat  <- state$Imat - deathI2
  return(list(state, deathI1, deathI2))
  }
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

  ## Set birth rate equal to the total death rate among all infected individuals and susceptible individuals from parasite virulence and background death
 birth_rate <- (sum(state$Imat - (state$Imat * (1 - d) * (1 - state$alpha))) + (sum(state$Svec) * d)) / sum(state$Svec)
  if (birth_rate != 0 & sum(state$Svec > 0)) {
    if (birth_rate < 1) {
    birth <- rbinom_mat(n = length(state$Svec), size = state$Svec, prob = birth_rate, nrow = nrow(state$Svec), ncol = ncol(state$Svec))
    } else {
  ## If the number of deaths is bigger than the number of S, keep birth marginally stochastic by adding births equal to the difference
   ## between death and birth, but at least double the number of S by allowing each S to reproduce
    birth <- state$Svec * floor(birth_rate) + rbinom_mat(n = length(state$Svec), size = state$Svec, prob = (birth_rate - floor(birth_rate)), nrow = nrow(state$Svec), ncol = ncol(state$Svec))
    }
  } else {
    birth <- matrix(rep(0, length(state$Svec)), nrow = 1)
  }
 birth

}
get_birth_det    <- function (deathS, deathI) {
  ## Total birth per S class is the sum of the number of deaths in S and death in the I classes
  rowSums(deathI[[3]]) + rowSums(deathI[[2]]) + deathS[[2]]
  }
## Function to wrap apply into returning a matrix (problems with always needing a matrix and reducing to a vector sometimes)
get_alpha_tol    <- function (state, numcol, numrow, mut_link_h, mut_link_p) {

  mut_link_p$linkinv(
  matrix(
      apply(mut_link_p$linkfun(state$alpha), 2, function (x) x - matrix(state$httraitvec, ncol = 1))
    , nrow = numrow
    , ncol = numcol
  )
  )

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
calc_startvals   <- function (alpha0, res0, tol0, gamma0, d, R0_init, N, power_c, power_exp, mut_link_h, mut_link_p) {

## alpha only considering host resistance
alpha_r    <- mut_link_p$linkinv(mut_link_p$linkfun(alpha0) - mut_link_h$linkfun(res0))

## alpha given host resistance and tolerance
alpha_rt   <- mut_link_p$linkinv(mut_link_p$linkfun(alpha_r) - mut_link_h$linkfun(tol0))

## From these alphas back calculate joint beta,
joint_beta <- R0_init * (gamma0 + d + alpha_rt) / N

## then the c on which the suboptimal strain resides,
needed_c   <- pt_calc_c(alpha = alpha_r, beta = joint_beta, curv = power_exp)

## then beta without the effects of host genotype
beta_p     <- power_tradeoff(alpha = alpha0, c = needed_c, curv = power_exp)

## then parasite intrinsic beta (proportion of optimal beta | alpha). Intrinsic beta is efficiency, however
 ## I decide to calculate it. If parasite tuning is introduced, then efficiency is defined by the matching
  ## of parasite
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
multlogit        <- function(minval = 0, maxval = 1, scale = 1) {
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
#' @param mut_host_sd_shift Host multiple of parasite sd
#' @param mut_host_mean_shift Host multiple of parasite mean (not set up yet)
#' @param mut_host_mu_shift Host multiple of parasite mean (not set up yet). On the wrong scale (muliple on prob scale) but ok for now because of tiny prob....
#' @param res0 ## starting host resistance value
#' @param tol0 ## starting host tolerance value
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
 , balance_birth       = TRUE
 , stochastic_birth    = TRUE
 , dt                  = 1
 , nt                  = 100000
 , rptfreq             = max(nt / 500, 1)
 , seed                = NULL
 , progress            = FALSE
 , debug               = FALSE
 , debug2              = FALSE
 , debug3              = FALSE
 , host_evo_delay      = FALSE
 , host_evo_delay_start= NULL
 , host_evo_delay_stop = NULL
 , agg_eff_adjust      = FALSE
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
    startvals  <- calc_startvals(alpha0, res0, tol0, gamma0, d, R0_init, N, power_c, power_exp, mut_link_h, mut_link_p)
    beta0      <- startvals$joint_beta
    ltraitvec  <- startvals$intrinsic_beta

    ## Alpha (intrinsic parasite mortality pressure)
    palphavec  <- mut_link_p$linkfun(alpha0)

    ## Initial trait vectors for the host genotypes. Assumes all hosts start with identical traits (for now)
    hrtraitvec <- rep(mut_link_h$linkfun(res0), length(res0))
    httraitvec <- rep(mut_link_h$linkfun(tol0), length(tol0))

    ## Alpha0 now becomes starting host mortality rate as a function of both host and parasite gentoype
    alpha0     <-  startvals$joint_alpha

    ## Parameters structure (parallel vectors), so these can
     ## Be modified via function and passed back ...
    state      <- list(
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
      NA, nrow = nrpt, ncol = 19
    , dimnames = list(
      NULL
    , c("time"
      , "num_S"
      , "num_S_strains"
      , "num_I"
      , "num_I_strains"
      , "pop_size"
      , paste0(c("mean_hl", "sd_hl"), "res")
      , paste0(c("mean_hl", "sd_hl"), "tol")
      , paste0(c("mean_pl", "sd_pl"), mut_var)
      , paste0(c("mean_pl", "sd_pl"), "alpha")
      , "mean_alpha"
      , "sd_alpha"
      , "mean_beta"
      , "sd_beta"
      , ifelse(mut_var == "beta", "gamma", "beta")  ##** list the non-evolving param
      ))))

    t_tot <- 0  ## for continuous model

    for (i in 1:nrpt) {
        ## cat("time ",i,"\n")

      if (host_evo_delay == TRUE & i == host_evo_delay_start) {
      mut_host_mu_shift <- 10
      }

      if (host_evo_delay == TRUE & i >= host_evo_delay_stop) {
      state$httraitvec <- state$httraitvec / 1.05
      mut_host_mu_shift <- 1000000000
      }

            for (j in 1:rptfreq) {
                ## cat("betavec:",betavec,"\n")
                ## cat("Imat:",Imat,"\n")

               ## Cut and paste to wherever there is a problem
                if (debug2 == TRUE) {
                  print(paste(i, j, sep = "  -  "))
                  print(str(state$Svec))
                  assign("state_check", state, .GlobalEnv)
                  if (i == 1 & j == 5) browser()  ## Fill these in manually with printed i and j
                }

                ## [Step 1]: Birth. Not accessible to death or infection until the next time step.

         #      if (i == 200 & j == 1) browser()

                ## Two options currently for birth. Balancing (birth = all death) and density dependent (decreasing with N)
                 ## S and I both reproduce here currrently: can change later so that investment in res and tol changes b (or I vs S)
              if (balance_birth == FALSE) {
                birth      <- get_birth_dd(state, N0 = N, b0 = b, decay = b_decay)
              } else {
                if (stochastic_birth == TRUE) {
                birth      <- get_birth_bal(state, d)
                }
              }

                ## [Step 2]: Death of S. Returning a list of updated state and death, so death can be used to calculate deterministic birth
                deathS     <- get_death_con(state, d, S = TRUE)
                state      <- deathS[[1]]

                ## [Step 3]: Infection.
                ## Prob of escaping infection completely
                uninf      <- rbinom_mat(n = length(state$Svec), size = state$Svec, prob = prod((1-state$beta)^state$Imat), nrow = nrow(state$Svec), ncol = ncol(state$Svec))
                ## Division of new infectives among strains
                 ## 'prob' is internally normalized
                newinf     <- get_inf(state$Svec, uninf, state$Imat, beta = state$beta)
                stopifnot(sum(rowSums(newinf) + uninf) == sum(state$Svec))

                ## [Step 4]: Death of I.
                ## Natural death + parasite induced death
                deathI     <- get_death_con(state, d, S = FALSE)
                state      <- deathI[[1]]

                ## Deterministic birth
                if (balance_birth == TRUE & stochastic_birth == FALSE) {
                birth      <- get_birth_det(deathS, deathI)
                }

                ## [Step 5]: Recovery of I. ##
                recover    <- rbinom_mat(n = length(c(state$Imat)), size = c(state$Imat), prob = c(state$gamma), nrow = nrow(state$Imat), ncol = ncol(state$Imat))

                ## [Step 6]: Mutation of new infections.
                ## Fraction of new infections -> mutation
               if (host_dyn_only == FALSE & sum(newinf) != 0) {
                  mutated  <- rbinom_mat(n = newinf, size = newinf, prob = mu, nrow = nrow(newinf), ncol = ncol(newinf))
               } else {
                  mutated  <- rbinom_mat(n = nrow(newinf), size = newinf, prob = mu, nrow = nrow(newinf), ncol = ncol(newinf))
               }

                if (debug3 == TRUE) { mutated[1,1] <- 2 }

                ## [Step 7]: Mutation of new hosts (during birth).
                ## For now assume that S hosts are the only hosts reproducing -- a common assumption, unclear if it should stay.
                 ## Probably not when costs of resistance and tolerance are included
                mutated_host <- rbinom(length(birth), size = birth, prob = mu/mut_host_mu_shift)
                ## Of the mutated hosts sort into resistance and tolerance mutants
                mutated_host_r <- rbinom(length(mutated_host), size = mutated_host, prob = mut_host_res_bias)
                mutated_host_t <- mutated_host - mutated_host_r

                if (debug3 == TRUE) { mutated_host[1] <- 2; mutated_host_r[1] <- 1; mutated_host_t[1] <- 1 }

                stopifnot(length(recover) == length(mutated))
                stopifnot(length(newinf)  == length(state$Imat))

                ## Updated birth from
                birth      <- birth - mutated_host

                ## [Step 8]: Update Infecteds.
                state$Imat <- state$Imat - recover + newinf - mutated

                ## cat("I,uninf, unmutated, mutated, recovered",
                 ## c(sum(Imat),uninf,sum(newinf-mutated),sum(mutated),sum(recover)),
                  ## "\n")

                dfun("before mutation")
                if (debug) print(mutated)

                ## [Step 9]: Find new phenotypes for mutated parasites and hosts.
                tot_mut <- sum(mutated)

                ## Original parasite taits that mutants arise from. Need this for later
                orig_trait <- list(
                   pos_trait = rep(state$ltraitvec, colSums(mutated))
                 , neg_trait = rep(state$palphavec, colSums(mutated)))

                ## Mut parasite first. Another choice of order that may/may not matter
                if (tot_mut > 0) {
                    state <- do_mut(
                      state
                    , mut_var        = mut_var
                    , orig_trait     = orig_trait  ## ^^Just care about intrinsic nature of a strain
                    , mut_mean       = mut_mean
                    , mut_sd         = mut_sd
                    , mut_host       = FALSE
                    , mut_type       = mut_type
                    , power_c        = power_c
                    , power_exp      = power_exp
                    , mut_link_p     = mut_link_p
                    , agg_eff_adjust = agg_eff_adjust)
                }

                ## Host mutations
                if (length(which(is.na(mutated_host)) > 0)) {
                browser()
                }

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
                  , orig_trait   = orig_trait
                  , power_c      = power_c
                  , power_exp    = power_exp
                  , mut_link_p   = mut_link_p
                  , mut_link_h   = mut_link_h
                  , mutated      = mutated
                  , mutated_host = mutated_host
                  , mut_var      = mut_var)
                }

                if (sum(state$Imat)==0 & host_dyn_only == FALSE) {
                    message(sprintf("system went extinct prematurely (t=%d)", i))
                    break
                }
                if (sum(state$Svec)==0 & host_dyn_only == TRUE) {
                    message(sprintf("system went extinct prematurely (t=%d)", i))
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
        lalpha_mean  <- sum(colSums(state$Imat)*state$palphavec)/num_I
        lalpha_sd    <- sqrt(sum(colSums(state$Imat)*(state$palphavec-ltrait_mean)^2)/num_I)
        num_S        <- sum(state$Svec)
        S_tot        <- length(state$Svec)
        lhres_mean   <- sum(state$Svec*state$hrtraitvec)/num_S
        lhres_sd     <- sqrt(sum(state$Svec*(state$hrtraitvec-lhres_mean)^2)/num_S)
        lhtol_mean   <- sum(state$Svec*state$httraitvec)/num_S
        lhtol_sd     <- sqrt(sum(state$Svec*(state$httraitvec-lhtol_mean)^2)/num_S)
        pop_size     <- num_I + num_S

        ## actual alpha and beta of all parasites in all host classes
        avg_alpha    <- mean(state$alpha)
        sd_alpha     <- mean(state$alpha)
        avg_beta     <- mean(state$beta)
        sd_beta      <- mean(state$beta)

        if (progress) cat(".")
        res[i,] <- c(
          i*rptfreq
        , num_S
        , S_tot
        , num_I
        , I_tot
        , pop_size
        , lhres_mean
        , lhres_sd
        , lhtol_mean
        , lhtol_sd
        , ltrait_mean
        , ltrait_sd
        , lalpha_mean
        , lalpha_sd
        , avg_alpha
        , sd_alpha
        , avg_beta
        , sd_beta
        , ifelse(mut_var == "beta", state$gamma[1,1], state$gamma[1,1]))

        ## DRY ...
        if (sum(state$Imat) == 0 & host_dyn_only == FALSE) {
            message(sprintf("system went extinct prematurely (t=%d)", i))
            break
        }
        if (sum(state$Svec) == 0 & host_dyn_only == TRUE) {
            message(sprintf("system went extinct prematurely (t=%d)", i))
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
