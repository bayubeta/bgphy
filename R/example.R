library(devtools)
library(ape)
library(PCMBase)
library(Rcpp)
library(PCMBaseCpp)
library(PCMFit)
library(progress)
library(bayesplot)

load("data/modelOU.rda")
load("data/tree.rda")
load("data/Xtips.rda")

source("R/prior.R")
source("R/varchanges.R")
source("R/posterior.R")
source("R/utils.R")
Rcpp::sourceCpp("src/mcmc.cpp")

# define priors
priors <- prior(modelOU)
priors$X0 <- prior.normal(0)
priors$Theta <- prior.normal(1)
priors$H <- priors$Sigma_x <- prior.invgamma(2)

# transform priors
priors_tr <- lapply(priors, varchange)


# define log unnormalized posterior
lupost <- function(p, model, X, tree, priors, tr){

  # priors: unbounded priors

  # change parameters to p (bounded)
  p_b <- tr$g(p)
  prop_model <- setParams(p_b, model)

  # log likelihood
  loglik <- PCMBase::PCMLik(X, tree, prop_model, metaI = PCMInfoCpp)

  # log priors
  # p = (X0, H, Theta, Sigma_x)

  lp_X0 <- priors$X0(p[1])
  lp_H <- priors$H(p[2])
  lp_Theta <- priors$Theta(p[3])
  lp_Sigma_x <- priors$Sigma_x(p[4])

  sum_log <- loglik + lp_X0 + lp_H + lp_Theta + lp_Sigma_x

  if (is.infinite(sum_log) || is.na(sum_log)){
    return(Inf)
  }

  return(sum_log[1])
}



rwm <- function(model, X, tree, priors, initial, nsteps, scale, progress = TRUE, ncheck = 10){

  # number of parameters
  d <- length(initial)

  # a vector of rejection rates, one for each step
  R <- log(runif(nsteps))

  # a matrix of parameters, each row is a set of parameter for each step
  P <- matrix(NA, ncol = d, nrow = nsteps)
  colnames(P) <- c("X0", "H", "Theta", "Sigma_x")


  # define tr, functions to transform parameters
  # tr$f: to unbounded
  # tr$g: from unbounded
  tr <- trfunc(priors)


  # convert initial values into unbounded space
  initial <- tr$f(initial)

  # assign the first row as the initial parameter
  P[1,] <- initial


  # simplify lupost as a function of parameters
  lu_post <- function(p){lupost(p, model, X, tree, priors, tr)}


  # ==================== for progress bar ====================
  if (progress){
    total <- nsteps
    pb <- progress::progress_bar$new(format = "[:bar] :current/:total (:percent)", total = total)
    pb$tick(0)
    pb$tick(1)
    Sys.sleep(1/1000)
  }
  # ==========================================================


  # =========================== for checkpoints ==============
  if (ncheck){
    icheck <- round(nsteps / ncheck)
    }
  # ==========================================================



  # record acceptance rate
  ar <- c(1)


  # ================= start RWMH algorithm ==================

  # save values for current state
  pars0 <- initial
  logpost0 <- lu_post(pars0)


  for (i in 2:nsteps){

    # random walk move
    pars1 <- mcmcmove(pars0, scale)

    # calculate log posterior on the new set of parameters
    logpost1 <- lu_post(pars1)

    # accept
    if (logpost1 - logpost0 > R[i]){
      P[i,] <- pars1
      pars0 <- pars1
      logpost0 <- logpost1

      ar[i] <- 1
    }
    else{
      # or reject
      P[i,] <- pars0

      ar[i] <- 0
    }

    # progress bar
    if (progress){
      pb$tick(1)
      Sys.sleep(1 / 1000)
    }

    # ================= check point =================
    if (ncheck && (i %% icheck == 0)){
      save(P, file = sprintf("check%02d.rda", as.integer(i/icheck)))
    }
  }


  class(P) <- "posterior"


  # ================= transform samples back =================

  P <- transform(P, priors_tr)


  attr(P, "accept") <- mean(ar)

  return(P)

}




set.seed(1234)
f1 <- rwm(model = modelOU, X = Xtips, tree = tree, priors = priors_tr,
           initial = c(0,2,2,1), nsteps = 1000, scale = 0.2)

mcmc_hist(f1, pars = c("X0", "H", "Theta", "Sigma_x"))
mcmc_trace(f1, pars = c("X0", "H", "Theta", "Sigma_x"))
mcmc_areas(f1, pars = c("X0", "H", "Theta", "Sigma_x"),
           prob_outer = 0.95, point_est = "mean")






# define priors
priors2 <- prior(modelOU)
priors2$X0 <- prior.uniform(min = -3, 3)
priors2$Theta <- priors2$H <- priors2$Sigma_x <- prior.uniform(min = 0, 5)

# transform priors
priors_tr2 <- lapply(priors2, varchange)

set.seed(1234)
f2 <- rwm(model = modelOU, X = Xtips, tree = tree, priors = priors_tr2,
           initial = c(0,2,2,1), nsteps = 100000, scale = 0.2, ncheck = 10)

mcmc_hist(f2, pars = c("X0", "H", "Theta", "Sigma_x"))
mcmc_trace(f2, pars = c("X0", "H", "Theta", "Sigma_x"))
mcmc_areas(f2, pars = c("X0", "H", "Theta", "Sigma_x"),
           prob_outer = 0.95, point_est = "mean")


tr <- trfunc(priors_tr2)
lupost(f2[32103,], modelOU, Xtips, tree, priors_tr2, tr)


p <- f2[32103,]
p_b <- tr$g(p)
prop_model <- setParams(p_b, model)

# log likelihood
loglik <- PCMBase::PCMLik(X, tree, prop_model, metaI = PCMInfoCpp)

# log priors
# p = (X0, H, Theta, Sigma_x)

lp_X0 <- priors_tr2$X0(p[1])
lp_H <- priors_tr2$H(p[2])
lp_Theta <- priors_tr2$Theta(p[3])
lp_Sigma_x <- priors_tr2$Sigma_x(p[4])

sum_log <- loglik + lp_X0 + lp_H + lp_Theta + lp_Sigma_x







set.seed(1234)
fitOU <- PCMFit(model = PCM(model = "OU", k = 1), tree = tree, X = Xtips,
                metaI = PCMBaseCpp::PCMInfoCpp)

fitOU$modelOptim



















