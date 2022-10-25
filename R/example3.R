library(devtools)
library(ape)
library(PCMBase)
library(Rcpp)
library(PCMBaseCpp)
library(PCMFit)
library(progress)
library(bayesplot)
library(TreeSim)
library(ggtree)

load("data/modelOU.rda")
load("data/tree.rda")
load("data/Xsim.rda")

source("R/prior.R")
source("R/varchanges.R")
source("R/posterior.R")
source("R/utils.R")
source("R/rwm.R")
Rcpp::sourceCpp("src/mcmc.cpp")



# define priors
# define priors
priors <- prior(modelOU)
priors$X0 <- priors$H <- priors$Theta <- prior.normal(mean = 0, sd = 5)
priors$Sigma_x <- prior.halfnormal(3)

# priors sampler
p_sampler <- prior_sampler(priors)
p_sampler

# transformed priors
priors_tr <- lapply(priors, varchange)

# define tr, functions to transform parameters
# tr$f: to unbounded
# tr$g: from unbounded
tr <- trfunc(priors_tr)


# sequence of exponents
K <- 50
y <- seq(0, 1, length.out = K)


MCMC_kernel <- function(p, llik_k){
  U = log(runif(1))
  npars <- length(p)
  p_new <- p + rnorm(npars, 0, sd = 0.01)

  alpha <- min(c(0, llik_k(p_new)-llik_k(p)))
  if (U < alpha){
    return(p_new)
  }
  else{
    return(p)
  }
}

# p = p_sampler(1)
# llik_k = function(p){lupost(p, modelOU, Xsim, tree, priors_tr, tr)}
# MCMC_kernel(p, llik_k)

# define tempered log posterior at time k
lp_k <- function(p,y){lupost_k(p, y, modelOU, Xsim, tree, priors_tr, tr)}


# number of particles
N <- 100

# initial values (step 0)
P <- t(apply(p_sampler(N), 1, tr$f)) # transform to unbounded
W <- rep(1/N, N)
logW <- log(W)

# logpost_k-1 values, at P_k-1
lp_1 <- apply(X, 1, lp_k, y = y[1])


for (k in 2:K){
  # logpost_k, at P_k-1
  lp <- apply(X, 1, lp_k, y = y[k])

  # reweight
  logW_tilde <- logW + lp - lp_1
  # and normalize
  logW <- logW_tilde - logsumexp(logW_tilde)

  # check ESS
  ESS <- 1/exp(logsumexp(2*logW))
  if (ESS <= 50){
    # if ESS is low, resample
    # particle index
    A <- sample(1:N, N, replace = TRUE, prob = exp(logW))
    P <- P[A,]
    logW <- log(rep(1/N, N))
  }

  # sample new values of P
  P <- t(apply(P, 1, MCMC_kernel, llik_k = function(p){lp_k(p, y[k])}))
  cat("k = ", k, "\n")
}






















