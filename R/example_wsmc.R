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
library(mvtnorm)

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

npars <- length(priors)


# sequence of exponents
Ttime <- 50
y <- seq(0, 1, length.out = Ttime)


MCMC_kernel <- function(p, yt, lp, em_var){
  U <- log(runif(1))
  npars <- length(p)

  p_new <- mvtnorm::rmvnorm(1, mean = p, sigma = diag(em_var))

  lp_prop <- lp(p_new)
  lp_now <- lp(p)

  alpha <- yt*(lp_prop - lp_now)

  if (alpha > U){return(c(p_new, lp_prop))}
  else{return(c(p, lp_now))}
}

# p = p_sampler(1)
# llik_k = function(p){lupost(p, modelOU, Xsim, tree, priors_tr, tr)}
# MCMC_kernel(p, llik_k)


model <- modelOU
X <- Xsim

# define log posterior as a function of p
lp <- function(p){lupost(p, model, X, tree, priors_tr, tr)}


# number of particles
N <- 1000


timestamp()
set.seed(1234)
# initial values (step 0)
P <- t(apply(p_sampler(N), 1, tr$f)) # transform to unbounded
W <- rep(1/N, N)
logW <- log(W)



##########################




cl <- parallel::makeCluster(parallel::detectCores())

PCMInfoCpp <- PCMBaseCpp::PCMInfoCpp
PCMParamLoadOrStore <- PCMBase::PCMParamLoadOrStore
parallel::clusterExport(cl,c("P", "lp", "PCMInfoCpp", "PCMParamLoadOrStore"))

envir_lp <- environment(lp)
parallel::clusterExport(cl, varlist = ls(envir_lp), envir = envir_lp)


# logpost values, at P_k-1
lp_t_1 <- parallel::parApply(cl, P, 1, lp)


for (t in 2:Ttime){
  start <- Sys.time()
  # reweight
  logW_tilde <- logW + (lp_t_1)*(y[t] - y[t-1])
  # and normalize
  logW <- logW_tilde - logsumexp(logW_tilde)


  # check ESS
  ESS <- 1/exp(logsumexp(2*logW))
  if (ESS <= (N/2)){
    # if ESS is low, resample
    # particle index
    A <- sample(1:N, N, replace = TRUE, prob = exp(logW))
    P <- P[A,]
    logW <- log(rep(1/N, N))
  }

  # empirical variance of the particles
  em_var <- 2*(diag(cov(P)))

  # sample new values of P, along with its logpost
  P_l <- t(parallel::parApply(cl, P, 1, MCMC_kernel, yt = y[t], lp = lp, em_var = em_var))


  # new samples
  P <- P_l[,1:npars]
  # lp for current particle
  lp_t <- P_l[,npars]


  (Sys.time() - start)[1]

  cat("t = ", t,"; ",
      as.numeric(difftime(Sys.time(), start, units = "secs"))," s\n")

  # assign current logpost values into old values
  lp_t_1 <- lp_t
}

parallel::stopCluster(cl)


# transform back
P_post <- t(apply(P, 1, tr$g))
colnames(P_post) <-  names(priors)

timestamp()
mcmc_dens(P_post)





