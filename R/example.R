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
# define priors
priors <- prior(modelOU)
priors$X0 <- priors$H <- priors$Theta <- prior.normal(mean = 0, sd = 5)
priors$Sigma_x <- prior.halfcauchy(4)

priors

timestamp()
set.seed(1234)
f2 <- rwm(model = modelOU, X = Xtips, tree = tree, priors = priors,
          initial = c(1,1,1,1), nsteps = 1000, scale = 0.1, ncheck = FALSE)
timestamp()

mcmc_hist(f2)
mcmc_trace(f2)
