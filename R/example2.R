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

load("data/modelMixed.rda")
load("data/treeMixed.rda")
load("data/Xmixed.rda")

source("R/prior.R")
source("R/varchanges.R")
source("R/posterior.R")
source("R/utils.R")
source("R/rwm.R")
Rcpp::sourceCpp("src/mcmc.cpp")


PCMTreePlot(treeMixed, ladderize = TRUE, size = 2) + geom_tiplab(size = 4, offset = 0.01)

print(PCMTable(modelMixed))


# define priors
# priors placeholder
priors <- prior(modelMixed)
# populate priors
priors$X0 <- prior.uniform(min = -3, max = 3)
priors$H_1 <- priors$H_2 <- prior.uniform(min = 0, max = 10)
priors$Theta_1 <- priors$Theta_2 <- prior.uniform(min = 0, max = 5)
priors$Sigma_x_1 <- priors$Sigma_x_2 <- prior.uniform(min = 0, max = 3)


# run rwm
# set.seed(1234)
# timestamp()
# f1 <- rwm(model = modelMixed, X = XMixed, tree = treeMixed, priors = priors,
#           initial = c(0,1,1,1,1,1,1), nsteps = 100000, scale = 0.2, ncheck = 10)
# timestamp()
# load("check10.rda")
# f1 <- t(apply(P, 1, tr$g))
# mcmc_hist(f1)
# mcmc_trace(f1)


# ==============================================================================


priors2 <- prior(modelMixed)
priors2$X0 <- prior.normal(mean = 0, sd = 5)
priors2$H_1 <- priors2$H_2 <- prior.gamma(shape = 3, rate = 1)
priors2$Theta_1 <- priors2$Theta_2 <- prior.gamma(shape = 1, rate = 0.1)
priors2$Sigma_x_1 <- priors2$Sigma_x_2 <- prior.halfnormal(sigma = 5)

timestamp()
set.seed(1234)
f2 <- rwm(model = modelMixed, X = XMixed, tree = treeMixed, priors = priors2,
          initial = c(0,1,1,1,1,1,1), nsteps = 10000, scale = 0.1, ncheck = FALSE)
timestamp()
#save(f2, file = "f2.rda")
attr(f2, "accept")

xgrid = seq(0.001, 10, length.out = 1000)
plot(xgrid, priors2$H_1(xgrid), type = "l")
plot(xgrid, priors2$Theta_1(xgrid), type = "l")

mcmc_hist(f2)
mcmc_trace(f2)








