
priorsamplers <- lapply(priors, priorsampler)

# likelihood given parameter values
PCMloglik <- function(X, tree, model, p){
  # model with parameters p
  prop_model <- setParams(p, model)

  # log likelihood calculation
  loglik <- PCMBase::PCMLik(X, tree, prop_model, metaI = PCMInfoCpp)

  return(loglik[1])
}


# marginal likelihood
logMargLik <- function(X, tree, model, priors, nsamples){

  # create g(theta), priors sampler
  g <- lapply(priors, priorsampler)

  # parameter names
  names <- names(g)

  # matrix of prior samples
  P <- matrix(unlist(lapply(g, function(f){f(nsamples)})), nrow = nsamples)
  colnames(P) <- names

  op <- options(PCMBase.Raise.Lik.Errors = FALSE)
  on.exit(options(op))



  loglikP <- apply(P, 1,
                    FUN = function(p){
                      PCMloglik(X = X, tree = tree, model = model, p)
                      })

  # apply logsumexp
  lmax <- max(loglikP)
  loglikP_shifted <- loglikP - lmax
  return(-log(nsamples) + lmax + log(sum(exp(loglikP_shifted))))
}



# marginal likelihood
logMargLik <- function(X, tree, model, priors, nsamples){

  # create g(theta), priors sampler
  g <- lapply(priors, priorsampler)

  # parameter names
  names <- names(g)

  # vector of log-likelihood
  loglikv <- numeric(nsamples)

  op <- options(PCMBase.Raise.Lik.Errors = FALSE)
  on.exit(options(op))

  for (i in 1:nsamples){
    repeat{
      # sample prior
      pr <- sapply(g, function(f){f(1)})
      # likelihood from priors
      loglik <- PCMloglik(X = X, tree = tree, model = model, pr)
      if (!is.na(loglik)){
        break()
      }
    }
    loglikv[i] <- loglik
  }


  # apply logsumexp
  lmax <- max(loglikv)
  loglikv_shifted <- loglikv - lmax
  return(-log(nsamples) + lmax + log(sum(exp(loglikv_shifted))))
}




set.seed(1234)
lm1 <- logMargLik(X, tree, model, priors, 10000)
lm2 <- logMargLik(X, tree2, model2, priors2, 10000)
lm3 <- logMargLik(X, tree2, model3, priors3, 10000)
LM <- c(lm1, lm2, lm3)
save(LM, file = "LM.rda")


lm1f <- logMargLik(Xf, tree, model, priors, 10000)
lm2f <- logMargLik(Xf, tree2, model2, priors2, 10000)
lm3f <- logMargLik(Xf, tree2, model3, priors3, 10000)
LMf <- c(lm1f, lm2f, lm3f)
save(LMf, file = "LMf.rda")

LM[1]-LM[2]
LM[1]-LM[3]
LM[2]-LM[3]

LMf[1]-LMf[2]
LMf[1]-LMf[3]
LMf[2]-LMf[3]
