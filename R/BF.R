

# naive Monte Carlo integral approximation

# marginal likelihood
#' @export
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

