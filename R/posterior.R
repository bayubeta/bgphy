
setParams <- function(p, model){
  # p: all parameters, excluding Sigmae_x
  # p := c(X0, H1, Theta1, Sigma_x1, ..., Hr, Thetar, Sigma_xr, Sigmae_x)
  # Load parameters into the model

  # create an empty vector of the size of total parameters count
  v <- numeric(PCMBase::PCMParamCount(model))

  # fill everything except the last with p
  v[-length(v)] <- p

  # load v into the model
  PCMBase::PCMParamLoadOrStore(model, v, offset = 0, load = TRUE)

  return(model)
}






lupost <- function(p, model, X, tree, priors_tr, tr){

  # p: vector of parameters, c(X0, H1, Theta1, Sigma_x1, ..., Hr, Thetar, Sigma_xr)
  # model: PCM model
  # tree: phylo object
  # priors_tr: list of priors pdf on unbounded space
  # tr: transformation functions

  # change parameters to p (bounded), because PCMLik needs parameter values from the original space
  p_b <- tr$g(p)

  # proposed model
  prop_model <- setParams(p_b, model)

  # do not raise warning
  op <- options(PCMBase.Raise.Lik.Errors = FALSE)
  on.exit(options(op))

  # log likelihood
  loglik <- PCMBase::PCMLik(X, tree, prop_model, metaI = PCMInfoCpp)

  # log priors
  log_p <- sum(mapply(function(f, x){f(x)}, priors_tr, p))

  # sum of log-likelihood and log priors
  sum_log <- loglik + log_p

  # if infinite, return a very low value
  if (is.infinite(sum_log) || is.na(sum_log)){
    return(-1e20)
  }

  return(sum_log[1])
}


lupost_k <- function(p, y, model, X, tree, priors_tr, tr){

  # p: vector of parameters, c(X0, H1, Theta1, Sigma_x1, ..., Hr, Thetar, Sigma_xr)
  # y: exponent
  # model: PCM model
  # tree: phylo object
  # priors_tr: list of priors pdf on unbounded space
  # tr: transformation functions

  # change parameters to p (bounded), because PCMLik needs parameter values from the original space
  p_b <- tr$g(p)

  # proposed model
  prop_model <- setParams(p_b, model)

  op <- options(PCMBase.Raise.Lik.Errors = FALSE)
  on.exit(options(op))

  # log likelihood
  loglik <- y*PCMBase::PCMLik(X, tree, prop_model, metaI = PCMInfoCpp)

  # log priors
  log_p <- sum(mapply(function(f, x){f(x)}, priors_tr, p))

  # sum of log-likelihood and log priors
  sum_log <- loglik + log_p

  # if too small or too large, return infinite
  if (is.infinite(sum_log) || is.na(sum_log)){
    return(Inf)
  }

  return(sum_log[1])
}










