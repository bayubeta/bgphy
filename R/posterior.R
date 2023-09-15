setParams <- function(p, model){
  # p: all parameters, excluding Sigmae_x
  # p := c(X0, H1, Theta1, Sigma_x1, ..., Hr, Thetar, Sigma_xr, Sigmae_x)
  # Load parameters into the model

  # create an empty vector of the size of total parameters count
  v <- numeric(attr(model, "p"))

  # fill everything except the last with p
  v[-length(v)] <- p

  # load v into the model
  model <- loadParams(model, v)

  return(model)
}


# log-likelihood given data, tree, and model with parameter values
PCMloglik <- function(X, tree, model, p){
  # model with parameters p
  prop_model <- setParams(p, model)

  # log likelihood calculation
  loglik <- PCMBase::PCMLik(X, tree, prop_model, metaI = PCMBaseCpp::PCMInfoCpp)

  return(loglik[1])
}

# log-unnormalized-posterior
lupost <- function(p, model, X, tree, priors_tr, tr){

  # p: vector of parameters, c(X0, H1, Theta1, Sigma_x1, ..., Hr, Thetar, Sigma_xr)
  # model: PCM model
  # tree: phylo object
  # priors_tr: list of priors pdf on unbounded space
  # tr: transformation functions, obtained from trfunc()

  if (is.null(dim(X))){X <- matrix(X, nrow = 1)}

  # change parameters to p_b(ounded),
  # because PCMLik (the likelihood) needs parameter values from the original space
  p_b <- tr$g(p)

  # do not raise warning
  op <- options(PCMBase.Raise.Lik.Errors = FALSE)
  on.exit(options(op))

  # log likelihood
  loglik <- PCMloglik(X, tree, model, p_b)

  # log priors
  # calculated in the unbounded space
  log_p <- sum(mapply(function(f, x){f(x)}, priors_tr, p))

  # sum of log-likelihood and log priors
  sum_log <- loglik + log_p

  # if infinite, return a very low value
  if (is.infinite(sum_log) || is.na(sum_log)){
    return(-1e20)
  }

  return(stats::setNames(c(sum_log[1],loglik), c("log_u_post", "loglik")))
}

