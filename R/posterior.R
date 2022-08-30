setParams = function(Pars, model){
  # Pars: a vector of 4 (X0, H, Theta, Sigma_x)

  prop_model <- PCM(model = class(model)[1], k = attr(model, "k"))
  prop_model$X0[1] <- Pars[1]
  prop_model$H[,,1] = Pars[2]
  prop_model$Theta[1] = Pars[3]
  prop_model$Sigma_x[,,1] = Pars[4]

  return(prop_model)
}


lupost_factory <- function(model, X, tree, prior){
  function(Pars){
    # log-likelihood
    loglik <- PCMBase::PCMLik(X, tree, setParams(Pars, model))

    # log priors
    lp_X0 <- prior$X0(Pars[1])
    lp_alpha <- prior$H(Pars[2])
    lp_theta <- prior$Theta(Pars[3])
    lp_sigma_x <- prior$Sigma_x(Pars[4])

    return(loglik + lp_X0 + lp_alpha + lp_theta + lp_sigma_x)
  }
}

#' @export
bpcm <- function(model, X, tree, prior = NULL, method = "MH"){
  UseMethod("bpcm")
}

#' @export
bpcm.PCM <- function(model, X, tree, prior = NULL, method = "MH"){
  if (is.null(prior)){
    prior <- prior(model)
  }

  lupost <- lupost_factory(model, X, tree, prior)

  pars_init <- c(1,1,1,1)
  out = metrop(lupost, pars_init, nbatch = 1000, scale = 1)

  return(out)
}


