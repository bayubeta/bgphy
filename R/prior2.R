


# set default priors given the model
defaultPriors <- function(model){
  # calculate tree height
  t_height <- max(ape::node.depth.edgelength(model$tree))

  # retrieve number of regimes
  rnames <- attr(model$model, "regimes") # get regime names
  r <- length(rnames) # get number of regimes

  # model types
  modeltypes <- attr(model, "modeltypes")

  # parameter names
  parnames <- getParamNames(model)

  # empty list of priors
  priors <- structure(sapply(parnames, function(x) NULL), class = "bgphy_prior")

  # ---- fill priors with default priors ----
  # for global X_0
  priors$X0 <- prior_normal(mean = 0, sd = 10)

  if (r == 1){
    # for one regime
    # alpha > 0, t_{1/2} = 0.05*t_height => log(2)/(0.05 * t_height)
    s_alpha <- log(2)/(0.05 * t_height * 2)
    priors$alpha <- prior_halfnormal(sigma = s_alpha)

    # theta, same scale as X_0
    priors$theta <- prior_normal(mean = 0, sd = 10)

    # sigma > 0, 3*tree_height
    priors$sigma <- prior_halft(scale = 3*t_height, nu = 1)
  }
  else{
    # for each regime
    for (i in 1:r){
      # alpha > 0, t_{1/2} = 0.05*t_height => log(2)/(0.05 * t_height)
      s_alpha <- log(2)/(0.05 * t_height * 2)
      priors[[paste0("alpha_", i)]] <- prior_halfnormal(sigma = s_alpha)

      # theta, same scale as X_0
      priors[[paste0("theta_", i)]] <- prior_normal(mean = 0, sd = 10)

      # sigma > 0, 3*tree_height
      s_sigma <- 3*t_height
      priors[[paste0("sigma_", i)]] <- prior_halft(scale = 3*t_height, nu = 1)
    }
  }

  return(priors)
}
