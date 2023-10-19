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
  priors <- structure(sapply(parnames, function(x) NULL), class = "bgphy_priors")

  # ---------------- fill priors with default priors ----------------
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

# print a list of priors of class bgphy_priors
#' @export
print.bgphy_priors <- function(bgphy_priors, ...){
  parnames <- names(bgphy_priors)

  for (par in parnames){
    cat(paste0(par))
    print(bgphy_priors[[par]])
    cat("\n")
  }
}


#============================= prior distributions PDF =============================

#' log Probability density functions.
#'
#' This is a collection of log probability density functions used for the prior object of class \code{mgpm_prior}.
#'
#' @return A log probability density function of class \code{priorpdf}.
#'
#' @examples
#' \dontrun{
#' modeltypes <- setNames(c("OU", "OU"), c("ancestral", "new"))
#' startNodes <- setNames(c(101, 135), c("ancestral", "new"))
#' OU <- setModel(tree = lizardTree, modeltypes = ("OU"))
#' # set custom priors
#' OU$priors$X0 <- prior_normal(mean = 0, sd = 2)
#' OU$priors$alpha <- prior_halfnormal(sigma = 3)
#' OU$priors$theta_1 <- prior_normal(mean = 2, sd = 2)
#' }
#'
#' @name priorpdf
#' @rdname priorpdf
#' @export
prior_uniform <- function(min = 0, max = 1){
  f <- function(x){
    stats::dunif(x, min = min, max = max, log = TRUE)
  }

  class(f) <- c("priorpdf", "uniform")
  attr(f, "bounds") <- c(min, max)
  attr(f, "params") <- stats::setNames(c(min, max), c("min", "max"))

  return(f)
}

#' @rdname priorpdf
#' @export
prior_normal <- function(mean = 0, sd = 1){
  f <- function(x){
    stats::dnorm(x, mean = mean, sd = sd, log = TRUE)
  }

  class(f) <- c("priorpdf", "normal")
  attr(f, "bounds") <- c(-Inf, Inf)
  attr(f, "params") <- stats::setNames(c(mean, sd), c("mean", "sd"))

  return(f)
}

#' @rdname priorpdf
#' @export
prior_gamma <- function(shape, rate = 1, scale = 1/rate){
  f <- function(x){
    stats::dgamma(x, shape, rate = rate, log = TRUE)
  }

  class(f) <- c("priorpdf", "gamma")
  attr(f, "bounds") <- c(0, Inf)
  attr(f, "params") <- stats::setNames(c(shape, rate, scale), c("shape", "rate", "scale"))

  return(f)
}

#' @rdname priorpdf
#' @export
prior_halfnormal <- function(sigma){
  f <- function(x){
    extraDistr::dhnorm(x, sigma = sigma, log = TRUE)
  }

  class(f) <- c("priorpdf", "halfnormal")
  attr(f, "bounds") <- c(0, Inf)
  attr(f, "params") <- stats::setNames(c(sigma), c("sigma"))

  return(f)
}

#' @rdname priorpdf
#' @export
prior_halfcauchy <- function(sigma){
  f <- function(x){
    extraDistr::dhcauchy(x, sigma = sigma, log = TRUE)
  }

  class(f) <- c("priorpdf", "halfcauchy")
  attr(f, "bounds") <- c(0, Inf)
  attr(f, "params") <- stats::setNames(c(sigma), c("sigma"))

  return(f)
}

#' @rdname priorpdf
#' @export
prior_halft <- function(scale, nu){
  f <- function(x){
    LaplacesDemon::dhalft(x, scale = scale, nu = nu)
  }

  class(f) <- c("priorpdf", "halft")
  attr(f, "bounds") <- c(0, Inf)
  attr(f, "params") <- stats::setNames(c(scale, nu), c("scale", "nu"))

  return(f)
}


# print single priorpdf
#' @export
print.priorpdf <- function(prior, ...){
  # print variable name and its distribution
  classes <- attr(prior, "class")
  type <- classes[length(classes)] # distribution name
  cat(" ~ ")
  cat(type, "(", sep = "")

  # print the parameter values
  params <- attr(prior, "params")
  parnames <- names(params)

  for (i in 1:length(parnames)){
    cat(parnames[i], " = ", params[[i]], sep = "")

    if (i != length(parnames)){
      cat(", ")
    }
    else{
      cat(")")
    }
  }
}





#============================= prior distributions RNG =============================
# ----------------------- To be used internally -----------------------

priorsampler <- function(priorpdf){
  UseMethod("priorsampler")
}

priorsampler.uniform <- function(priorpdf){
  pars <- attr(priorpdf, "bounds")
  g <- function(n){
    stats::runif(n, min = pars[1], max = pars[2])
  }

  attr(g, "params") <- attr(priorpdf, "bounds")
  class(g) <- c("sampler", attr(priorpdf, "class"))

  return(g)
}

priorsampler.normal <- function(priorpdf){
  pars <- attr(priorpdf, "params")
  g <- function(n){
    stats::rnorm(n, mean = pars[1], sd = pars[2])
  }

  attr(g, "params") <- attr(priorpdf, "params")
  class(g) <- c("sampler", attr(priorpdf, "class"))

  return(g)
}

priorsampler.gamma <- function(priorpdf){
  pars <- attr(priorpdf, "params")
  g <- function(n){
    stats::rgamma(n, shape = pars[1], rate = pars[2])
  }

  attr(g, "params") <- attr(priorpdf, "params")
  attr(g, "class") <- c("sampler", attr(priorpdf, "class"))

  return(g)
}

priorsampler.halfnormal <- function(priorpdf){
  pars <- attr(priorpdf, "params")
  g <- function(n){
    extraDistr::rhnorm(n, sigma = pars)
  }

  attr(g, "params") <- attr(priorpdf, "params")
  attr(g, "class") <- c("sampler", attr(priorpdf, "class"))

  return(g)
}

priorsampler.halfcauchy <- function(priorpdf){
  pars <- attr(priorpdf, "params")
  g <- function(n){
    extraDistr::rhcauchy(n, sigma = pars)
  }

  attr(g, "params") <- attr(priorpdf, "params")
  attr(g, "class") <- c("sampler", attr(priorpdf, "class"))
  return(g)
}

priorsampler.halft <- function(priorpdf){
  pars <- attr(priorpdf, "params")
  g <- function(n){
    LaplacesDemon::rhalft(n, scale = pars[1], nu = pars[2])
  }

  attr(g, "params") <- attr(priorpdf, "params")
  attr(g, "class") <- c("sampler", attr(priorpdf, "class"))
  return(g)
}


prior_sampler <- function(priors){
  p_sampler <- lapply(priors, priorsampler)
  p <- length(p_sampler)
  function(n){
    mapply(function(f,n){f(n)}, p_sampler, n)
  }
}

