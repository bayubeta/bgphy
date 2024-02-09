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
    # sigma > 0, 3*tree_height
    priors$sigma <- prior_halft(nu = 1, sigma = 3*t_height)

    if (modeltypes == "OU"){
      # alpha > 0, t_{1/2} = 0.05*t_height => log(2)/(0.05 * t_height)
      s_alpha <- log(2)/(0.05 * t_height * 2)
      priors$alpha <- prior_halfnormal(sigma = s_alpha)

      # theta, same scale as X_0
      priors$theta <- prior_normal(mean = 0, sd = 10)
    }
  }
  else{
    # for each regime
    for (i in 1:r){
      # sigma > 0, 3*tree_height
      s_sigma <- 3*t_height
      priors[[paste0("sigma_", i)]] <- prior_halft(nu = 1, sigma = 3*t_height)

      if (modeltypes[i] == "OU"){
        # alpha > 0, t_{1/2} = 0.05*t_height => log(2)/(0.05 * t_height)
        s_alpha <- log(2)/(0.05 * t_height * 2)
        priors[[paste0("alpha_", i)]] <- prior_halfnormal(sigma = s_alpha)

        # theta, same scale as X_0
        priors[[paste0("theta_", i)]] <- prior_normal(mean = 0, sd = 10)
      }
    }
  }

  return(priors)
}

# print a list of priors of class bgphy_priors
#' Print priors information
#'
#' Prints objects of class `bgphy_priors`, which contains objects of class `priorpdf`.
#'
#' @param x An object of class `bgphy_priors`.
#' @param ... Other arguments to be passed to `print` or `print.priorpdf`
#'
#' @examples
#' \dontrun{
#' # global model
#' OU1 <- setModel(tree = lizardTree, regime_names = "Regime1", modeltypes = "OU")
#' OU1$priors
#'
#' # mixed, BM to OU
#' BMOU <- setModel(tree = lizardTree,
#'                  regime_names = c("Ancestral", "New"),
#'                  modeltypes = c("BM", "OU"),
#'                  startNodes = list(Ancestral = c(101), New = c(135)))
#' BMOU$priors
#' }
#'
#' @export
print.bgphy_priors <- function(x, ...){
  parnames <- names(x)

  for (par in parnames){
    cat(paste0(par))
    print(x[[par]], ...)
    cat("\n")
  }
}


#============================= prior distributions PDF =============================

#' log Probability density functions.
#'
#' This is a collection of log probability density functions used for the prior object of class \code{bgphy_priors}.
#'
#' @return A log probability density function of class \code{priorpdf}.
#'
#' @param min Lower bound.
#' @param max Upper bound.
#' @param mean Mean.
#' @param sd Standard deviation.
#' @param shape Shape parameter.
#' @param rate Rate parameter.
#' @param sigma Scale parameter.
#' @param nu Degrees of freedom.
#'
#'
#' @examples
#' \dontrun{
#' # global model
#' OU1 <- setModel(tree = lizardTree, regime_names = "Regime1", modeltypes = "OU")
#' # set custom priors
#' OU1$priors$X0 <- OU1$priors$theta <- prior_normal(mean = 0, sd = 3)
#' OU1$priors$alpha <- OU1$priors$sigma <- prior_halfnormal(sigma = 2)
#'
#' # mixed, BM to OU
#' BMOU <- setModel(tree = lizardTree,
#'                  regime_names = c("Ancestral", "New"),
#'                  modeltypes = c("BM", "OU"),
#'                  startNodes = list(Ancestral = c(101), New = c(135)))
#' # set custom priors
#' BMOU$priors$X0 <- prior_normal(mean = 0, sd = 3)
#' BMOU$priors$sigma_1 <- prior_halft(nu = 1, sigma = 2)
#' BMOU$priors$alpha_2 <- prior_halfnormal(sigma = 5)
#' BMOU$priors$theta_2 <- prior_normal(mean = 3, sd = 5)
#' BMOU$priors$sigma_2 <- prior_halfnormal(sigma = 3)
#' }
#'
#' @name priorpdf
#' @rdname priorpdf
#' @export
prior_uniform <- function(min = 0, max = 1){
  stopifnot("Parameter values are invalid" = is.numeric(c(min, max)) & (min < max))
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
  stopifnot("Parameter values are invalid" = is.numeric(c(mean, sd)) & (sd>0))
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
prior_gamma <- function(shape, rate = 1){
  stopifnot("Parameter values are invalid" = is.numeric(c(shape, rate)) &
              all(c(shape, rate)>0))
  f <- function(x){
    stats::dgamma(x, shape, rate = rate, log = TRUE)
  }

  class(f) <- c("priorpdf", "gamma")
  attr(f, "bounds") <- c(0, Inf)
  attr(f, "params") <- stats::setNames(c(shape, rate), c("shape", "rate"))

  return(f)
}

#' @rdname priorpdf
#' @export
prior_halfnormal <- function(sigma){
  stopifnot("Parameter values are invalid" = is.numeric(sigma) & (sigma>0))
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
  stopifnot("Parameter values are invalid" = is.numeric(sigma) & (sigma>0))
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
prior_halft <- function(nu, sigma){
  stopifnot("Parameter values are invalid" = is.numeric(c(nu, sigma)) & all(c(nu, sigma)>0))
  f <- function(x){
    extraDistr::dht(x, nu = nu, sigma = sigma, log = TRUE)
  }

  class(f) <- c("priorpdf", "halft")
  attr(f, "bounds") <- c(0, Inf)
  attr(f, "params") <- stats::setNames(c(nu, sigma), c("nu", "sigma"))

  return(f)
}


# print single priorpdf
#' Print single prior information
#'
#' Prints objects of class `priorpdf`.
#'
#' @param x An object of class `priorpdf`.
#' @param ... Other arguments to be passed to `print` or `print.priorpdf`
#'
#' @examples
#' \dontrun{
#' # global model
#' OU1 <- setModel(tree = lizardTree, regime_names = "Regime1", modeltypes = "OU")
#' OU1$priors$alpha
#'
#' # mixed, BM to OU
#' BMOU <- setModel(tree = lizardTree,
#'                  regime_names = c("Ancestral", "New"),
#'                  modeltypes = c("BM", "OU"),
#'                  startNodes = list(Ancestral = c(101), New = c(135)))
#' BMOU$priors$sigma_2
#' }
#'
#' @export
print.priorpdf <- function(x, ...){
  # print variable name and its distribution
  classes <- attr(x, "class")
  type <- classes[length(classes)] # distribution name
  cat(" ~ ")
  cat(type, "(", sep = "")

  # print the parameter values
  params <- attr(x, "params")
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

#' RNG samplers for prior distribution.
#'
#' Converts a `priorpdf` object into an RNG sampler.
#'
#' @param priorpdf An object of class `priorpdf`.
#'
#' @return A sampler of class `sampler` with subclass that depends on the type
#'         of distribution of the `priorpdf` object. This sampler only takes number
#'         of samples `n` as its argument, while the distributional parameters are
#'         fully specified by the `priorpdf` object.
#'
#' @examples
#' \dontrun{
#' OU1 <- setModel(tree = lizardTree, regime_names = "Regime1", modeltypes = "OU")
#' ralpha <- priorsampler(OU1$priors$alpha)
#' ralpha(100)
#' }
#'
#' @name priorsampler
#' @rdname priorsampler
#' @export
priorsampler <- function(priorpdf){
  UseMethod("priorsampler")
}

#' @rdname priorsampler
#' @export
priorsampler.uniform <- function(priorpdf){
  pars <- attr(priorpdf, "bounds")
  g <- function(n){
    stats::runif(n, min = pars[1], max = pars[2])
  }

  attr(g, "params") <- attr(priorpdf, "bounds")
  class(g) <- c("sampler", attr(priorpdf, "class"))

  return(g)
}

#' @rdname priorsampler
#' @export
priorsampler.normal <- function(priorpdf){
  pars <- attr(priorpdf, "params")
  g <- function(n){
    stats::rnorm(n, mean = pars[1], sd = pars[2])
  }

  attr(g, "params") <- attr(priorpdf, "params")
  class(g) <- c("sampler", attr(priorpdf, "class"))

  return(g)
}

#' @rdname priorsampler
#' @export
priorsampler.gamma <- function(priorpdf){
  pars <- attr(priorpdf, "params")
  g <- function(n){
    stats::rgamma(n, shape = pars[1], rate = pars[2])
  }

  attr(g, "params") <- attr(priorpdf, "params")
  attr(g, "class") <- c("sampler", attr(priorpdf, "class"))

  return(g)
}

#' @rdname priorsampler
#' @export
priorsampler.halfnormal <- function(priorpdf){
  pars <- attr(priorpdf, "params")
  g <- function(n){
    extraDistr::rhnorm(n, sigma = pars)
  }

  attr(g, "params") <- attr(priorpdf, "params")
  attr(g, "class") <- c("sampler", attr(priorpdf, "class"))

  return(g)
}

#' @rdname priorsampler
#' @export
priorsampler.halfcauchy <- function(priorpdf){
  pars <- attr(priorpdf, "params")
  g <- function(n){
    extraDistr::rhcauchy(n, sigma = pars)
  }

  attr(g, "params") <- attr(priorpdf, "params")
  attr(g, "class") <- c("sampler", attr(priorpdf, "class"))
  return(g)
}

#' @rdname priorsampler
#' @export
priorsampler.halft <- function(priorpdf){
  pars <- attr(priorpdf, "params")
  g <- function(n){
    extraDistr::rht(n, nu = pars[1], sigma = pars[2])
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

