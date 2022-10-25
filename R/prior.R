#' #' export
#' set_priors <- function(model, ...){
#'   UseMethod("set_priors")
#' }
#'
#' #' export
#' set_priors.PCM <- function(model){
#'   # set priors for a PCM object
#'
#'   # number of traits (dimension)
#'   k <- attr(model, "k")
#'
#'   # number of regimes
#'   r <- length(attr(model, "r"))
#'
#'   # number of parameters
#'   p <- attr(model, "p")
#'
#'   parnames <- PCMGetParamNames(model)
#'
#'   priors <- setNames(vector())
#'
#'   class(pars) <- "prior"
#'   return(pars)
#' }


#' @export
prior <- function(model, ...){
  UseMethod("prior")
}

prior.default <- function(model, ...){

  parnames <- PCMGetParamNames(model)

  structure(sapply(parnames, function(x) NULL), class = "prior")

}



#=================== pdfs of prior distributions ===================

prior.uniform <- function(min = 0, max = 1){
  f <- function(x){
    stats::dunif(x, min = min, max = max, log = TRUE)
  }

  class(f) <- c("uniform", "prior")
  attr(f, "bounds") = c(min, max)

  return(f)
}

prior.normal <- function(mean = 0, sd = 1){
  f <- function(x){
    stats::dnorm(x, mean = mean, sd = sd, log = TRUE)
  }

  class(f) <- c("normal", "prior")
  attr(f, "params") <- setNames(c(mean, sd), c("mean", "sd"))

  return(f)
}


prior.gamma <- function(shape, rate = 1, scale = 1/rate){
  f <- function(x){
    stats::dgamma(x, shape, rate = rate, log = TRUE)
  }

  class(f) <- c("gamma", "prior")
  attr(f, "bounds") <- c(0, Inf)
  attr(f, "params") <- setNames(c(shape, rate, scale), c("shape", "rate", "scale"))

  return(f)
}


prior.halfnormal <- function(sigma){
  f <- function(x){
    extraDistr::dhnorm(x, sigma = sigma, log = TRUE)
  }

  class(f) <- c("halfnormal", "prior")
  attr(f, "bounds") <- c(0, Inf)
  attr(f, "params") <- setNames(c(sigma), c("sigma"))

  return(f)
}

prior.halfcauchy <- function(sigma){
  f <- function(x){
    extraDistr::dhcauchy(x, sigma = sigma, log = TRUE)
  }

  class(f) <- c("halfcauchy", "prior")
  attr(f, "bounds") <- c(0, Inf)
  attr(f, "params") <- setNames(c(sigma), c("sigma"))

  return(f)
}



#=================== samplers of prior distributions ===================



priorsampler <- function(prior){
  UseMethod("priorsampler")
}

priorsampler.normal <- function(f){
  pars <- attr(f, "params")
  g <- function(n){
    stats::rnorm(n, mean = pars[1], sd = pars[2])
  }

  attr(g, "params") <- attr(f, "params")
  class(g) <- c("sampler", attr(f, "class"))

  return(g)
}

priorsampler.gamma <- function(f){
  pars <- attr(f, "params")
  g <- function(n){
    stats::rgamma(n, shape = pars[1], rate = pars[2])
  }

  attr(g, "params") <- attr(f, "params")
  attr(g, "class") <- c("sampler", attr(f, "class"))

  return(g)
}

priorsampler.halfnormal <- function(f){
  pars <- attr(f, "params")
  g <- function(n){
    extraDistr::rhnorm(n, sigma = pars)
  }

  attr(g, "params") <- attr(f, "params")
  attr(g, "class") <- c("sampler", attr(f, "class"))

  return(g)
}

priorsampler.halfcauchy <- function(f){
  pars <- attr(f, "params")
  g <- function(n){
    extraDistr::rhcauchy(n, sigma = pars)
  }

  attr(g, "params") <- attr(f, "params")
  attr(g, "class") <- c("sampler", attr(f, "class"))
  return(g)
}


prior_sampler <- function(model, ...){
  UseMethod("prior_sampler")
}

prior_sampler.prior <- function(priors){
  p_sampler <- lapply(priors, priorsampler)
  p <- length(p_sampler)
  function(n){
    mapply(function(f,n){f(n)}, p_sampler, n)
  }
}

