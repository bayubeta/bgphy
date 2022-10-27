#' @export
prior <- function(model, ...){
  UseMethod("prior")
}

#' @export
prior.default <- function(model, ...){

  parnames <- PCMGetParamNames(model)

  structure(sapply(parnames, function(x) NULL), class = "prior")

}

#' @export
print.prior <- function(priors){

  varnames <- names(priors)

  for (i in 1:length(varnames)){

    # print variable name and its distribution
    type <- attr(priors[[i]], "class")[1]
    cat(varnames[i], " ~ ", type, "(", sep = "")

    # print the parameter values
    params <- attr(priors[[i]], "params")
    parnames <- names(params)
    for (j in 1:length(parnames)){
      cat(parnames[j], " = ", params[[j]], sep = "")

      if (j != length(parnames)){
        cat(", ")
      }
      else{
        cat(")")
      }
    }
    cat("\n")
  }
}


#=================== pdfs of prior distributions ===================

#' @export
prior_uniform <- function(min = 0, max = 1){
  f <- function(x){
    stats::dunif(x, min = min, max = max, log = TRUE)
  }

  class(f) <- c("uniform", "priorpdf")
  attr(f, "bounds") = c(min, max)
  attr(f, "params") <- stats::setNames(c(min, max), c("min", "max"))

  return(f)
}

#' @export
prior_normal <- function(mean = 0, sd = 1){
  f <- function(x){
    stats::dnorm(x, mean = mean, sd = sd, log = TRUE)
  }

  class(f) <- c("normal", "prior")
  attr(f, "params") <- stats::setNames(c(mean, sd), c("mean", "sd"))

  return(f)
}

#' @export
prior_gamma <- function(shape, rate = 1, scale = 1/rate){
  f <- function(x){
    stats::dgamma(x, shape, rate = rate, log = TRUE)
  }

  class(f) <- c("gamma", "prior")
  attr(f, "bounds") <- c(0, Inf)
  attr(f, "params") <- stats::setNames(c(shape, rate, scale), c("shape", "rate", "scale"))

  return(f)
}


#' @export
prior_halfnormal <- function(sigma){
  f <- function(x){
    extraDistr::dhnorm(x, sigma = sigma, log = TRUE)
  }

  class(f) <- c("halfnormal", "prior")
  attr(f, "bounds") <- c(0, Inf)
  attr(f, "params") <- stats::setNames(c(sigma), c("sigma"))

  return(f)
}


#' @export
prior_halfcauchy <- function(sigma){
  f <- function(x){
    extraDistr::dhcauchy(x, sigma = sigma, log = TRUE)
  }

  class(f) <- c("halfcauchy", "prior")
  attr(f, "bounds") <- c(0, Inf)
  attr(f, "params") <- stats::setNames(c(sigma), c("sigma"))

  return(f)
}



#=================== samplers of prior distributions ===================


#' @export
priorsampler <- function(prior){
  UseMethod("priorsampler")
}

#' @export
priorsampler.uniform <- function(prior){
  pars <- attr(prior, "bounds")
  g <- function(n){
    log(stats::runif(n, min = pars[1], max = pars[2]))
  }

  attr(g, "params") <- attr(prior, "bounds")
  class(g) <- c("sampler", attr(prior, "class"))

  return(g)
}


#' @export
priorsampler.normal <- function(prior){
  pars <- attr(prior, "params")
  g <- function(n){
    stats::rnorm(n, mean = pars[1], sd = pars[2])
  }

  attr(g, "params") <- attr(prior, "params")
  class(g) <- c("sampler", attr(prior, "class"))

  return(g)
}

#' @export
priorsampler.gamma <- function(prior){
  pars <- attr(prior, "params")
  g <- function(n){
    stats::rgamma(n, shape = pars[1], rate = pars[2])
  }

  attr(g, "params") <- attr(prior, "params")
  attr(g, "class") <- c("sampler", attr(prior, "class"))

  return(g)
}

#' @export
priorsampler.halfnormal <- function(prior){
  pars <- attr(prior, "params")
  g <- function(n){
    extraDistr::rhnorm(n, sigma = pars)
  }

  attr(g, "params") <- attr(prior, "params")
  attr(g, "class") <- c("sampler", attr(prior, "class"))

  return(g)
}

#' @export
priorsampler.halfcauchy <- function(prior){
  pars <- attr(prior, "params")
  g <- function(n){
    extraDistr::rhcauchy(n, sigma = pars)
  }

  attr(g, "params") <- attr(prior, "params")
  attr(g, "class") <- c("sampler", attr(prior, "class"))
  return(g)
}


#' @export
prior_sampler <- function(priors){
  UseMethod("prior_sampler")
}


#' @export
prior_sampler.prior <- function(priors){
  p_sampler <- lapply(priors, priorsampler)
  p <- length(p_sampler)
  function(n){
    mapply(function(f,n){f(n)}, p_sampler, n)
  }
}

