#' @export
set_priors <- function(model, ...){
  UseMethod("set_priors")
}

#' @export
set_priors.PCM <- function(model){
  # set priors for a PCM object

  # number of traits (dimension)
  k <- attr(model, "k")

  # number of regimes
  r <- length(attr(model, "r"))

  # number of parameters
  p <- attr(model, "p")


  # 1-dimensional case
  if (k == 1){
    # if there is only a single (global) regime
    if (r == 1){
      # make a list whose elements mirror the model's parameters
      pars <- sapply(names(model)[-p], function(x) NULL)
      A <- 10
      pars$X0 <- pars$Theta <- prior.uniform(min = -A, max = A)
      pars$H <- pars$Sigma_x <- prior.uniform(min = 0, max = 2*A)
    }
  }

  class(pars) <- "prior"
  return(pars)
}



# uniform prior

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

prior.invgamma <- function(shape, rate = 1, scale = 1/rate){
  f <- function(x){
    invgamma::dinvgamma(x, shape, rate = rate, scale = scale, log = TRUE)
  }

  class(f) <- c("invgamma", "prior")
  attr(f, "bounds") <- c(0, Inf)
  attr(f, "params") <- setNames(c(shape, rate, scale), c("shape", "rate", "scale"))

  return(f)
}









