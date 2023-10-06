#' @export
prior <- function(model, ...){
  UseMethod("prior")
}

#' Set a prior object given a model of class \code{PCM}.
#'
#' This function takes a model of class \code{PCM}, created from \code{PCMBase},
#'   and outputs a blank prior object of class \code{mgpm_prior}.
#'   Users can then modify the prior accordingly using functions of class \code{priorpdf}.
#'   Modifying an \code{mgpm_prior} object is similar to modifying a list,
#'   in which the users can replace the values in the list with \code{priorpdf} functions
#'   provided by the \code{bgphy} package or by the users.
#'
#' @param model A \code{PCM} object.
#'
#' @return A prior object of class \code{mgpm_prior}, which is required for the main function, \code{mgpm}.
#'
#' @examples
#' OU <- PCMBase::PCM("OU")
#' priorOU <- setPriors(OU)
#' priorOU$X0 <- prior_normal(mean = 0, sd = 1)
#' priorOU$H <- prior_halfnormal(sigma = 1)
#' priorOU$Theta <- prior_normal(mean = 0, sd = 1)
#' priorOU$Sigma_x <- prior_halfnormal(sigma = 1)
#' @export
# setPriors <- function(model){
#   parnames <- PCMGetParamNames(model)
#
#   structure(sapply(parnames, function(x) NULL), class = "mgpm_prior")
# }


#' @export
prior.default <- function(model, ...){

  parnames <- PCMGetParamNames(model)

  structure(sapply(parnames, function(x) NULL), class = "prior")

}


#' @export
print.mgpm_prior <- function(priors, ...){

  varnames <- names(priors)

  for (i in 1:length(varnames)){
    cat("  ", varnames[i], " ~ ", sep = "")
    print(priors[[i]], unit = FALSE)
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
#' OU <- PCMBase::PCM("OU")
#' priorOU <- setPriors(OU)
#' priorOU$X0 <- prior_normal(mean = 0, sd = 1)
#' priorOU$H <- prior_halfnormal(sigma = 1)
#' priorOU$Theta <- prior_normal(mean = 0, sd = 1)
#' priorOU$Sigma_x <- prior_halfnormal(sigma = 1)
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


#' @export
priorPDF <- function(f, bounds){
  fname <- deparse(substitute(f))
  class(f) <- c("priorpdf", "custom")
  attr(f, "bounds") <- bounds
  attr(f, "fname") <- fname
  return(f)
}


# print priorpdf class
#' @export
print.priorpdf <- function(priorpdf, unit = TRUE, ...){
  # print variable name and its distribution
  classes <- attr(priorpdf, "class")
  type <- classes[length(classes)]

  if (unit){
    cat("  ~ ")
  }

  if (type == "custom"){
    fname <- attr(priorpdf, "fname")
    cat(fname)
  }else{
    cat(type, "(", sep = "")

    # print the parameter values
    params <- attr(priorpdf, "params")
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
}






#============================= prior distributions RNG =============================


#' @export
priorsampler <- function(priorpdf){
  UseMethod("priorsampler")
}

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


#' @export
prior_sampler <- function(priors){
  UseMethod("prior_sampler")
}


#' @export
prior_sampler.mgpm_prior <- function(priors){
  p_sampler <- lapply(priors, priorsampler)
  p <- length(p_sampler)
  function(n){
    mapply(function(f,n){f(n)}, p_sampler, n)
  }
}

