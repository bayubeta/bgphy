#' @export
prior <- function(model, ...){
  UseMethod("prior")
}

#' @export
prior.PCM <- function(model){
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
      pars$X0 <- pars$Theta <- function(x) stats::dunif(x, -A, A, log = TRUE)
      pars$H <- pars$Sigma_x <- function(x) stats::dunif(x, 0, 2*A, log = TRUE)
    }
  }

  class(pars) <- "prior"
  return(pars)
}







