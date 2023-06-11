#' Bayesian inference of mixed Gaussian phylogenetic models.
#'
#' The main function to run Bayesian inference of mixed Gaussian phylogenetic models.
#'
#' @param model An object of class \code{PCM}.
#' @param X Measurements on tips.
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @param priors A prior object of class \code{mgpm_prior}.
#' @param nsample Number of random draws to be made.
#' @param scale The scale parameter to determine the scale of the approximated covariance matrix
#'   obtained from the Laplace's approximation of the posterior distribution. The default is 1.
#' @param parallel Utilization of parallel processing. Defaults to \code{TRUE}.
#'
#' @return A list of class \code{mgpm_posterior} which contains:
#' * `E`:   Expected value of the posterior distribution.
#' * `Q`:   Values used to calculate the expected value using Importance Sampling.
#' * `W`:   Normalized weights for each row of `Q`.
#' * `WAIC`:   Widely Available Criterion Score.
#'
#' @examples
#' OU <- PCMBase::PCM("OU")
#' priorOU <- setPriors(OU)
#' priorOU$X0 <- prior_normal(mean = 0, sd = 1)
#' priorOU$H <- prior_halfnormal(sigma = 1)
#' priorOU$Theta <- prior_normal(mean = 0, sd = 1)
#' priorOU$Sigma_x <- prior_halfnormal(sigma = 1)
#'
#' post <- mgpm(OU, XOU[1,], lizardTree, priorOU, 100)
#'
#' @export
mgpm <- function(model, X, tree, priors, nsample = 10000, initial = NULL, scale = 1, parallel = TRUE){

  # initial position, sample from priors if not provided
  if (is.null(initial)){
    initial <- prior_sampler(priors)(1)
  }

  if (is.null(dim(X))){X <- matrix(X, nrow = 1)}

  res <- IS(model = model, X = X, tree = tree, priors = priors,
            initial = initial, nsample = nsample, scale = scale, parallel = parallel)

  class(res) <- "mgpm_posterior"

  return(res)
}
