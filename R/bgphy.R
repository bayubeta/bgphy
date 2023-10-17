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
#' # OU <- PCMBase::PCM("OU")
#' # priorOU <- setPriors(OU)
#' # priorOU$X0 <- prior_normal(mean = 0, sd = 1)
#' # priorOU$H <- prior_halfnormal(sigma = 1)
#' # priorOU$Theta <- prior_normal(mean = 0, sd = 1)
#' # priorOU$Sigma_x <- prior_halfnormal(sigma = 1)
#' #
#' # post <- mgpm(OU, XOU[1,], lizardTree, priorOU, 100)
#'
#' @export
bgphy <- function(model, X, nsample = 10000, scale = 1, parallel = TRUE){

  # run inference using importance sampling
  res <- IS(model = model, X = X, nsample = nsample, scale = scale, parallel = parallel)

  # matrix of normalized weights
  MW <- matrix(rep(res$W, dim(res$Q)[2]), ncol = dim(res$Q)[2])

  # compute mean
  res$mean <- colSums(res$Q*MW)

  # compute std error of posterior parameters mean
  # matrix of means in each row
  Mm <- matrix(rep(res$mean, nsample), nrow = dim(res$Q)[1], ncol = dim(res$Q)[2], byrow = TRUE)
  var_lp <- colSums((MW^2)*((res$Q - Mm)^2))
  res$std_error <- stats::setNames(sqrt(var_lp), names(res$mean))

  # compute (estimated) standard deviation
  var_hat <- colSums(res$Q^2 * MW) - (res$mean)^2
  res$std <- stats::setNames(sqrt(var_hat), names(res$mean))

  # matrix of quantiles
  Mq <- est_quantiles(res$Q, res$W)
  rownames(Mq) <- names(res$mean)
  res$quantiles <- Mq

  # change order of the list
  res <- res[c("mean", "std_error", "std", "quantiles", "WAIC", "Q", "W")]


  class(res) <- "bgphy_posterior"

  return(res)
}


#' @export
print.bgphy_posterior <- function(post, ...){
  # get parameter names
  par_names <- names(post$mean)

  # get quantity names
  q_names <- c(names(post)[1:3], c("2.5%", "50%", "97.5%"))

  # put information into a matrix
  M <- matrix(unlist(post[1:4]), nrow = length(par_names))
  colnames(M) <- q_names
  rownames(M) <- par_names

  print(M)

  cat("\n")

  cat(paste0("WAIC: ", post$WAIC))
}

