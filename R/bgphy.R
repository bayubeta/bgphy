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
  attr(res, "model") <- model
  attr(res, "X") <- X

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


# posterior predictive check
#' @export
post_pred_check <- function(post, nsample = 100, ...){
  # check if input is of class bgphy_posterior
  stopifnot("Object is not of class bgphy_posterior" = class(post) == "bgphy_posterior")

  N <- length(post$W)

  # retrieve model
  model <- attr(post, "model")

  # retrieve data
  X <- attr(post, "X")

  # number of tips of the tree
  ntips <- ape::Ntip(model$tree)

  # sample according to weights
  post_par <- post$Q[sample(1:N, nsample, prob = post$W, replace = TRUE),]

  # empty matric for storing Xrep (posterior predictive values)
  Xrep <- matrix(nrow = nsample, ncol = ncol(X))

  for (i in 1:nsample){
    # assign posterior params to model
    post_model <- setParams(post_par[i,], model$model)

    # simulate data using PCMBase, take only the tips
    Xrep[i,] <- PCMBase::PCMSim(model$tree, post_model, post_model$X0)[,1:ntips]
  }

  # ====================== plotting routines ======================
  # calculate maximum range for y
  max_y <- max(apply(Xrep, 1, function(x){density(x, n = 10)$y}))


  # base
  par(mar = c(3, 1, 1, 1))
  plot(density(X[1,]), main = "", xlab = "", ylab = "", yaxt = "n",
       col = rgb(0,0,0,0), ylim = c(0, 1.1*max_y))
  grid()

  # posterior predictive densities
  for (i in 1:nsample){
    lines(density(Xrep[i,]), col = rgb(0.2,0.2,1,0.07))
  }

  # data
  lines(density(X[1,]), lwd = 3, col = rgb(0,0,0.5,1))

  legend("topleft", legend = c("Data", "Simulated data"),
         col = c(rgb(0,0,0.5,1), rgb(0.2,0.2,1,0.5)),
         lty = c(1,1), lwd = c(2,1), seg.len = .75,
         bty = "n")

}

