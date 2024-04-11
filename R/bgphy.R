#' Bayesian inference of mixed Gaussian phylogenetic models.
#'
#' The main function to run Bayesian inference of mixed Gaussian phylogenetic models.
#'
#' @param model An object of class \code{bgphy_model}, which is created by [bgphy::setModel()] function.
#' @param X Measurements on tips. \code{X} needs to be a matrix of dimension \code{(1 x N)}, where \code{N} is the number of tips on the tree in \code{model$tree}.
#'   The column names of \code{X} needs to match the names on the tips.
#' @param nsample Number of random draws to be made.
#' @param scale The scale parameter to determine the scale of the approximated covariance matrix
#'   obtained from the Laplace's approximation of the posterior distribution. The default is 1.
#' @param parallel Utilization of parallel processing. Defaults to \code{TRUE}.
#'
#' @return A list of class \code{mgpm_posterior} which contains:
#' * `mean`:   Expected value of the posterior distribution.
#' * `std_error`:   Estimated standard error of the posterior mean.
#' * `std`:   Estimated standard deviation of the posterior distribution.
#' * `quantiles`:  Estimated quantiles of marginal posterior distribution.
#' * `Q`:   Drawn values from the proposal distribution (multivariate Normal).
#' * `W`:   Normalized weights for each row of `Q`.
#' * `loss`:   Estimated posterior predictive loss.
#' * `ESS`:   Estimated effective sample size.
#'
#' @examples
#' \dontrun{
#  # global model
#' OU1 <- setModel(tree = lizardTree, regime_names = "Regime1", modeltypes = "OU")
#' post_OU <- bgphy(OU1,
#'                  matrix(XOU[1,], nrow = 1,
#'                         dimnames = list(NULL, colnames(XOU))))
#'
#' # mixed, OU to OU
#' OUOU <- setModel(tree = lizardTree,
#'                  regime_names = c("Ancestral", "New"),
#'                  modeltypes = c("OU", "OU"),
#'                  startNodes = list(Ancestral = c(101), New = c(135)))
#' post_OUOU <- bgphy(BMOU,
#'                    matrix(XMixedOU[1,], nrow = 1,
#'                           dimnames = list(NULL, colnames(XMixedOU))))
#' }
#'
#' @export
bgphy <- function(model, X, nsample = 10000, scale = 1, parallel = TRUE){

  # ------------------------- perform checks on inputs -------------------------
  parnames <- getParamNames(model)
  priors_class <- sapply(model$priors, function(x){class(x)[1]})
  # check if all priors are of class priorpdf
  stopifnot("Priors are not of class priorpdf." = all(priors_class == "priorpdf"))
  # check if there is a prior for each parameter
  # stopifnot("Priors are not defined for all parameters." = length(parnames) == length(priors_class))
  # check if X is a numeric matrix
  stopifnot("X is not a numeric matrix." = is.numeric(X) & is.matrix(X))
  # check if X has the correct dimension
  # number of tips
  ntips <- ape::Ntip(model$tree)
  stopifnot("X is not of a correct dimension." = all(dim(X) == c(1, ntips)))
  # match tips' names with X names
  stopifnot("Names on the tips do not match column names of X." = all(colnames(X) == model$tree$tip.label))
  # make sure nsample integer > 0
  stopifnot("Invalid value of nsample." = is.numeric(nsample) & ((abs(nsample - round(nsample)) < .Machine$double.eps^0.5)) & (nsample > 0))
  # make sure scale > 0
  stopifnot("Invalid value of scale" = is.numeric(scale) & (scale > 0))
  # make sure parallel is boolean
  stopifnot("Invalid value of parallel" = is.logical(parallel))


  # ---------------------------- start inference -------------------------------
  # run inference using importance sampling
  res <- IS(model = model, X = X, nsample = nsample, scale = scale, parallel = parallel)
  # ----------------------------------------------------------------------------

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

  # matrix of marginal quantiles
  Mq <- est_quantiles(res$Q, res$W)
  rownames(Mq) <- names(res$mean)
  res$quantiles <- Mq

  # calculate posterior predictive loss
  ppred <- ppred_loss(res, X, model)
  res$loss <- ppred$loss

  # calculate ESS
  res$ESS <- 1/sum(res$W^2)

  # change order of the list
  res <- res[c("mean", "std_error", "std", "quantiles", "loss", "ESS", "Q", "W")]

  class(res) <- "bgphy_posterior"
  attr(res, "model") <- model
  attr(res, "X") <- X
  attr(res, "Xppred") <- ppred$Xppred

  return(res)
}



#' Print posterior information
#'
#' This function prints the information of a posterior distribution on the console.
#'
#' @param x An object of class \code{bgphy_posterior}.
#' @param ... Further argument to be passed to `print`.
#'
#' @examples
#' \dontrun{
#' # global model
#' OU1 <- setModel(tree = lizardTree, regime_names = "Regime1", modeltypes = "OU")
#' post_OU <- bgphy(OU1,
#'                  matrix(XOU[1,], nrow = 1,
#'                         dimnames = list(NULL, colnames(XOU))))
#' print(post_OU)
#'
#' # mixed, OU to OU
#' OUOU <- setModel(tree = lizardTree,
#'                  regime_names = c("Ancestral", "New"),
#'                  modeltypes = c("OU", "OU"),
#'                  startNodes = list(Ancestral = c(101), New = c(135)))
#' post_OUOU <- bgphy(BMOU,
#'                    matrix(XMixedOU[1,], nrow = 1,
#'                           dimnames = list(NULL, colnames(XMixedOU))))
#' print(post_OUOU)
#' }
#'
#' @export
print.bgphy_posterior <- function(x, ...){
  # get parameter names
  par_names <- names(x$mean)

  # get quantity names
  q_names <- c(names(x)[1:3], c("2.5%", "50%", "97.5%"))

  # put information into a matrix
  M <- matrix(unlist(x[1:4]), nrow = length(par_names))
  colnames(M) <- q_names
  rownames(M) <- par_names

  print(M, ...)

  cat("\n")

  cat(paste0("Posterior predictive loss: ", round(x$loss, 5), "\n"))

  cat(paste0("ESS: ", round(x$ESS, 5), "\n"))
}




#' Posterior predictive check
#'
#' Plots a number of density curves drawn from the posterior predictive distribution and the data density curve.
#'
#' @param post An object of class \code{bgphy_posterior}.
#' @param ... Other arguments to be passed on to the `plot` function.
#'
#' @examples
#' \dontrun{
#  # global model
#' OU1 <- setModel(tree = lizardTree, regime_names = "Regime1", modeltypes = "OU")
#' post_OU <- bgphy(OU1,
#'                  matrix(XOU[1,], nrow = 1,
#'                         dimnames = list(NULL, colnames(XOU))))
#' post_pred_check(post_OU)
#'
#' # mixed, OU to OU
#' OUOU <- setModel(tree = lizardTree,
#'                  regime_names = c("Ancestral", "New"),
#'                  modeltypes = c("OU", "OU"),
#'                  startNodes = list(Ancestral = c(101), New = c(135)))
#' post_OUOU <- bgphy(BMOU,
#'                    matrix(XMixedOU[1,], nrow = 1,
#'                           dimnames = list(NULL, colnames(XMixedOU))))
#' post_pred_check(post_OUOU)
#' }
#'
#' @export
post_pred_check <- function(post, ...){
  # check if input is of class bgphy_posterior
  stopifnot("Object is not of class bgphy_posterior." = class(post) == "bgphy_posterior")

  # since bgphy_posterior contains posterior predictive samples,
  # just pick the first min(100, nsamples)

  N <- min(c(100, length(post$W)))

  X <- attr(post, "X")
  Xrep <- attr(post, "Xppred")[1:N,]

  # ====================== plotting routines ======================
  # calculate maximum range for y
  max_y <- max(apply(Xrep, 1, function(x){stats::density(x, n = 10)$y}))


  # base
  plot(stats::density(X[1,]), main = "",
       xlab = "Data at the tips", ylab = "Density", yaxt = "n",
       col = grDevices::rgb(0,0,0,0), ylim = c(0, 1.1*max_y), ...)
  graphics::grid()

  # posterior predictive densities
  for (i in 1:N){
    graphics::lines(stats::density(Xrep[i,]), col = grDevices::rgb(0.2,0.2,1,0.07))
  }

  # data
  graphics::lines(stats::density(X[1,]), lwd = 3, col = grDevices::rgb(0,0,0.5,1))

  graphics::legend("topleft", legend = c("Data", "Simulated data"),
                   col = c(grDevices::rgb(0,0,0.5,1), grDevices::rgb(0.2,0.2,1,0.5)),
                   lty = c(1,1), lwd = c(2,2), seg.len = .75,
                   bty = "n")

}







# simulate from the posterior distribution using Sampling Importance Resampling

#' Simulate posterior samples
#'
#' Generates sample from the posterior distribution using sampling-importance resampling.
#'
#' @param post An object of class \code{bgphy_posterior}.
#' @param nsample Number of samples or draws from the posterior distribution.
#'
#' @examples
#' \dontrun{
#  # global model
#' OU1 <- setModel(tree = lizardTree, regime_names = "Regime1", modeltypes = "OU")
#' post_OU <- bgphy(OU1,
#'                  matrix(XOU[1,], nrow = 1,
#'                         dimnames = list(NULL, colnames(XOU))))
#' post_simulate(post_OU)
#'
#' # mixed, OU to OU
#' OUOU <- setModel(tree = lizardTree,
#'                  regime_names = c("Ancestral", "New"),
#'                  modeltypes = c("OU", "OU"),
#'                  startNodes = list(Ancestral = c(101), New = c(135)))
#' post_OUOU <- bgphy(BMOU,
#'                    matrix(XMixedOU[1,], nrow = 1,
#'                           dimnames = list(NULL, colnames(XMixedOU))))
#' post_simulate(post_OUOU)
#' }
#'
#' @export
post_simulate <- function(post, nsample = 10000){
  # check if input is of class bgphy_posterior
  stopifnot("Object is not of class bgphy_posterior." = class(post) == "bgphy_posterior")
  # check if nsample is valid
  stopifnot("Invalid value of nsample." = is.numeric(nsample) & ((abs(nsample - round(nsample)) < .Machine$double.eps^0.5)) & (nsample > 0))

  # number of particles
  npart <- length(post$W)

  # sample indices by weights
  ind <- sample(1:npart, nsample, replace = TRUE, prob = post$W)

  post_sample <- post$Q[ind,]

  return(post_sample)
}






