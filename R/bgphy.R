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
#' * `quantiles`:  Estimated quantiles.
#' * `Q`:   Drawn values from the proposal distribution (Multivariate Normal).
#' * `W`:   Normalized weights for each row of `Q`.
#' * `WAIC`:   Widely Available Criterion Score.
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


  # ---------------------------- start inference ----------------------------
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



#' Print posterior information
#'
#' This function prints the information of a posterior distribution on the console.
#'
#' @param post An object of class \code{bgphy_posterior}.
#'
#' @examples
#' \dontrun{
#  # global model
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

  cat(paste0("WAIC: ", post$WAIC, "\n"))
}




#' Posterior predictive check
#'
#' Plots a number of density curves drawn from the posterior predictive distribution and the data density curve.
#'
#' @param post An object of class \code{bgphy_posterior}.
#' @param nsim Number of simulations, or draws from the posterior predictive distribution.
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
post_pred_check <- function(post, nsim = 100){
  # check if input is of class bgphy_posterior
  stopifnot("Object is not of class bgphy_posterior." = class(post) == "bgphy_posterior")
  # check if nsample is valid
  stopifnot("Invalid value of nsample." = is.numeric(nsim) & ((abs(nsim - round(nsim)) < .Machine$double.eps^0.5)) & (nsim > 0))

  N <- length(post$W)

  # retrieve model
  model <- attr(post, "model")

  # retrieve data
  X <- attr(post, "X")

  # number of tips of the tree
  ntips <- ape::Ntip(model$tree)

  # sample according to weights
  post_par <- post$Q[sample(1:N, nsim, prob = post$W, replace = TRUE),]

  # empty matric for storing Xrep (posterior predictive values)
  Xrep <- matrix(nrow = nsim, ncol = ncol(X))

  for (i in 1:nsim){
    # assign posterior params to model
    post_model <- setParams(post_par[i,], model$model)

    # simulate data using PCMBase, take only the tips
    Xrep[i,] <- PCMBase::PCMSim(model$tree, post_model, post_model$X0)[,1:ntips]
  }

  # ====================== plotting routines ======================
  # calculate maximum range for y
  max_y <- max(apply(Xrep, 1, function(x){stats::density(x, n = 10)$y}))


  # base
  graphics::par(mar = c(3, 1, 1, 1))
  plot(stats::density(X[1,]), main = "", xlab = "", ylab = "", yaxt = "n",
       col = grDevices::rgb(0,0,0,0), ylim = c(0, 1.1*max_y))
  graphics::grid()

  # posterior predictive densities
  for (i in 1:nsim){
    graphics::lines(stats::density(Xrep[i,]), col = grDevices::rgb(0.2,0.2,1,0.07))
  }

  # data
  graphics::lines(stats::density(X[1,]), lwd = 3, col = grDevices::rgb(0,0,0.5,1))

  graphics::legend("topleft", legend = c("Data", "Simulated data"),
                   col = c(grDevices::rgb(0,0,0.5,1), grDevices::rgb(0.2,0.2,1,0.5)),
                   lty = c(1,1), lwd = c(2,1), seg.len = .75,
                   bty = "n")

}

