#' @export
IS <- function(model, X, tree, priors, initial, nsample, scale = 1, parallel = TRUE){
  #========== Laplace approximation + Importance sampling ==========
  # number of parameters
  d <- length(initial)

  # ================== change of variables routine (to log space)
  ############
  # priors_tr: f_X to f_Y
  # tr: X to Y and vice versa
  # apply change of variables to priors into unbounded space
  priors_tr <- varchange(priors)

  ############
  # define tr, functions to transform parameters
  # tr$f: to unbounded
  # tr$g: from unbounded
  tr <- trfunc(priors_tr)

  # convert initial values into unbounded space
  initial <- tr$f(initial)

  # simplify lupost as a function of parameters
  lu_post <- function(p){bmgpm::lupost(p, model, X, tree, priors_tr, tr)[[1]]}

  # ================== begin the Laplace approximation routine
  # search posterior mode
  optRes <- optim(par = initial, fn = lu_post, method = "BFGS",
                  control = list(fnscale=-1), hessian = TRUE)
  # posterior mode (log space)
  post_mode <- optRes$par
  # approximated covariance, scale to focus on the important area around mode
  appr_cov <- solve(-optRes$hessian)*scale

  # ================== begin the (self-normalized) Importance sampling routine
  # sample from the normal distribution
  q <- mvtnorm::rmvnorm(nsample, mean = post_mode, sigma = appr_cov)

  if (parallel){
    cl <- parallel::makeCluster(parallel::detectCores(),"PSOCK")
    parallel::clusterExport(cl, varlist = c("PCMloglik", "setParams", "loadParams"))
    # log-unnormalized posterior & loglik
    lup_ll <- t(parallel::parApply(cl, q, 1, bmgpm::lupost,
                                   model, X, tree, priors_tr, tr))
    logp <- lup_ll[,"log_u_post"] # log-unnormalized posterior
    loglik <- lup_ll[,"loglik"] # log-likelihood
    # log proposal
    logq <- parallel::parApply(cl, q, 1, mvtnorm::dmvnorm,
                               mean = post_mode, sigma = appr_cov, log = TRUE)
    on.exit(parallel::stopCluster(cl))
  }else{
    # log-unnormalized posterior & loglik
    lup_ll <- t(apply(q, 1, bmgpm::lupost, model, X, tree, priors_tr, tr)) # log-unnormalized posterior
    logp <- lup_ll[,"log_u_post"] # log-unnormalized posterior
    loglik <- lup_ll[,"loglik"] # log-likelihood
    logq <- apply(q, 1, mvtnorm::dmvnorm,
                  mean = post_mode, sigma = appr_cov, log = TRUE) # log proposal
  }

  logW <- logp - logq # log of weights
  W <- exp(logW) # weights

  logWst <- logW - logsumexp(logW) # log of normalized weights
  Wst <- exp(logWst) # normalized weights
  MWst <- matrix(rep(Wst, d), ncol = d) # matrix of normalized weights

  # compute expected values
  Epars_ls <- colSums(q*MWst) # expected value (in log space)
  Epars <-tr$g(Epars_ls) # transform the expected value to the original space

  # compute WAIC
  meanlp <- sum(loglik*Wst) # mean of log predictive dist
  varlp <- sum((Wst^2)*((loglik - meanlp)^2)) # variance of log predictive dist
  lppd <- logsumexp(loglik + logWst) # log pointwise predictive distribution
  WAIC <- -2*(lppd - varlp) # WAIC

  Q <- t(apply(q, 1, tr$g)) # proposed values on the original space

  return(list(E = Epars, Q = Q, W = Wst, WAIC = WAIC))
}