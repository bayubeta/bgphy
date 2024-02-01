# function to run the Importance Sampling method
IS <- function(model, X, nsample, scale = 1, parallel = TRUE){
  #================ Laplace approximation + Importance sampling ================

  # number of parameters
  d <- length(model$priors)

  # ================== change of variables routine (to log space)
  ############
  # priors_tr: f_X to f_Y
  # tr: X to Y and vice versa
  # apply change of variables to priors into unbounded space
  priors_tr <- lapply(model$priors, prior_transform)

  ############
  # define tr, functions to transform parameters
  # tr$f: to unbounded
  # tr$g: from unbounded
  tr <- trfunc(priors_tr)


  repeat{ # in case of singular Hessian, repeat from different starting point
    # find initial position for optim(), by the mean of the prior distribution
    init <- apply(prior_sampler(model$priors)(1000), 2, stats::median) # set initial value as median of samples drawn from priors

    # transform to unbounded space
    initial <- tr$f(init)

    # simplify lupost as a function of parameters
    lu_post <- function(p){lupost(p, model$model, X, model$tree, priors_tr, tr)[[1]]}

    # ================== begin the Laplace approximation routine
    # search posterior mode
    optRes <- stats::optim(par = initial, fn = lu_post, method = "BFGS",
                           control = list(fnscale=-1), hessian = TRUE)
    # posterior mode (log space)
    post_mode <- optRes$par
    # approximated covariance, scale to focus on the important area around mode
    appr_cov <- try(solve(-optRes$hessian)*scale, silent = TRUE)

    if (all(class(appr_cov) != "try-error")){
      break()
    }
  }


  # ================== begin the (self-normalized) Importance sampling routine
  # sample from the normal distribution
  q <- mvtnorm::rmvnorm(nsample, mean = post_mode, sigma = appr_cov)

  if (parallel){
    cl <- parallel::makeCluster(parallel::detectCores(),"PSOCK")
    parallel::clusterExport(cl, varlist = c("PCMloglik", "setParams", "loadParams"), envir = environment())
    # log-unnormalized posterior & loglik
    lup_ll <- t(parallel::parApply(cl, q, 1, lupost,
                                   model$model, X, model$tree, priors_tr, tr)) # log-unnormalized posterior & likelihood
    logp <- lup_ll[,1] # log-unnormalized posterior
    loglik <- lup_ll[,2] # log-likelihood
    # log proposal
    logq <- parallel::parApply(cl, q, 1, mvtnorm::dmvnorm,
                               mean = post_mode, sigma = appr_cov, log = TRUE) # log density of normal
    on.exit(parallel::stopCluster(cl))

  }else{
    # log-unnormalized posterior & loglik
    lup_ll <- t(apply(q, 1, lupost, model$model, X, model$tree, priors_tr, tr)) # log-unnormalized posterior & likelihood
    logp <- lup_ll[,1] # log-unnormalized posterior
    loglik <- lup_ll[,2] # log-likelihood
    logq <- apply(q, 1, mvtnorm::dmvnorm,
                  mean = post_mode, sigma = appr_cov, log = TRUE) # log density of normal
  }

  logW <- logp - logq # log of weights
  W <- exp(logW) # weights

  logWst <- logW - logsumexp(logW) # log of normalized weights
  Wst <- exp(logWst) # normalized weights
  MWst <- matrix(rep(Wst, d), ncol = d) # matrix of normalized weights

  # compute expected values
  Epars_ls <- colSums(q*MWst) # expected value (in log space)
  Epars <- tr$g(Epars_ls) # transform the expected value to the original space

  # compute WAIC
  meanlp <- sum(loglik*Wst) # mean of log predictive dist
  varlp <- sum((Wst^2)*((loglik - meanlp)^2)) # variance of log predictive dist
  lppd <- logsumexp(loglik + logWst) # log pointwise predictive distribution
  WAIC <- -2*(lppd - varlp) # WAIC

  Q <- t(apply(q, 1, tr$g)) # proposed values on the original space

  return(list(Q = Q, W = Wst, WAIC = WAIC))
}


# function to estimate quantiles given sample and weights
est_quantiles <- function(Q, W, probs = c(0.025, 0.5, 0.975)){

  # number of params
  d <- dim(Q)[2]

  # matrix of quantiles
  qM <- matrix(nrow = d, ncol = length(probs))

  for (i in 1:d){
    # sort the parameter values
    qSort <- sort(Q[,i], index.return = TRUE)
    # sorted parameter values
    q_sorted <- qSort$x
    # indices
    idx <- qSort$ix

    # compute quantiles
    qtls <- sapply(probs, function(p){min(which(cumsum(W[idx]) >= p))})

    # store to matrix
    qM[i,] <- q_sorted[qtls]
  }

  colnames(qM) <- probs


  return(qM)
}


